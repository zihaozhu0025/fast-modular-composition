#include <gmp.h>
#include<math.h>
#include "flint/flint.h"
#include "flint/nmod_vec.h"
#include "flint/nmod_poly.h"
#include "flint/nmod_mat.h"
#include "flint/ulong_extras.h"
mp_ptr poly_ceil(mp_ptr f, slong len, slong h){
//	mp_limb_t m = nmod_poly_modulus(f);
	mp_ptr res;
//	nmod_poly_init(res, m);
//	nmod_poly_set(res,f);
//	nmod_poly_truncate(res,h);
	slong rlen = FLINT_MIN(len,h);
	res=_nmod_vec_init(rlen);
	_nmod_vec_set(res,f,rlen);
	return res;
}
mp_ptr poly_floor(mp_ptr f, slong len,slong l){
	//mp_limb_t m = nmod_poly_modulus(f);
	mp_ptr res;
	//nmod_poly_init(res,m);
	res = _nmod_vec_init(len);
	if(l>=len){
		_nmod_vec_zero(res,len);
		return res;
	}
	_nmod_poly_shift_right(res,f,len,l);	
	return res;
}

mp_ptr transpose_kar_mul(slong n,mp_ptr a, slong len1, mp_ptr c, slong len2, nmod_t mod){
	slong m = len1-1;//degree of a
	//printf("n is %li m is %li deg(c) is %li\n",n,m, len2-1);
	//printf("c is");for(slong i=len2-1;i>=0;i--){printf(" %lix^%li",c[i],i);}printf("\n");
	//printf("a is");for(slong i=len1-1;i>=0;i--){printf(" %lix^%li",a[i],i);}printf("\n");
	if (FLINT_MAX(len2-1,0) > m + n){
		flint_abort();
	}	
	mp_ptr r,s,t;
	if (n==0){
		//printf("base1\n");
		r = _nmod_vec_init(1);
		slong k;
		slong coeff = 0;
		for (k=0;k<=m;k++){
			coeff = nmod_add(coeff, nmod_mul(a[k], c[k],mod),mod);
		}
		r[0] = coeff;
		//nmod_poly_set_coeff_ui(r, 0, coeff);
		return r;
	}
	if (m<=0){
		//printf("base2\n");
		r = _nmod_vec_init(len2);
		//_nmod_poly_mul(r,a,len1, c,len2,mod);
		return r;
	}
	slong mu,vu,h,rho;
	mu = m/2+1;
	vu = n/2+1;
	h = FLINT_MAX(mu,vu);
	rho = FLINT_MIN(h,len1)-1; // degree(ceil(a,h))
	//printf("rho is %li \n", rho);
	//rho = nmod_poly_degree(poly_ceil(a,h));
	mp_ptr u,v;
	slong lu,lv;
	if (mu > n){
		//printf("secret 1\n");
		lu = FLINT_MIN(len1,mu) + FLINT_MIN(len2,n + rho);
		lv = FLINT_MIN(len2,mu) + FLINT_MAX(0,len2-mu);
		u = transpose_kar_mul(n,poly_ceil(a,len1,mu), FLINT_MIN(len1,mu), poly_ceil(c,len2,n+rho),FLINT_MIN(len2,n+rho),mod);
		v = transpose_kar_mul(n,poly_ceil(a,len2,mu),FLINT_MIN(len2,mu),poly_floor(c,len2,mu),FLINT_MAX(0,len2-mu),mod);
		_nmod_poly_add(r,u,lu,v,lv,mod);
		return r;
	}
	if (vu > m){
		//printf("secret 2\n");
		lu = len1 + FLINT_MIN(m+vu,len2)-1;
		lv = len1 + FLINT_MIN(vu,len2)-1;
		u = transpose_kar_mul(vu-1,a,len1,poly_ceil(c,len2,m+vu),FLINT_MIN(m+vu,len2),mod);
		v = transpose_kar_mul(n-vu,a,len1,poly_ceil(c,len2,vu),FLINT_MIN(vu,len2),mod);
		_nmod_poly_add(r,u,lu,v,lv,mod);
		return r;
	}
	slong lambda;
	//lambda = degree(nmod_poly_add(poly_ceil(a,h),poly_floor(a,h)));
	lambda = FLINT_MAX(h,len1-h)-1;
	//printf("lambda is %li \n",lambda);
	mp_ptr p1,p2,p3,p4,p5,p6;
	slong l1,l2,l3,l4,l5,l6,lr,ls,lt;
	p1 = _nmod_vec_init(n);
	p2 = _nmod_vec_init(n);
	p3 = _nmod_vec_init(n);
	p4 = _nmod_vec_init(n);
	p5 = _nmod_vec_init(n);
	p6 = _nmod_vec_init(n);
	_nmod_poly_sub(p1,poly_ceil(c,len2,rho+h),FLINT_MIN(len2,rho+h),poly_floor(poly_ceil(c,len2,rho+2*h),len2,h),FLINT_MAX(0,FLINT_MIN(len2,rho+2*h)-h),mod);
	l1=FLINT_MIN(len2,rho+h);
	while(l1 && p1[l1-1] == WORD(0)){l1--;}
	_nmod_poly_sub(p2,poly_floor(c,len2,2*h),FLINT_MAX(0,len2-2*h),poly_floor(poly_ceil(c,len2,m+n-h+1),FLINT_MIN(len2,m+n-h+1),h),FLINT_MAX(0,FLINT_MIN(len2,m+n-h+1)-h),mod);
	l2=FLINT_MAX(0,FLINT_MIN(len2,rho+2*h)-h);
	while(l2 && p2[l2-1] == WORD(0)){l2--;}
	_nmod_poly_add(p3,poly_ceil(a,len1,h),FLINT_MIN(len1,h),poly_floor(a,len1,h),FLINT_MAX(0,len1-h),mod);
	l3=FLINT_MAX(FLINT_MIN(len1,h),FLINT_MAX(0,len1-h));
	while(l3 && p3[l3-1] == WORD(0)){l3--;}
	r = transpose_kar_mul(h-1, poly_ceil(a,len1,h),FLINT_MIN(len1,h),p1,l1,mod);
	s = transpose_kar_mul(n-h, poly_floor(a,len1,h), FLINT_MAX(len1-h,0),p2,l2,mod);
	t = transpose_kar_mul(h-1,p3,l3,poly_floor(poly_ceil(c,len2,lambda + 2*h),FLINT_MIN(len2,lambda+2*h),h),FLINT_MAX(0,FLINT_MIN(len2,lambda+2*h)-h),mod);
	lr = FLINT_MIN(len1,h) + l1-1;
	ls = FLINT_MAX(len1-h,0)+l2;
	lt = l3 + FLINT_MIN(len2,lambda+2*h);
	_nmod_poly_add(p4, s, ls, poly_ceil(t,n,n-h+1),FLINT_MIN(lt,n-h+1),mod);
	_nmod_poly_shift_left(p6, p4,FLINT_MAX(ls, FLINT_MIN(lt,n-h+1)),h);
	_nmod_poly_add(p5,t,lt,p6,FLINT_MAX(ls,FLINT_MIN(lt,n-h+1))+h,mod);
	_nmod_poly_add(r,r,lr,p5,FLINT_MAX(ls,FLINT_MIN(lt,n-h+1))+h,mod);
	printf("p4");for(slong i=FLINT_MAX(ls,FLINT_MIN(lt,n-h+1))-1;i>=0;i--){printf(" %lix^%li",p4[i],i);}printf("\n");
	printf("p6");for(slong i=FLINT_MAX(ls,FLINT_MAX(ls,FLINT_MIN(lt,n-h+1)))+h-1;i>=0;i--){printf(" %lix^%li",p6[i],i);}printf("\n");
	printf("p5");for(slong i=FLINT_MAX(lt,FLINT_MAX(ls,FLINT_MIN(lt,n-h+1)))-1;i>=0;i--){printf(" %lix^%li",p5[i],i);}printf("\n");
	return r;
}
mp_ptr newton_to_monomial_geom(mp_ptr f, mp_ptr qs, slong N,nmod_t mod){
	mp_ptr q,v,u,w,z;
	q = _nmod_vec_init(N);
	v = _nmod_vec_init(N);
	u = _nmod_vec_init(N);
	w = _nmod_vec_init(N);
	z = _nmod_vec_init(N);
	q[0]=1;
	v[0]=f[0];
	u[0]=1;
	w[0]=f[0];
	z[0]=1;
	slong i,lenu,lenw;
	for (i=1;i<N;i++){
		q[i] = nmod_mul(q[i-1],qs[i-1],mod);
		v[i] = nmod_mul((i%2==0?1:mod.n-1),nmod_mul(f[i],q[i],mod),mod);
		//u[i] = u[i-1]*qs[i]/(1-qs[i]);
		u[i] = nmod_mul(u[i-1],nmod_div(qs[i],((1-qs[i])+(mod.n))%(mod.n),mod),mod);
		//w[i] = v[i]/u[i];
		w[i] = nmod_div(v[i],u[i],mod);
		z[i] = nmod_mul((i%2==0?1:mod.n-1),nmod_div(u[i],q[i],mod),mod);
		if(u[i]!=0){
			lenu=i+1;
		}
		if(w[i]!=0){
			lenw=i+1;
		}
		//z[i] = (i%2==0?1:-1)*u[i]/q[i];
	}
	mp_ptr G = transpose_kar_mul(N,u,lenu,w,lenw,mod);
	printf("G is");
	for(i=0;i<N;i++){printf(" %lix^%li",G[i],i);}printf("\n");
	for (i=0;i<N;i++){
		G[i]=nmod_mul(G[i],z[i],mod);
	}
	return G;
}
mp_ptr newton_interp_geom(mp_ptr v, mp_ptr qs, slong N,nmod_t mod){
	mp_ptr u,w,q,r;
	u=_nmod_vec_init(N);
	w=_nmod_vec_init(N);
	q=_nmod_vec_init(N);
	r = _nmod_vec_init(N);
	u[0]=1;
	w[0]=v[0];
	q[0]=1;
	slong i;
	for (i=1;i<N;i++){
		q[i]=nmod_mul(q[i-1],qs[i-1],mod);
		u[i] = nmod_mul(u[i-1],qs[i]-1,mod);
		w[i] = nmod_div(v[i],u[i],mod); 
	}
	nmod_poly_t f,g,G;
	char *x = "x";
	nmod_poly_init(f,mod.n);
	nmod_poly_init(g,mod.n);
	nmod_poly_init(G,mod.n);
	for (i=0;i<N;i++){
		nmod_poly_set_coeff_ui(f,i,w[i]);
		nmod_poly_set_coeff_ui(g,i,nmod_mul((i%2==0?1:mod.n-1),nmod_div(q[i],u[i],mod),mod));
	}
	nmod_poly_mul(G,f,g);
	nmod_poly_print_pretty(G,x);printf("\n");
	for(i=0;i<N;i++){
		r[i]=nmod_div(nmod_poly_get_coeff_ui(G,i),q[i],mod);
		flint_printf("%li ",r[i]);
	}
	printf("\n");
	//nmod_poly_print(r); flint_printf("\n");
	return r;
}
void
_nmod_poly_compose_mod_brent_kung_fast(mp_ptr res, mp_srcptr poly1, slong len1, 
                            mp_srcptr poly2,
                            mp_srcptr poly3, slong len3, nmod_t mod)
{

    nmod_mat_t A, B, C, L;
    mp_ptr t, h, qs,hi,li,coeffs;
    slong i, j,n, m, N, q;

    n = len3 - 1;

    if (len3 == 1)
        return;

    if (len1 == 1)
    {
        res[0] = poly1[0];
        return;
    }

    if (len3 == 2)
    {
        res[0] = _nmod_poly_evaluate_nmod(poly1, len1, poly2[0], mod);
        return;
    }
	N = 2*n;
	qs = _nmod_vec_init(N);
	qs[0]=1;
	for(q = 2; q < mod.n; q++){
		slong k = 1;
		slong qk = q;
		while(qk != 1 && k < N){
			qs[k]=qk;
			qk = ((qk*q)% mod.n);
			k++;
		}
		if(qk!=1){
			break;
		}
	}
	printf("q is %li \n", q);
	if(q==N){
		flint_abort();
	}
    m = ceil(n_sqrt(n));

	nmod_mat_init(L, m, N, mod.n);
    nmod_mat_init(A, m, N, mod.n);
    nmod_mat_init(B, m, m, mod.n);
    nmod_mat_init(C, m, N, mod.n);

    h = _nmod_vec_init(n);
    t = _nmod_vec_init(n);
	li = _nmod_vec_init(n);
	hi = _nmod_vec_init(n);
    /* Set rows of B to the segments of poly1 */
    for (i = 0; i < len1 / m; i++)
        _nmod_vec_set(B->rows[i], poly1 + i*m, m);

    _nmod_vec_set(B->rows[i], poly1 + i*m, len1 % m);
	char * x = "x";
    /* Set rows of A to powers of poly2 */
	for(j=0;j<N;j++){
    	A->rows[0][j] = UWORD(1);
	}
	//A->rows[0][0]=UWORD(1);
	_nmod_vec_set(hi, poly2, n);
    //_nmod_vec_set(A->rows[1], poly2, n);
    for (i = 1; i < m; i++){
		_nmod_poly_evaluate_nmod_vec_iter(A->rows[i], hi,n, qs, N,mod);
		_nmod_poly_mulmod(hi,hi, n, poly2, n, poly3, len3,mod);
        //_nmod_poly_mulmod(A->rows[i], A->rows[i-1],
        //    n, poly2, n, poly3, len3, mod);
    }
    nmod_mat_mul(C, B, A);

    /* Evaluate block composition using the Horner scheme */
    _nmod_vec_set(res, C->rows[m - 1], n);
	//_nmod_poly_mulmod(h, A->rows[m - 1], n, poly2, n, poly3, len3, mod);
	for(j=0;j<N;j++)
		L->rows[0][j] = UWORD(1);
	_nmod_vec_set(li, hi, n);
	for (i=1;i<m;i++){
		_nmod_poly_evaluate_nmod_vec_iter(L->rows[i], li, n, qs, N, mod);
		_nmod_poly_mulmod(li,li,n,hi,n,poly3,len3,mod);
	}
	mp_ptr r;
	r = _nmod_vec_init(N);
	_nmod_vec_zero(r,N);
	for(i=0;i<m;i++){
		for(j=0;j<N;j++){
			//r[j] += C->rows[i][j] * L->rows[i][j];
			r[j] = nmod_add(r[j],nmod_mul(C->rows[i][j] , L->rows[i][j], mod),mod);
		}
	}
	coeffs =  newton_to_monomial_geom(newton_interp_geom(r,qs,N,mod),qs,N,mod);
    //for (i = m - 2; i >= 0; i--)
    //{
    //    _nmod_poly_mulmod(t, res, n, h, n, poly3, len3, mod);
    //    _nmod_poly_add(res, t, n, C->rows[i], n, mod);
    //}

    _nmod_vec_clear(h);
    _nmod_vec_clear(t);

    nmod_mat_clear(A);
    nmod_mat_clear(B);
    nmod_mat_clear(C);
}

void
nmod_poly_compose_mod_brent_kung_fast(nmod_poly_t res, 
                    const nmod_poly_t poly1, const nmod_poly_t poly2,
                    const nmod_poly_t poly3)
{
    slong len1 = poly1->length;
    slong len2 = poly2->length;
    slong len3 = poly3->length;
    slong len = len3 - 1;

    mp_ptr ptr2;

    if (len3 == 0)
    {
        flint_printf("Exception (nmod_poly_compose_mod_brent_kung). Division by zero.\n");
        flint_abort();
    }

    if (len1 >= len3)
    {
        flint_printf("Exception (nmod_poly_compose_brent_kung). The degree of the \n"
               "first polynomial must be smaller than that of the modulus.\n");
        flint_abort();
    }

    if (len1 == 0 || len3 == 1)
    {
        nmod_poly_zero(res);
        return;
    }

    if (len1 == 1)
    {
        nmod_poly_set(res, poly1);
        return;
    }

    if (res == poly3 || res == poly1)
    {
        nmod_poly_t tmp;
        nmod_poly_init_preinv(tmp, res->mod.n, res->mod.ninv);
        nmod_poly_compose_mod_brent_kung(tmp, poly1, poly2, poly3);
        nmod_poly_swap(tmp, res);
        nmod_poly_clear(tmp);
        return;
    }

    ptr2 = _nmod_vec_init(len);

    if (len2 <= len)
    {
        flint_mpn_copyi(ptr2, poly2->coeffs, len2);
        flint_mpn_zero(ptr2 + len2, len - len2);
    }
    else
    {
        _nmod_poly_rem(ptr2, poly2->coeffs, len2,
                             poly3->coeffs, len3, res->mod);
    }

    nmod_poly_fit_length(res, len);
    _nmod_poly_compose_mod_brent_kung_fast(res->coeffs,
        poly1->coeffs, len1, ptr2, poly3->coeffs, len3, res->mod);
    res->length = len;
    _nmod_poly_normalise(res);

    _nmod_vec_clear(ptr2);
}
int main(){
	char *x = "x";
	nmod_poly_t f,g,h,r;
	nmod_poly_init(g,109);//initializes a polynomial array modulo 109
	nmod_poly_init(f,109);
	nmod_poly_init(h,109);
	nmod_poly_init(r,109);
	nmod_poly_read(g);
	nmod_poly_read(h);
	nmod_poly_read(f);
	nmod_poly_print_pretty(g,x);flint_printf("\n");
	//mp_ptr gt = poly_floor(poly_ceil(g->coeffs,g->length, 14),g->length,5);
	//flint_printf("ceiled");for(slong i=0;i<5;i++){printf(" %lix^%li",gt[i],i);}printf("\n");

	nmod_poly_print_pretty(h,x);flint_printf("\n");
	nmod_poly_print_pretty(f,x);flint_printf("\n");
	nmod_poly_compose_mod_brent_kung_fast(r,g,h,f);
	nmod_poly_print_pretty(r,x); flint_printf("\n");
	//nmod_poly_truncate(r,4);
	//nmod_poly_print_pretty(r,x); flint_printf("\n");
	nmod_poly_clear(r);
	nmod_poly_clear(f);
	nmod_poly_clear(g);
	nmod_poly_clear(h);
	return 0;
}
