#include <gmp.h>
#include "flint/flint.h"
#include "flint/nmod_vec.h"
#include "flint/nmod_poly.h"
#include "flint/nmod_mat.h"
#include "flint/ulong_extras.h"
void p_ceil(nmod_poly_t res,nmod_poly_t f, slong h){
	nmod_poly_set(res,f);
	nmod_poly_truncate(res,h);
}
void p_floor(nmod_poly_t res, nmod_poly_t f,slong l){
	nmod_poly_shift_right(res,f,l);
}
void t_kar_mul(nmod_poly_t res, slong n, nmod_poly_t a, nmod_poly_t c,slong acc){
	slong m = nmod_poly_degree(a);
	if(nmod_poly_degree(c) > n + m){
		printf("ERROR after %li calls\n",acc);
		//flint_abort();
		nmod_poly_zero(res);	
	}
	nmod_t mod = a->mod; 
	nmod_poly_t r,s,t;
	if(n==0){
		slong k,coeff;
		coeff = 0;
		for(k=0;k<=m;k++){
			coeff = nmod_add(coeff, nmod_mul(nmod_poly_get_coeff_ui(a,k), nmod_poly_get_coeff_ui(c,k),mod),mod);
		}
		nmod_poly_set_coeff_ui(res,0,coeff);
		return;
	}
	if(m <=0){
		nmod_poly_scalar_mul_nmod(res, c , nmod_poly_get_coeff_ui(a,0));
		return;
	}
	slong mu,vu,h,rho;
	mu = m/2+1;
	vu = n/2+1;
	h = FLINT_MAX(mu,vu);
	rho = FLINT_MIN(h,a->length) - 1;
	nmod_poly_t u,v;
	nmod_poly_init(u,mod.n);
	nmod_poly_init(v,mod.n);
	nmod_poly_t p1,p2,p3;
	nmod_poly_init(p1,mod.n);
	nmod_poly_init(p2,mod.n);
	nmod_poly_init(p3,mod.n);
	if(mu > n){
		p_ceil(p1,a,mu);
		p_ceil(p2,c,n+rho);
		t_kar_mul(u,n,p1,p2,acc+1);
		p_floor(p1,a,mu);
		p_floor(p2,c,mu);
		t_kar_mul(v,n,p1,p2,acc+1);
		nmod_poly_add(res,u,v);
		nmod_poly_clear(p1);
		nmod_poly_clear(p2);
		nmod_poly_clear(p3);
		nmod_poly_clear(u);
		nmod_poly_clear(v);
		return;
	}
	if(vu > m){
		p_ceil(p1,c,m+vu);
		t_kar_mul(u,vu-1,a,p1,acc+1);
		p_floor(p1,c,vu);
		t_kar_mul(v,n-vu,a,p1,acc+1);
		nmod_poly_add(res,u,v);
		nmod_poly_clear(p1);
		nmod_poly_clear(p2);
		nmod_poly_clear(p3);
		nmod_poly_clear(u);
		nmod_poly_clear(v);
		return;
	}
	nmod_poly_init(r,mod.n);
	nmod_poly_init(s,mod.n);
	nmod_poly_init(t,mod.n);
	p_ceil(p1,a,h);
	p_floor(p2,a,h);
	nmod_poly_add(p1,p1,p2);
	slong lambda = nmod_poly_degree(p1);
	p_ceil(p1,a,h);
	p_ceil(p2,c,rho+h);
	p_ceil(p3,c,rho+2*h);
	p_floor(p3,p3,h);
	nmod_poly_sub(p2,p2,p3);
	t_kar_mul(r,h-1,p1,p2,acc+1);

	p_floor(p1,a,h);
	p_floor(p2,c,2*h);
	p_ceil(p3,c,m+n-h+1);
	p_floor(p3,p3,h);
	nmod_poly_sub(p2,p2,p3);
	t_kar_mul(s,n-h,p1,p2,acc+1);

	p_ceil(p1,a,h);
	p_floor(p2,a,h);
	nmod_poly_add(p1,p1,p2);
	p_ceil(p2,c,lambda+2*h);
	p_floor(p2,p2,h);
	t_kar_mul(t,h-1,p1,p2,acc+1);
	
	p_ceil(res,t,n-h+1);
	nmod_poly_add(res,s,res);
	nmod_poly_shift_left(res,res,h);
	nmod_poly_add(res,t,res);
	nmod_poly_add(res,r,res);
	nmod_poly_clear(p1);
	nmod_poly_clear(p2);
	nmod_poly_clear(p3);
	nmod_poly_clear(r);
	nmod_poly_clear(s);
	nmod_poly_clear(t);
	nmod_poly_clear(u);
	nmod_poly_clear(v);
	return;
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
	nmod_poly_t up,wp;
	nmod_poly_init(up,mod.n);
	nmod_poly_init(wp,mod.n);
	nmod_poly_set_coeff_ui(up,0,u[0]);
	nmod_poly_set_coeff_ui(wp,0,w[0]);
	for (i=1;i<N;i++){
		q[i] = nmod_mul(q[i-1],qs[i-1],mod);
		v[i] = nmod_mul((i%2==0?1:mod.n-1),nmod_mul(f[i],q[i],mod),mod);
		u[i] = nmod_mul(u[i-1],nmod_div(qs[i],((1-qs[i])+(mod.n))%(mod.n),mod),mod);
		w[i] = nmod_div(v[i],u[i],mod);
		z[i] = nmod_mul((i%2==0?1:mod.n-1),nmod_div(u[i],q[i],mod),mod);
		if(u[i]!=0){
			lenu=i+1;
		}
		if(w[i]!=0){
			lenw=i+1;
		}
		nmod_poly_set_coeff_ui(up,i,u[i]);
		nmod_poly_set_coeff_ui(wp,i,w[i]);
	}
	nmod_poly_t G;
	nmod_poly_init(G, mod.n);
	t_kar_mul(G,N,up,wp,0);
	char * x = "x";
	mp_ptr Gp;
	Gp = _nmod_vec_init(N);
	_nmod_vec_zero(Gp,N);
	for (i=0;i<=nmod_poly_degree(G);i++){
		//G[i]=nmod_mul(G[i],z[i],mod);
		Gp[i]=nmod_mul(nmod_poly_get_coeff_ui(G,i),z[i],mod);
	}
	nmod_poly_clear(G);
	_nmod_vec_clear(q);
	_nmod_vec_clear(v);
	_nmod_vec_clear(u);
	_nmod_vec_clear(w);
	_nmod_vec_clear(z);
	nmod_poly_clear(up);
	nmod_poly_clear(wp);
	return Gp;
}
mp_ptr newton_interp_geom(mp_ptr v, mp_ptr qs, slong N,nmod_t mod){
	mp_ptr u,w,q,r;
	u=_nmod_vec_init(N);
	w=_nmod_vec_init(N);
	q=_nmod_vec_init(N);
	r=_nmod_vec_init(N);
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
	for(i=0;i<N;i++){
		r[i]=nmod_div(nmod_poly_get_coeff_ui(G,i),q[i],mod);
	}
	_nmod_vec_clear(u);
	_nmod_vec_clear(w);
	_nmod_vec_clear(q);
	nmod_poly_clear(f);
	nmod_poly_clear(g);
	nmod_poly_clear(G);
	return r;
}
void
_nmod_poly_compose_mod_brent_kung_fast(mp_ptr res, mp_srcptr poly1, slong len1, 
                            mp_srcptr poly2,
                            mp_srcptr poly3, slong len3, nmod_t mod)
{
    nmod_mat_t A, B, C, L;
    mp_ptr t, h, qs,hi,li;
    slong i,j,n, m,N,q;

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
    m = n_sqrt(n)+1;
	//set N
	N = 2*n-1;
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
	if(q==N){
		flint_abort();
	}
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
    /* Set rows of A to powers of poly2 */
	for(j=0;j<N;j++){
    	A->rows[0][j] = 1;
	}
	//A->rows[0][0]=UWORD(1);
	_nmod_vec_set(hi, poly2, n);
    //_nmod_vec_set(A->rows[1], poly2, n);
    for (i = 1; i < m; i++){
		_nmod_poly_evaluate_nmod_vec_iter(A->rows[i], hi,n, qs, N,mod);
		_nmod_poly_mulmod(hi,hi, n, poly2, n, poly3, len3,mod);
    }

    nmod_mat_mul(C, B, A);

    /* Evaluate block composition using the Horner scheme */
    _nmod_vec_set(res, C->rows[m - 1], n);
	for(j=0;j<N;j++)
		L->rows[0][j] = 1;
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
			r[j] = nmod_add(r[j],nmod_mul(C->rows[i][j] , L->rows[i][j], mod),mod);
		}
	}
	//mp_ptr tmp;	
	mp_ptr ng = newton_interp_geom(r,qs,N,mod);
	mp_ptr nm = newton_to_monomial_geom(ng,qs,N,mod);
	_nmod_vec_set(r,nm,N);
	mp_ptr blank;
	blank = _nmod_vec_init(N);
	_nmod_poly_divrem(blank, res,r,N,poly3,len3,mod);
	_nmod_vec_clear(ng);
	_nmod_vec_clear(nm);
	_nmod_vec_clear(blank);
	_nmod_vec_clear(r);
    _nmod_vec_clear(h);
    _nmod_vec_clear(t);
	_nmod_vec_clear(qs);
	_nmod_vec_clear(li);
	_nmod_vec_clear(hi);

    nmod_mat_clear(A);
    nmod_mat_clear(B);
    nmod_mat_clear(C);
	nmod_mat_clear(L);
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
	nmod_poly_t f,g,h,r,p;
	nmod_poly_init(g,9001);//initializes a polynomial array modulo 109
	nmod_poly_init(f,9001);
	nmod_poly_init(h,9001);
	nmod_poly_init(r,9001);
	nmod_poly_init(p,9001);
	flint_rand_t state;
	//flint_randinit(state);
	//nmod_poly_randtest_monic(g,state,i-1);
	//nmod_poly_randtest_monic(h,state,i-1);
	//nmod_poly_randtest_monic(f,state,i);
	nmod_poly_read(g);
	nmod_poly_read(h);
	nmod_poly_read(f);
	printf("g ");nmod_poly_print_pretty(g,x);flint_printf("\n");
	printf("h ");nmod_poly_print_pretty(h,x);flint_printf("\n");
	printf("f ");nmod_poly_print_pretty(f,x);flint_printf("\n");
	long fails=0;
	int T = 100;
	for(slong i=0;i<T;i++){
		printf("n is %li \n",nmod_poly_degree(f));
		nmod_poly_set_coeff_ui(g,nmod_poly_degree(g)+1,i);
		nmod_poly_set_coeff_ui(h,nmod_poly_degree(h)+1,i);
		nmod_poly_set_coeff_ui(f,nmod_poly_degree(f)+1,i);
		nmod_poly_compose_mod_brent_kung(r,g,h,f);
		//nmod_poly_print_pretty(r,x); flint_printf("\n");
		nmod_poly_compose_mod_brent_kung_fast(p,g,h,f);
		//printf("result is ");nmod_poly_print_pretty(p,x); flint_printf("\n");
		if(!nmod_poly_equal(r,p)){
			printf("failed for\n");
			printf("g ");nmod_poly_print_pretty(g,x);flint_printf("\n\n");
			printf("h ");nmod_poly_print_pretty(h,x);flint_printf("\n\n");
			printf("f ");nmod_poly_print_pretty(f,x);flint_printf("\n\n");
			nmod_poly_print_pretty(r,x);flint_printf("\n\n");
			nmod_poly_print_pretty(p,x);flint_printf("\n\n");
			fails++;
		}
	}
	printf("%li/%d passed \n", T-fails,T);
	nmod_poly_clear(p);
	nmod_poly_clear(r);
	nmod_poly_clear(f);
	nmod_poly_clear(g);
	nmod_poly_clear(h);
	return 0;
}
