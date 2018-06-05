#include<gmp.h>
#include"flint/flint.h"
#include"flint/nmod_vec.h"
#include"flint/nmod_poly.h"
#include"flint/nmod_mat.h"
#include"flint/ulong_extras.h"

void p_ceil(nmod_poly_t res,nmod_poly_t f, slong h){
	nmod_poly_set(res,f);
	nmod_poly_truncate(res,h);
}
void p_floor(nmod_poly_t res, nmod_poly_t f,slong l){
	nmod_poly_shift_right(res,f,l);
}
void t_mul(nmod_poly_t res, slong n, nmod_poly_t a, nmod_poly_t c){
	slong m = nmod_poly_degree(a);
	if(nmod_poly_degree(c) > m + n){
		printf("ERROR\n");
		nmod_poly_zero(res);
	}
	nmod_poly_reverse(res,a,m+1);
	nmod_poly_mul(res,res,c);
	p_ceil(res,res,n+m+1);
	p_floor(res,res,m);
}
void t_kar_mul(nmod_poly_t res, slong n, nmod_poly_t a, nmod_poly_t c){
	slong m = nmod_poly_degree(a);
	//char * x = "x";
	if(nmod_poly_degree(c) > n + m){
		printf("ERROR\n");
		//nmod_poly_print_pretty(a,x);printf("\n");
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
	mu = m/2 + 1;
	vu = n/2 + 1;
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
		t_kar_mul(u,n,p1,p2);
		p_floor(p1,a,mu);
		p_floor(p2,c,mu);
		t_kar_mul(v,n,p1,p2);
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
		t_kar_mul(u,vu-1,a,p1);
		p_floor(p1,c,vu);
		t_kar_mul(v,n-vu,a,p1);
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
	t_kar_mul(r,h-1,p1,p2);

	p_floor(p1,a,h);
	p_floor(p2,c,2*h);
	p_ceil(p3,c,m+n-h+1);
	p_floor(p3,p3,h);
	nmod_poly_sub(p2,p2,p3);
	t_kar_mul(s,n-h,p1,p2);

	p_ceil(p1,a,h);
	p_floor(p2,a,h);
	nmod_poly_add(p1,p1,p2);
	p_ceil(p2,c,lambda+2*h);
	p_floor(p2,p2,h);
	t_kar_mul(t,h-1,p1,p2);
	
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
	t_mul(G,N,up,wp);
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
