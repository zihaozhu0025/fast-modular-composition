#include"fast_mod_com_precompute_preinv.h"
#include<time.h>
int main(){
	char *x = "x";
	nmod_poly_t f,g,h,r,p,finv;
	nmod_mat_t A,L,A1;
	int T;
	scanf("%d",&T);
	nmod_poly_init(g,9001);//initializes a polynomial array modulo 109
	nmod_poly_init(f,9001);
	nmod_poly_init(h,9001);
	nmod_poly_init(r,9001);
	nmod_poly_init(p,9001);
	nmod_poly_init(finv, 9001);
	nmod_poly_read(g);
	nmod_poly_read(h);
	nmod_poly_read(f);
	nmod_t mod = f->mod;
	slong n = nmod_poly_degree(f);
	slong m = n_sqrt(n)+1;
	slong N = n*2 - 1;
	nmod_mat_init(A,m,N,mod.n);
	nmod_mat_init(L,m,N,mod.n);
	nmod_mat_init(A1,m,n,mod.n);
	long fails=0;
	mp_ptr qs;
	nmod_poly_reverse(finv, f , f->length);
	nmod_poly_inv_series(finv, finv, f->length);
	clock_t c1,c2,d;
	c2=clock();
	nmod_poly_precompute_matrix(A1,h,f,finv);
	c2=clock()-c2;
	c1=clock();
	qs = _nmod_vec_init(N);
	_precompute_geometric_series(qs, mod, N);	
	_nmod_poly_precompute_matrix_interp(A,L,h->coeffs,h->length,f->coeffs,f->length,
			finv->coeffs, finv->length, f->mod, qs);
	c1=clock() - c1;
	for(slong i=1;i<=T;i++){
		d=clock();
		nmod_poly_compose_mod_brent_kung_precomp_preinv(r,g,A1,f,finv);
		d=clock()-d;
		c2+=d;
		d=clock();
		nmod_poly_compose_mod_brent_kung_fast_precomp_preinv(p,g,A,L,f, finv,qs);
		d=clock()-d;
		c1+=d;
		if(!nmod_poly_equal(r,p)){
			printf("failed for\n");
			printf("g ");nmod_poly_print_pretty(g,x);flint_printf("\n\n");
			printf("h ");nmod_poly_print_pretty(h,x);flint_printf("\n\n");
			printf("f ");nmod_poly_print_pretty(f,x);flint_printf("\n\n");
			nmod_poly_print_pretty(r,x);flint_printf("\n\n");
			nmod_poly_print_pretty(p,x);flint_printf("\n\n");
			fails++;
		}
		else{
			printf("passed %li/%d\n",i,T);
		}
		nmod_poly_read(g);
	}
	printf("%li/%d passed \n", T-fails,T);
	double t1,t2;
	t1 = ((double)c1)/CLOCKS_PER_SEC;
	t2 = ((double)c2)/CLOCKS_PER_SEC;
	printf("c1 is %f, c2 is %f\n",t1,t2);
	nmod_poly_clear(p);
	nmod_poly_clear(r);
	nmod_poly_clear(f);
	nmod_poly_clear(g);
	nmod_poly_clear(h);
	return 0;
}
