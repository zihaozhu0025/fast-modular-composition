#include "geometric_interpolation.h"
#include<gmp.h>
#include"flint/flint.h"
#include"flint/nmod_vec.h"
#include"flint/nmod_poly.h"
#include"flint/nmod_mat.h"
#include"flint/ulong_extras.h"

void _precompute_geometric_series(mp_ptr qs, nmod_t mod, slong N){
	slong q;
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
}

void _nmod_poly_precompute_matrix_interp(nmod_mat_t A, nmod_mat_t L, const mp_srcptr poly1, slong len1,
	   	const mp_srcptr poly2, slong len2, const mp_srcptr poly2inv, slong len2inv, nmod_t mod, const mp_ptr qs)
{
	slong i, n, N, m;
	n = len2 - 1;
	N = 2*n-1;
	m = n_sqrt(n) + 1;
	for(i=0;i<N;i++){
		A->rows[0][i]=1;
		//L->rows[0][i]=1;
	}
	L->rows[0]=A->rows[0];
	mp_ptr hi,hm;
	hi = _nmod_vec_init(n);
	hm = _nmod_vec_init(n);
	_nmod_vec_set(hi,poly1,len1);
	for(i=1;i<m;i++){
		_nmod_poly_evaluate_nmod_vec_fast(A->rows[i], hi, n, qs, N, mod);
		_nmod_poly_mulmod_preinv(hi,hi,n,poly1,len1,poly2,len2,poly2inv,len2inv,mod);
	}
	_nmod_vec_set(hm,hi,n);
	for(i=1;i<m;i++){
		_nmod_poly_evaluate_nmod_vec_fast(L->rows[i], hi, n, qs, N, mod);
		_nmod_poly_mulmod_preinv(hi,hi,n,hm,n,poly2,len2,poly2inv, len2inv, mod);
	}

}
void
_nmod_poly_compose_mod_brent_kung_fast_precomp_preinv(mp_ptr res, const mp_srcptr poly1, slong len1,
                            const nmod_mat_t A, const nmod_mat_t L,
                            const mp_srcptr poly3, slong len3, const mp_srcptr poly3inv, const mp_ptr qs, nmod_t mod)
{
    nmod_mat_t B, C;
    mp_ptr t, h,hi,li;
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
        res[0] = _nmod_poly_evaluate_nmod(poly1, len1, A->rows[1][0], mod);
        return;
    }
    m = n_sqrt(n)+1;
	//set N
	N = 2*n-1; 
	nmod_mat_init(B, m, m, mod.n);
    nmod_mat_init(C, m, N, mod.n);

    h = _nmod_vec_init(n);
    t = _nmod_vec_init(n);
    /* Set rows of B to the segments of poly1 */
    for (i = 0; i < len1 / m; i++)
        _nmod_vec_set(B->rows[i], poly1 + i*m, m);

    _nmod_vec_set(B->rows[i], poly1 + i*m, len1 % m);
    nmod_mat_mul(C, B, A);
	mp_ptr r;
	r = _nmod_vec_init(N);
	_nmod_vec_zero(r,N);
	for(i=0;i<m;i++){
		for(j=0;j<N;j++){
			r[j] = nmod_add(r[j],nmod_mul(C->rows[i][j] , L->rows[i][j], mod),mod);
		}
	}
	mp_ptr ng = newton_interp_geom(r,qs,N,mod);
	mp_ptr nm = newton_to_monomial_geom(ng,qs,N,mod);
	//_nmod_vec_set(r,nm,N);
	_nmod_poly_rem(res, nm, N, poly3, len3, mod);
	_nmod_vec_clear(ng);
	_nmod_vec_clear(nm);
	_nmod_vec_clear(r);
    _nmod_vec_clear(h);
    _nmod_vec_clear(t);
    nmod_mat_clear(B);
    nmod_mat_clear(C);
}
void
nmod_poly_compose_mod_brent_kung_fast_precomp_preinv(nmod_poly_t res,
                    const nmod_poly_t poly1, const nmod_mat_t A, const nmod_mat_t L,
                    const nmod_poly_t poly3, const nmod_poly_t poly3inv, const mp_ptr qs)
{
    slong len1 = poly1->length;
    slong len3 = poly3->length;
    slong len = len3 - 1;

    if (len3 == 0)
    {
        flint_printf("Exception (nmod_poly_compose_mod_brent_kung_precomp_preinv). Division by zero.\n");
        flint_abort();
    }

    if (len1 >= len3)
    {
        flint_printf("Exception (nmod_poly_compose_mod_brent_kung_precomp_preinv). The degree of the \n"
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
        nmod_poly_compose_mod_brent_kung_fast_precomp_preinv(tmp, poly1, A, L,
                                                poly3,poly3inv,qs);
        nmod_poly_swap(tmp, res);
        nmod_poly_clear(tmp);
        return;
    }

    nmod_poly_fit_length(res, len);
    _nmod_poly_compose_mod_brent_kung_fast_precomp_preinv(res->coeffs,
        poly1->coeffs, len1, A, L, poly3->coeffs, len3, poly3inv->coeffs, qs, res->mod);
    res->length = len;
    _nmod_poly_normalise(res);

}
