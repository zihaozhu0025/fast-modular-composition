#include<gmp.h>
#include"flint/flint.h"
#include"flint/nmod_vec.h"
#include"flint/nmod_poly.h"
#include"flint/nmod_mat.h"
#include"flint/ulong_extras.h"

void _precompute_geometric_series(mp_ptr, nmod_t, slong);
void _nmod_poly_precompute_matrix_interp(nmod_mat_t, nmod_mat_t, const mp_srcptr, slong, const mp_srcptr, slong, const mp_srcptr, slong, nmod_t, const mp_ptr);
void nmod_poly_compose_mod_brent_kung_fast_precomp_preinv(nmod_poly_t, const nmod_poly_t, const nmod_mat_t, const nmod_mat_t,
		const nmod_poly_t, const nmod_poly_t, const mp_ptr qs);
