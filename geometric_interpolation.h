#include<gmp.h>
#include"flint/flint.h"
#include"flint/nmod_vec.h"
#include"flint/nmod_poly.h"
#include"flint/nmod_mat.h"
#include"flint/ulong_extras.h"
void p_ceil(nmod_poly_t, nmod_poly_t, slong);
void p_floor(nmod_poly_t, nmod_poly_t, slong);
void t_mul(nmod_poly_t, slong, nmod_poly_t, nmod_poly_t);
mp_ptr newton_to_monomial_geom(mp_ptr, mp_ptr, slong, nmod_t);
mp_ptr newton_interp_geom(mp_ptr, mp_ptr, slong, nmod_t);

