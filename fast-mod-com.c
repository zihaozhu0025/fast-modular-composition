#include <gmp.h>
#include "flint/flint.h"
#include "flint/nmod_vec.h"
#include "flint/nmod_poly.h"
#include "flint/nmod_mat.h"
#include "flint/ulong_extras.h"

void
_nmod_poly_compose_mod_brent_kung_fast(mp_ptr res, mp_srcptr poly1, slong len1, 
                            mp_srcptr poly2,
                            mp_srcptr poly3, slong len3, nmod_t mod)
{
    nmod_mat_t A, B, C;
    mp_ptr t, h;
    slong i, n, m;

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

    m = n_sqrt(n) + 1;

    nmod_mat_init(A, m, n, mod.n);
    nmod_mat_init(B, m, m, mod.n);
    nmod_mat_init(C, m, n, mod.n);

    h = _nmod_vec_init(n);
    t = _nmod_vec_init(n);

    /* Set rows of B to the segments of poly1 */
    for (i = 0; i < len1 / m; i++)
        _nmod_vec_set(B->rows[i], poly1 + i*m, m);

    _nmod_vec_set(B->rows[i], poly1 + i*m, len1 % m);

    /* Set rows of A to powers of poly2 */
    A->rows[0][0] = UWORD(1);
    _nmod_vec_set(A->rows[1], poly2, n);
    for (i = 2; i < m; i++)
        _nmod_poly_mulmod(A->rows[i], A->rows[i-1],
            n, poly2, n, poly3, len3, mod);

    nmod_mat_mul(C, B, A);

    /* Evaluate block composition using the Horner scheme */
    _nmod_vec_set(res, C->rows[m - 1], n);
    _nmod_poly_mulmod(h, A->rows[m - 1], n, poly2, n, poly3, len3, mod);

    for (i = m - 2; i >= 0; i--)
    {
        _nmod_poly_mulmod(t, res, n, h, n, poly3, len3, mod);
        _nmod_poly_add(res, t, n, C->rows[i], n, mod);
    }

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
    _nmod_poly_compose_mod_brent_kung(res->coeffs,
        poly1->coeffs, len1, ptr2, poly3->coeffs, len3, res->mod);
    res->length = len;
    _nmod_poly_normalise(res);

    _nmod_vec_clear(ptr2);
}
int main(){
	nmod_poly_t f,g,h,r;
	nmod_poly_init(g,29);//initializes a polynomial array modulo 29
	nmod_poly_init(f,29);
	nmod_poly_init(h,29);
	nmod_poly_init(r,29);
	nmod_poly_set_coeff_ui(h,4,6);
	nmod_poly_set_coeff_ui(h,2,3);
	nmod_poly_set_coeff_ui(g,3,5);
	nmod_poly_set_coeff_ui(g,0,6);
	nmod_poly_mul(f,g,g);
	//sets h = 6x^4 + 3x^2, g = 5x^3 + 6, f = g^2 = 25x^6+60x^3+36
	nmod_poly_compose_mod_brent_kung_fast(r,g,h,f);
	nmod_poly_print(r); flint_printf("\n");
	nmod_poly_clear(r);
	nmod_poly_clear(f);
	nmod_poly_clear(g);
	nmod_poly_clear(h);
	return 0;
}
