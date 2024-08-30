#include <NTL/ZZ_pEX.h>
#include <NTL/ZZ_pEXFactoring.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZ_pXFactoring.h>
#include <cassert>
#include <optional>
#include "Fp2k.hpp"

size_t myroots(NTL::ZZ_pE *roots, NTL::ZZ_pEX const &f)
{
    //XXX This is only marginally faster than NTL::CanZass().
    //XXX Is it worth the extra (potentially buggy) code?
    assert(NTL::IsOne(NTL::LeadCoeff(f)));
    NTL::ZZ_pEX x; NTL::SetX(x);
    auto sfd = NTL::SquareFreeDecomp(f);
    size_t idx = 0;
    for (auto const &[g,m]: sfd) {
        NTL::ZZ_pEXModulus F;
        NTL::build(F, g);
        auto h = NTL::PowerXMod(NTL::ZZ_pE::cardinality(), F);
        auto g1 = NTL::GCD(g, h - x);
        for (auto const &r: NTL::FindRoots(g1))
            for (ssize_t i = 0; i < m; ++i)
                roots[idx++] = r;
    }
    return idx;
}

std::optional<NTL::ZZ_pE> sqrt(NTL::ZZ_pE const &alpha) 
// Take EXTREME caution to only call this when ZZ_pE agrees with alpha 
// This is the same as with all arithmetic in field extensions...
{
    NTL::ZZ_pEX f;
    NTL::ZZ_pE roots[2];
    SetCoeff(f, 2);
    f[0] = -alpha;
    if (myroots(roots, f) == 0) {
        return {};
    }
    return roots[0];
}

std::optional<NTL::ZZ_pE> sqrt(Fp2k Fext, NTL::ZZ_pE const &alpha) 
{
    NTL::ZZ_pEPush Push(Fext.F);
    return sqrt(alpha);
}

NTL::ZZ_pX _construct_f_from_f2(NTL::ZZ_pEX const &f2) {
    NTL::ZZ_pEX f2_conj;
    size_t deg = NTL::deg(f2);
    auto p = NTL::ZZ_p::modulus();

    NTL::SetCoeff(f2_conj, deg);
    for (size_t i = 0 ; i < deg ; i++) {
        f2_conj[i] = NTL::power(f2[i], p);
    }

    NTL::ZZ_pEX f_Fp2 = f2*f2_conj;
    NTL::ZZ_pX f;
    size_t deg_f = NTL::deg(f_Fp2);

    NTL::SetCoeff(f, deg_f);
    for (size_t i = 0 ; i < deg_f ; i++) {
        f[i] = NTL::rep(f_Fp2[i])[0];
    }

    return f;
}

NTL::mat_ZZ_p _compute_frobenius_matrix(unsigned degree) {
    NTL::mat_ZZ_p frob_action;
    NTL::ZZ p = NTL::ZZ_p::modulus();
    frob_action.SetDims(degree, degree);
    for (unsigned i = 0 ; i < degree ; i++) {
        NTL::ZZ_pX i_th_basis_poly;
        NTL::SetCoeff(i_th_basis_poly, i);
        NTL::ZZ_pE i_th_basis;
        NTL::conv(i_th_basis, i_th_basis_poly);
        NTL::ZZ_pE i_th_image = NTL::power(i_th_basis, p);
        NTL::ZZ_pX image_poly = NTL::rep(i_th_image);
        for (unsigned j = 0 ; j < degree ; j++) {
            if (NTL::deg(image_poly) >= (long) j) {
                frob_action[i][j] = NTL::rep(i_th_image)[j];
            } else {
                frob_action[i][j] = 0;
            }
        }
    }
    return frob_action;
}

// ZZ_pE must be Fp2
Fp2k::Fp2k(unsigned degree)
{
    assert(degree > 0);
    assert(NTL::ZZ_pE::degree() == 2);

    this->k = degree;

    NTL::ZZ_pEContext F;
    if (degree == 1) {
        F.save();
        NTL::ZZ_pX Fp2_gen_poly;
        NTL::ZZ_pE Fp2_gen;
        SetCoeff(Fp2_gen_poly, 1);
        conv(Fp2_gen, Fp2_gen_poly);
        this->F = std::move(F);
        this->Fp2_gen = Fp2_gen;
        this->Fp2_gen_nonzero = 1;
        this->frob_action = _compute_frobenius_matrix(2);
    }

    if (degree > 1) {
        NTL::ZZ_pEX f_2;
        NTL::ZZ_pX f;
        bool found = false;
        for (size_t i = 0 ; i < 1000 ; i++) {
            NTL::BuildIrred(f_2, degree);
            f = _construct_f_from_f2(f_2);
            
            if (NTL::SquareFreeDecomp(f)[0].b == 1) { //<- check if irreducible
                this->mod_Fp2 = f_2;
                found = true;
                break;
            }
        }
        if (!(found)) {
            throw;
        }

        auto Fp2_modulus = NTL::ZZ_pE::modulus();
        assert(NTL::deg(Fp2_modulus) == 2);

        NTL::ZZ_pE roots[2];
        { //init Fp2k temp
            NTL::ZZ_pEPush push(f);
            F.save();
            this->frob_action = _compute_frobenius_matrix(2*k);

            NTL::ZZ_pEX mod;

            // Convert Fp2 modulus to a polynomial over Fp2k
            NTL::conv(mod, Fp2_modulus);

            // Solve to find a generator
            myroots(roots, mod);

            this->F = std::move(F);

            // Have to pick the right one w.r.t the modulus we chose earlier
            auto r1 = roots[0];
        }

        NTL::ZZ_pEX r0;
        NTL::ZZ_pX i_poly;
        NTL::SetCoeff(i_poly, 1);
        NTL::ZZ_pE i;
        NTL::conv(i, i_poly);

        NTL::conv(r0, NTL::rep(roots[0]));
        NTL::ZZ_pEX r0_minpoly = NTL::MinPolyMod(r0 % this->mod_Fp2, this->mod_Fp2);
        if (r0_minpoly[0] == i) {
            this->Fp2_gen = roots[1]; // something got turned upside down, but w/e it works now
        } else {
            this->Fp2_gen = roots[0];
        }
        for (Fp2_gen_nonzero = 1; NTL::IsZero(NTL::rep(roots[0])[Fp2_gen_nonzero]); ++Fp2_gen_nonzero);
    }

    std::pair<NTL::ZZ_pE, NTL::ZZ_pE> starting_curve_coeffs;
    NTL::ZZ_pEX xmap_up;
    NTL::ZZ_pEX xmap_down;
    NTL::ZZ_pEX ymap_up;
    NTL::ZZ_pEX ymap_down;
    {
        NTL::ZZ p = NTL::ZZ_p::modulus();
        NTL::ZZ_pEPush push(this->F);
        if (p % 4 == 3) {
            this->iota_degree = 1;    
            this->starting_a = 1;
            this->starting_b = 0;
            NTL::ZZ_pE w = this->Fp2_gen;
            assert (w*w == -1);

            

            NTL::SetCoeff(xmap_up, 1);
            NTL::SetCoeff(xmap_down, 0);
            NTL::SetCoeff(ymap_up, 0);
            NTL::SetCoeff(ymap_down, 0);

            xmap_up[1] = -1;
            ymap_up[0] = w;

            this->maximal = true;
        } else if (p % 3 == 2) {
            this->iota_degree = 1;
            this->starting_a = 0;
            this->starting_b = 1;
            NTL::ZZ_pE w = this->Fp2_gen;
            assert (w*w == -3);

            NTL::SetCoeff(xmap_up, 1);
            NTL::SetCoeff(xmap_down, 0);
            NTL::SetCoeff(ymap_up, 0);
            NTL::SetCoeff(ymap_down, 0);

            xmap_up[1] = (w-1)/2;

            this->maximal = false;
        } else if (p % 7 == 3 || p % 7 == 5 || p % 7 == 6) {
            this->iota_degree = 2;
            this->starting_a = -35;
            this->starting_b = 98;
            NTL::ZZ_pE w = this->Fp2_gen;
            assert (w*w == -7);
            NTL::SetCoeff(xmap_up, 2);
            NTL::SetCoeff(xmap_down, 1);
            NTL::SetCoeff(ymap_up, 2);
            NTL::SetCoeff(ymap_down, 2);
            this->maximal = false;
            
            xmap_up[0] = 35*w - 63;
            xmap_up[0] /= 8;
            xmap_up[1] = w + 7;
            xmap_up[1] /= 4;
            xmap_up[2] = -w - 3;
            xmap_up[2] /= 8;
            xmap_down[0] = w - 7;
            xmap_down[0] /= 2;
            xmap_down[1] = 1;
            ymap_up[0] = -21*w - 119;
            ymap_up[0] /= 16;
            ymap_up[1] = -3*w + 7;
            ymap_up[1] /= 4;
            ymap_up[2] = w - 5;
            ymap_up[2] /= 16;
            ymap_down[0] = -7*w + 21;
            ymap_down[0] /= -2;
            ymap_down[1] = -w + 7;
            ymap_down[2] = -1;
        } else if (p % 11 == 2 || p % 11 == 6 || p % 11 == 7 || p % 11 == 8 || p % 11 == 10) {
            this->iota_degree = 3;
            this->starting_a = -1056;
            this->starting_b = 13552;
            NTL::ZZ_pE w = this->Fp2_gen;
            assert (w*w == -11);
            NTL::SetCoeff(xmap_up, 3);
            NTL::SetCoeff(xmap_down, 2);
            NTL::SetCoeff(ymap_up, 3);
            NTL::SetCoeff(ymap_down, 3);
            this->maximal = false;

            xmap_up[0] = 20768*w + 73568;
            xmap_up[0] /= 9;
            xmap_up[1] = -352*w - 1936;
            xmap_up[1] /= 3;
            xmap_up[2] = -4*w + 44;
            xmap_up[2] /= 3;
            xmap_up[3] = w - 5;
            xmap_up[3] /= 18;
            xmap_down[0] = 88*w + 440;
            xmap_down[1] = -4*w - 44;
            xmap_down[2] = 1;
            ymap_up[0] = -24640*w + 15488;
            ymap_up[0] /= 27;
            ymap_up[1] = 88*w - 2024;
            ymap_up[1] /= 9;
            ymap_up[2] = 10*w + 22;
            ymap_up[2] /= 3;
            ymap_up[3] = -w - 4;
            ymap_up[3] /= 27;
            ymap_down[0] = 2816*w + 7744;
            ymap_down[1] = -264*w - 1320;
            ymap_down[2] = -6*w - 66;
            ymap_down[2] /= -1;
            ymap_down[3] = -1;

        } else if (p % 19 == 2 || p % 19 == 3 || p % 19 == 8 || p % 19 == 10 || p % 19 == 12 || p % 19 == 13 || p % 19 == 14 || p % 19 == 15 || p % 19 == 18)  {
            this->iota_degree = 6;
            this->starting_a = -152;
            this->starting_b = 722;
            NTL::ZZ_pE w = this->Fp2_gen;
            assert (w*w == -19);
            NTL::SetCoeff(xmap_up, 5);
            NTL::SetCoeff(xmap_down, 4);
            NTL::SetCoeff(ymap_up, 6);
            NTL::SetCoeff(ymap_down, 6);
            this->maximal = false;

            xmap_up[0] = 46208*w + 219488;
            xmap_up[0] /= 5;
            xmap_up[1] = -17328*w - 103968;
            xmap_up[1] /= 5;
            xmap_up[2] = 2128*w + 17328;
            xmap_up[2] /= 5;
            xmap_up[3] = -76*w - 1216;
            xmap_up[3] /= 5;
            xmap_up[4] = -2*w + 38;
            xmap_up[4] /= 5;
            xmap_up[5] = w - 9;
            xmap_up[5] /= 50;
            xmap_down[0] = 31768*w + 147288;
            xmap_down[0] /= 25;
            xmap_down[1] = -456*w - 2888;
            xmap_down[1] /= 1;
            xmap_down[2] = 266*w + 2546;
            xmap_down[2] /= 5;
            xmap_down[3] = -2*w - 38;
            xmap_down[3] /= 1;
            xmap_down[4] = 1;
            xmap_down[4] /= 1;
            ymap_up[0] = 877952*w + 3072832;
            ymap_up[0] /= 125;
            ymap_up[1] = -600704*w - 658464;
            ymap_up[1] /= 125;
            ymap_up[2] = 112632*w - 147288;
            ymap_up[2] /= 125;
            ymap_up[3] = -608*w + 54872;
            ymap_up[3] /= 125;
            ymap_up[4] = -1672*w - 5852;
            ymap_up[4] /= 125;
            ymap_up[5] = 27*w + 57;
            ymap_up[5] /= 25;
            ymap_up[6] = -2*w - 7;
            ymap_up[6] /= 125;
            ymap_down[0] = 18875968*w + 38629888;
            ymap_down[0] /= -125;
            ymap_down[1] = -450528*w - 1316928;
            ymap_down[1] /= -5;
            ymap_down[2] = 528504*w + 2174664;
            ymap_down[2] /= -25;
            ymap_down[3] = -2432*w - 14440;
            ymap_down[3] /= -1;
            ymap_down[4] = 684*w + 6384;
            ymap_down[4] /= -5;
            ymap_down[5] = -3*w - 57;
            ymap_down[5] /= -1;
            ymap_down[6] = -1;
        } else {
            throw; //Have to be extremely unlucky
        }
    }

    this->iota_x_num = xmap_up;
    this->iota_x_denom = xmap_down;
    this->iota_y_num = ymap_up;
    this->iota_y_denom = ymap_down;
    std::cout << "Done with one!" << std::endl;
}

NTL::ZZ_pE lift(NTL::ZZ_pE const &alpha, Fp2k const &Fext)
/*
Takes an element alpha of Fp2, and an extension Fp2k,
returns the of alpha in Fp2k
*/
{
    if (Fext.k == 1) {
        return alpha;
    }

    if (alpha == 0) {
        return alpha;
    }
    
    NTL::ZZ_pE alpha_bar;

    NTL::ZZ_p a = NTL::rep(alpha)[0], b;
    if (NTL::deg(NTL::rep(alpha)) == 1) {
        b = NTL::rep(alpha)[1];
    } else {
        b = 0;
    }

    {
        NTL::ZZ_pEPush push(Fext.F);
        alpha_bar = a + Fext.Fp2_gen*b;
    }
    return alpha_bar;
}



NTL::ZZ_pE coerce(NTL::ZZ_pE const &alpha, Fp2k const &Fext) //super confusing that the first ZZ_pE here refers to Fp2, and the second refers to Fp2k...
/*
Takes an element alpha of Fp2k, and an extension Fp2k,
returns the of alpha in Fp2k
*/
{
    if (Fext.k == 1){ // Nothing to be done in this case.
        return alpha;
    }

    NTL::ZZ_pX alpha_poly;
    NTL::ZZ_pE alpha_low;
    SetCoeff(alpha_poly, 1);
    {
        NTL::ZZ_pEPush push(Fext.F);
        alpha_poly[0] = trace(alpha) / (2*Fext.k);
        auto const idx = Fext.Fp2_gen_nonzero;
        alpha_poly[1] = NTL::coeff(NTL::rep(alpha), idx) / NTL::rep(Fext.Fp2_gen)[idx];
    }
    conv(alpha_low, alpha_poly);

    //assert (lift(alpha_low, Fext) == alpha);

    return alpha_low;
}

NTL::ZZ_pEX lift(NTL::ZZ_pEX const &f, Fp2k const &Fext)
/*
Takes a polynomial over Fp2, and an extension Fp2k,
returns the polynomial in Fp2k.
*/
{
    if (Fext.k == 1){ // Nothing to be done in this case.
        return f;
    }
    if (f == 0){ // 0 is zero
        return f;
    }
    NTL::ZZ_pEPush push(Fext.F);

    NTL::ZZ_pEX f_high;
    size_t d = NTL::deg(f);

    NTL::SetCoeff(f_high, d);
    
    f_high[0] = lift(coeff(f, 0), Fext);
    for (size_t i = 1; i <= d; ++i){
        f_high[i] = lift(coeff(f, i), Fext);
    }

    return f_high;
}


NTL::ZZ_pEX coerce(NTL::ZZ_pEX const &f, Fp2k const &Fext)
/*
Takes a polynomial over Fp2k, and an extension Fp2k,
returns the polynomial in Fp2. Assumes it exists.
*/
{
    if (Fext.k == 1){ // Nothing to be done in this case.
        return f;
    }

    NTL::ZZ_pEX f_low;
    size_t d = NTL::deg(f);

    NTL::SetCoeff(f_low, d);
    
    f_low[0] = coerce(coeff(f, 0), Fext);
    for (size_t i = 1; i <= d; ++i){
        f_low[i] = coerce(coeff(f, i), Fext);
    }

    return f_low;
}

NTL::ZZ_pE Fp2k::frob(NTL::ZZ_pE alpha) const
{ // Lots of conversion a bit annoying...
    NTL::vec_ZZ_p alpha_vec;
    NTL::conv(alpha_vec, NTL::rep(alpha));
    alpha_vec.SetLength(this->k*2);

    NTL::vec_ZZ_p image_alpha_vec = alpha_vec*this->frob_action;

    NTL::ZZ_pE image_alpha;
    NTL::ZZ_pX image_alpha_poly;
    NTL::conv(image_alpha_poly, image_alpha_vec);
    NTL::conv(image_alpha, image_alpha_poly);

    return image_alpha;
}

NTL::ZZ_pEX MinPoly(NTL::ZZ_pE const &alpha, Fp2k const &Fext)
/*
Computes the min poly of alpha \in Fext over Fp2
*/
{
    if (Fext.k == 1) {
        NTL::ZZ_pEX f;
        NTL::SetCoeff(f, 1);
        f[0] = -alpha;
        return f;
    }

    NTL::ZZ_pEX a_min_poly;

    NTL::ZZ_pEX alpha_over_Fp2;
    NTL::conv(alpha_over_Fp2, NTL::rep(alpha));

    a_min_poly = NTL::MinPolyMod(alpha_over_Fp2 % Fext.mod_Fp2, Fext.mod_Fp2);

    return a_min_poly;
}
