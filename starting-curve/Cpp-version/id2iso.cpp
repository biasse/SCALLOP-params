#include <NTL/ZZ_pE.h>

#include "quaternions.hpp"
#include "id2iso.hpp"
#include "isog.hpp"
#include "utils.hpp"
#include "Fp2k.hpp"

int torsionToFieldDegree(Integer ell_e) {
    //TODO: Optimize this, and account for twist. (Can just take the smallest k s.t. p^k = \pm 1 in Z_ell^e)
    auto const &p = NTL::ZZ_p::modulus();
    NTL::ZZ cur = p;
    for (size_t k = 1; k < 100; k++ ) {
        NTL::ZZ ord_k = cur - (k % 2 ? -1 : +1);
        if (ord_k % ell_e == 0)
            return k;
        cur = cur * p % ell_e;
    }
    // std::cout << "Tried all k's, returning -1" << std::endl;
    // throw; // MARIA: I changed this to return -1 because it kept throwing errors for me. I tried to increase k but it still couldnt find one, so decided to change it to return -1 in these cases and simply discard the ell_e. a TODO is to make this more robust.
    return -1;
}

ecp endo_i(ecp const &P)
{
    Fp2k F = P.field();
    assert ((P.curve().a() == F.starting_a) && (P.curve().b() == F.starting_b));  // only makes sense for the starting curve

    NTL::ZZ_pEPush push(P.field().F);

    NTL::ZZ_pE xmap_up = NTL::coeff(F.iota_x_num, 0);
    NTL::ZZ_pE xmap_down = NTL::coeff(F.iota_x_denom, 0);
    NTL::ZZ_pE ymap_up = NTL::coeff(F.iota_y_num, 0);
    NTL::ZZ_pE ymap_down = NTL::coeff(F.iota_y_denom, 0);

    NTL::ZZ_pE xi(1);

    P.normalize();
    NTL::ZZ_pE x = P.get_x();
    NTL::ZZ_pE y = P.get_y();
    for (size_t i = 1; i <= F.iota_degree; i++) { //NTL::coeff returns zero when i > degree
        xi *= x;
        xmap_up += xi*NTL::coeff(F.iota_x_num, i);
        xmap_down += xi*NTL::coeff(F.iota_x_denom, i);
        ymap_up += xi*NTL::coeff(F.iota_y_num, i);
        ymap_down += xi*NTL::coeff(F.iota_y_denom, i);
    }


    ecp wP = ecp(P.curve_ptr(), P.field(), xmap_up/xmap_down, y*ymap_up/ymap_down);

    if (F.maximal) {
        return wP;
    } else {
        return 2 * wP + P; //This could be optimised to instead write alpha in this "better" basis (1 + w)/2
    }
}

ecp evalEndo(quat const &alpha, ecp P, Integer order_P)
{
    if ((alpha.alg.q != 1) && (alpha.alg.q != 3) && (alpha.alg.q != 7) && (alpha.alg.q != 11) && (alpha.alg.q != 19) && (alpha.alg.q != 43))
        throw;  // not currently implemented
    assert(alpha.alg.p == NTL::ZZ_p::modulus());

    assert(!(order_P*P));

    ecp iP = endo_i(P),
        jP = P.frob(),
        kP = endo_i(jP);

    assert (-endo_i(iP) == alpha.alg.q*P);
    assert (-iP.frob() = kP);
    assert (alpha[4] == 1);

    ecp alpha_P = alpha[0] % order_P * P;
    alpha_P += alpha[1] % order_P * iP;
    alpha_P += alpha[2] % order_P * jP;
    alpha_P += alpha[3] % order_P * kP;

    return alpha_P;
}

std::vector<std::pair<ecp,std::pair<int, int>>> idealToKernel(quat const &alpha, Integer const &N, ec const &E0, std::map<unsigned,Fp2k> &FieldExtensions, std::map<NTL::ZZ,std::pair<ecp,ecp>> &TorsionBases) {
    // I = O0<alpha, N>
    // Field extensions should maybe be something like an unordered map instead?
    std::unordered_map<int, int> facN = factor(N);

    std::vector<std::pair<ecp,std::pair<int, int>>> kerGens;
    for (const auto& [ell,e] : facN) {

        auto ell_e_m_1 = NTL::power(NTL::ZZ(ell), e-1);
        auto ell_e = ell_e_m_1 * ell;

        assert((2*alpha.alg.q) % alpha[4] == 0);
        bool extra = (alpha[4] % 2 == 0 && ell == 2) || (alpha[4] % alpha.alg.q == 0 && ell == alpha.alg.q);
        auto ell_ext = ell_e * (extra ? ell : 1);

        quat endo = alpha.conjugate() * alpha[4];

        auto it = TorsionBases.find(ell_ext);
        if (it == TorsionBases.end()) {
            unsigned k = torsionToFieldDegree(ell_ext);
            std::cout << "Constructing torsion basis of E[" << ell << "^" << e << "] over Fp2k, k = " << k << std::endl;
            auto jt = FieldExtensions.find(k);
            assert(jt != FieldExtensions.end());  //TODO
//            if (jt == FieldExtensions.end())
//                jt = FieldExtensions.emplace(std::make_pair(k, Fp2k {k})).first;
            assert(jt->first == jt->second.k);
            auto bas = E0.torsionBasis(jt->second, ell, e + extra);
            it = TorsionBases.emplace(std::make_pair(ell_ext, bas)).first;
        }
        std::cout << "Finding kernel generator of order" << ell << "^" << e << std::endl;
        auto const &Basis = it->second;
        assert(!(ell_ext * Basis.first) && ell_ext/ell * Basis.first);
        assert(!(ell_ext * Basis.second) && ell_ext/ell * Basis.second);

        auto const project = [&](ecp const &pt) {
            return evalEndo(endo, pt, ell_ext);
        };

        ecp K = project(Basis.first);
        assert(!(ell_e*K));
        if (!(ell_e_m_1*K))
            K = project(Basis.second);
        assert(!(ell_e*K));
        assert((ell_e_m_1)*K);

        std::cout << "Done!" << std::endl;
        kerGens.push_back(std::pair<ecp,std::pair<int, int>>(K, {ell,e}));
    }
    return kerGens;
}

isog_chain idealToIsogeny(quat const &alpha, Integer const &N, ec const &E0, std::map<unsigned,Fp2k> &FieldExtensions, std::map<NTL::ZZ,std::pair<ecp,ecp>> &TorsionBases) {
    auto kerGens = idealToKernel(alpha, N, E0, FieldExtensions, TorsionBases);
    std::sort(kerGens.begin(), kerGens.end(), [](auto &left, auto &right) { //It helps alot to sort the kernel generators like this to avoid evaluating expensive isogenies as much as possible
        return left.second.first > right.second.first;
    });
    //std::cout << "Computing isogeny chain from kernel generators" << std::endl;
    return isog_chain(kerGens);
}
