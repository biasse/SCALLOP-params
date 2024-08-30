#pragma once

#include <optional>

#include <NTL/ZZ_pE.h>
#include "Fp2k.hpp"

class ecp;

class ec
{
    private:
        NTL::ZZ_pE _a, _b;
    public:
        ec() {};
        ec(NTL::ZZ_pE const &a, NTL::ZZ_pE const &b);

        NTL::ZZ_pE const &a() const { return this->_a; }
        NTL::ZZ_pE const &b() const { return this->_b; }
        // Kohel: b2 = 0, b4 = 2*a, b6 = 4*b
    
    std::optional<ecp> lift_x(Fp2k const &Fext, NTL::ZZ_pE const &x) const;

    ecp random_point(Fp2k const &Fext) const;
    ecp random_point_of_order(Fp2k const &Fext, NTL::ZZ cof, int ell, int k) const;

    bool operator==(ec const &other) const { return this == &other || (this->_a == other._a && this->_b == other._b); }
    bool operator!=(ec const &other) const { return !(*this == other); }

    ecp operator()(Fp2k const &Fext, NTL::ZZ_pE new_x, NTL::ZZ_pE new_y) const;

    NTL::ZZ_pE j_invariant() const;

    ec static from_j(NTL::ZZ_pE const &j);

    std::pair<ecp, ecp> const torsionBasis(Fp2k const &Fext, int ell, int k) const;

    friend std::ostream& operator<<(std::ostream& o, ec const &E) { return o << "{y^2 = x^3 + (" << NTL::rep(E._a) << ")*x + (" << NTL::rep(E._b) << ")}"; }

};
