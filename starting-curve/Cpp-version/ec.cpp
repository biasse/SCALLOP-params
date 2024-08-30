#include <cassert>
#include <optional>

#include <NTL/ZZ_pE.h>
#include <NTL/ZZ.h>
#include "Fp2k.hpp"
#include "ec.hpp"
#include "ecp.hpp"
#include "utils.hpp"

ec::ec(NTL::ZZ_pE const &a, NTL::ZZ_pE const &b) : _a{a}, _b{b}
{
    if ((4*this->_a*this->_a*this->_a + 27*this->_b*this->_b) == 0)
        throw std::logic_error("curve is singular");
}

ecp ec::operator()(Fp2k const &Fext, NTL::ZZ_pE new_x, NTL::ZZ_pE new_y) const {
    return ecp(std::make_shared<const ec>(*this), Fext, new_x, new_y);
}

std::optional<ecp> ec::lift_x(Fp2k const &Fext, NTL::ZZ_pE const &x) const
{
    NTL::ZZ_pEPush(Fext.F);

    NTL::ZZ_pE rhs = (x*x + lift(this->_a, Fext))*x + lift(this->_b, Fext);
    auto y = sqrt(rhs);
    if (y)
        return ecp(std::make_shared<const ec>(*this), Fext, x, *y);
    return {};
}

ecp ec::random_point(Fp2k const &Fext) const
{

    std::cout << "Running random_point" << std::endl;
    NTL::ZZ_pEPush Push(Fext.F);
    while (true) {
        auto pt = this->lift_x(Fext, NTL::random_ZZ_pE());
        if (pt) {
            std::cout << "Done random_point" << std::endl;
            return *pt;
        }
    }
}

NTL::ZZ_pE ec::j_invariant() const
{
    auto &a = _a, &b = _b;
    auto a3 = NTL::ZZ_pE(4)*a*a*a;
    auto b2 = NTL::ZZ_pE(27)*b*b;
    return NTL::ZZ_pE(1728) * a3/(a3 + b2);
}

ec ec::from_j(NTL::ZZ_pE const &j)
{
    if (j == 0) { return ec(NTL::ZZ_pE(0), NTL::ZZ_pE(1));}
    NTL::ZZ_pE const f(1728);
    if (j == f) { return ec(NTL::ZZ_pE(1), NTL::ZZ_pE(0));}
    NTL::ZZ_pE two(2), three(3);
    auto j2 = j*j;
    auto j3 = j2*j;
    NTL::ZZ_pE a = three * (f*j - j2);
    NTL::ZZ_pE b = two * (-j3 + two*f*j2 - f*f*j);
    auto E = ec(a, b);
    assert(E.j_invariant() == j);
    return E;
}

ecp ec::random_point_of_order(Fp2k const &Fext, NTL::ZZ cof, int ell, int e) const {
    ecp P = this->random_point(Fext);
    P = cof*P;
    while (!(NTL::power(NTL::ZZ(ell), e-1)*P)) {
        P = this->random_point(Fext);
        P = cof*P;
    }

    assert (NTL::power(NTL::ZZ(ell), e-1)*P);
    assert (!(NTL::power(NTL::ZZ(ell), e)*P));
    return P;
}

std::pair<ecp, ecp> const ec::torsionBasis(Fp2k const &Fext, int ell, int e) const
{
    // What if we are on the twist? 
    NTL::ZZ cof = NTL::power(NTL::ZZ_p::modulus(), long(Fext.k)) - NTL::power_long(long(-1), long(Fext.k % 2));

    NTL::ZZ ellcof = NTL::power(NTL::ZZ(ell), e-1);
    NTL::ZZ le = ellcof * ell;
    assert (cof % le == 0);
    cof /= le;

    ecp P = this->random_point_of_order(Fext, cof, ell, e);

    std::cout << "Found P, looking for Q..." << std::endl;
    auto ellP = ellcof*P;

    while (true) {
        ecp Q = this->random_point_of_order(Fext, cof, ell, e);
        std::cout << "Found potential Q..." << std::endl;
        if (DLP(ellcof*Q, ellP, ell, 1) == NTL::ZZ(-1))
            return {P, Q};
    }
}
