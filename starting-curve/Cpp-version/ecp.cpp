#include <cassert>

#include <optional>

#include <NTL/ZZ_pE.h>
#include "Fp2k.hpp"
#include "ec.hpp"
#include "ecp.hpp"
#include "utils.hpp"


ecp::ecp(std::shared_ptr<const ec> curve, Fp2k const &Fext, NTL::ZZ_pE const &xx, NTL::ZZ_pE const &yy)
    : E{curve}, Fext{Fext}, x{xx}, y{yy}, z{NTL::ZZ_pE(1)}
{
#ifndef NDEBUG
    NTL::ZZ_pEPush Push(Fext.F);
    ec E = *curve.get();
    assert(y*y == x*x*x + lift(E.a(), Fext)*x + lift(E.b(), Fext));
#endif
}

ecp::ecp(std::shared_ptr<const ec> curve, Fp2k const &Fext, NTL::ZZ_pE const &xx, NTL::ZZ_pE const &yy, NTL::ZZ_pE const &zz)
    : E{curve}, Fext{Fext}, x{xx}, y{yy}, z{zz}
{
#ifndef NDEBUG
    NTL::ZZ_pEPush Push(Fext.F);
    assert(NTL::deg(NTL::ZZ_pE::modulus()) == 2*Fext.k);
    ec E = *curve.get();
    assert(y*y*z == x*x*x + lift(E.a(), Fext)*x*z*z + lift(E.b(), Fext)*z*z*z);
#endif
}

ecp operator*(Integer k, ecp P) {
    if (k < 0) {
        return (-k)*(-P);
    }

    ecp result(P.E, P.field());
    ecp Q = P;

    while (k > 0) {
        if (k % 2 == 1) {
            result += Q;
        }
        Q += Q;
        k >>= 1;
    }

    return result;
}


void ecp::normalize() const
{
    NTL::ZZ_pEPush push(field().F);
    if (NTL::IsZero(z)) {
        NTL::set(y);
    }
    else {
        auto invz = NTL::inv(z);
        x *= invz;
        y *= invz;
        NTL::set(z);
    }
}

std::pair<NTL::ZZ_pE, NTL::ZZ_pE> ecp::_lambda(ecp const &other) const
{
    //Notice this is (should?) only called from +, where ZZ_pEPush is already done, so we dont have to do it again
    auto const &P = *this, &Q = other;
    NTL::ZZ_pE const &x1 = P.x, &y1 = P.y, &z1 = P.z;
    NTL::ZZ_pE const &x2 = Q.x, &y2 = Q.y, &z2 = Q.z;

    if (x1*z2 == x2*z1) {
        NTL::ZZ_pE x1sq = x1 * x1;
        NTL::ZZ_pE z1sq = z1 * z1;
        return {x1sq+x1sq+x1sq + lift(this->curve().a(), P.field())*z1sq, (y1+y1) * z1};
    }

    return {y2*z1 - y1*z2, x2*z1 - x1*z2};
}

ecp ecp::operator+(ecp other) const
{
    if (this->curve() != other.curve())
        throw std::logic_error("points not on the same curve");

    if (!other) return *this;
    if (!*this) return other;


    if (&this->field() != &other.field())
        throw std::logic_error("points not over the same field"); 
    
    NTL::ZZ_pEPush push(this->field().F);

    NTL::ZZ_pE const &x1 = this->x, &y1 = this->y, &z1 = this->z;
    NTL::ZZ_pE const &x2 = other.x, &y2 = other.y, &z2 = other.z;

    if (x1*z2 == x2*z1 && y1*z2 != y2*z1) {
        assert(y2*z1 == -y1*z2);
        return {this->E, this->Fext};
    }

    auto lam = this->_lambda(other);

    NTL::ZZ_pE const &u = lam.first, &v = lam.second;

    NTL::ZZ_pE uu = u*u, vv = v*v;
    NTL::ZZ_pE zz = z1 * z2;

    NTL::ZZ_pE z3_ = vv*zz;
    NTL::ZZ_pE x3_ = uu*zz - vv*(x1*z2 + x2*z1);

    NTL::ZZ_pE vz1 = v*z1;
    NTL::ZZ_pE y3 = u*(x1*z3_ - x3_*z1) - y1*v*z3_;
    NTL::ZZ_pE x3 = x3_*vz1;
    NTL::ZZ_pE z3 = z3_*vz1;
    return ecp(this->E, this->Fext, x3,y3,z3);
}
