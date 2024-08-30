#pragma once

#include <cassert>
#include <optional>
#include <memory> // Include the memory header for shared_ptr

#include <NTL/ZZ_pE.h>
#include "Fp2k.hpp"
#include "ec.hpp"
#include "utils.hpp"

class ecp {
    private:
        std::shared_ptr<const ec> E;
        std::reference_wrapper<const Fp2k> Fext;
        mutable NTL::ZZ_pE x, y, z;
    public:
        // Given only a curve and field extension, it creates point at infinity
        // ecp();
        ecp(std::shared_ptr<const ec> curve, Fp2k const &Fext) : E{curve}, Fext{Fext}, x{0}, y{1}, z{0} {};
        ecp(std::shared_ptr<const ec> curve, Fp2k const &Fext, NTL::ZZ_pE const &x, NTL::ZZ_pE const &y, NTL::ZZ_pE const &z);
        ecp(std::shared_ptr<const ec> curve, Fp2k const &Fext, NTL::ZZ_pE const &x, NTL::ZZ_pE const &y);

        NTL::ZZ_pE const &get_x() const {return this->x;}
        NTL::ZZ_pE const &get_y() const {return this->y;}
        NTL::ZZ_pE const &get_z() const {return this->z;}

        NTL::ZZ_pE const aff_x() const {if (this->is_identity()) {return NTL::ZZ_pE(0);} NTL::ZZ_pEPush push(this->field().F); return this->x/this->z;}
        NTL::ZZ_pE const aff_y() const {if (this->is_identity()) {return NTL::ZZ_pE(1);} NTL::ZZ_pEPush push(this->field().F); return this->y/this->z;}

        bool is_identity() const {return NTL::IsZero(this->z);}

    friend std::ostream& operator<<(std::ostream& o, ecp const &P)
    {
        return o << "(" << NTL::rep(P.x)
                 << " : " << NTL::rep(P.y)
                 << " : " << NTL::rep(P.z)
                 << ")";
    }

    ec const &curve() const { return *E; }
    std::shared_ptr<const ec> const &curve_ptr() const { return E; }
    Fp2k const &field() const { return this->Fext.get(); }

    void normalize() const;
    ecp const &normalized() { this->normalize(); return *this;}

    // For addition methods: Should there be some kind of promotion functionality
    //...when adding a point with another point in a subfield? // composite of fields will have degree lcm
    //...hmm complicated, but maybe when adding a point over Fp2k with one over Fp2? That should be easy ish. // agreed
    std::pair<NTL::ZZ_pE,NTL::ZZ_pE> _lambda(ecp const &other) const;
    ecp operator+(ecp other) const;
    friend ecp operator*(Integer k, ecp P);
    friend ecp operator*(int k, ecp P) { return Integer(k)*P; }

    ecp const frob() const
    {
        auto const &fld = field();
        NTL::ZZ_pEPush push(fld.F);
        return {this->E, Fext, fld.frob(x), fld.frob(y), fld.frob(z)};
    }


    bool operator==(ecp const &other) const
    {
        return (this->curve() == other.curve() 
        && this->aff_x() == other.aff_x() 
        && this->aff_y() == other.aff_y()
        && this->is_identity() == other.is_identity()); 
    }

    ecp &operator+=(ecp const &other) { return *this = *this + other; }
    ecp operator-() const { ecp Q = *this; Q.y = -Q.y; return Q; }
    ecp operator-(ecp const &other) const { return *this + (-other); }
    ecp &operator-=(ecp const &other) { return *this = *this - other; }
    operator bool() const { return !(this->is_identity());}
};


namespace std {
    //Is there any more fancy requirements on this? I have no idea what Im doing here.
    template<>
    struct hash<NTL::ZZ_pE>
    {
        size_t operator()(NTL::ZZ_pE const &a) const {
            size_t h = 0xb00b;
            if (a == 0) {
                return h;
            }
            auto apoly = NTL::rep(a);
            size_t d = NTL::deg(apoly);
            for (size_t i = 0; i <= d; i++) {
                size_t c;
                conv(c, NTL::coeff(apoly, i));
                h ^= c;
            }
            return h;
        }
    };
}

namespace std {
    //Is there any more fancy requirements on this?
    template<>
    struct hash<ecp>
    {
        size_t operator()(ecp const &P) const {
            size_t h = 0xec77;
            h += std::hash<NTL::ZZ_pE>()(P.curve().a());
            h += std::hash<NTL::ZZ_pE>()(P.curve().b());

            h ^= std::hash<NTL::ZZ_pE>()(P.aff_x());
            h *= 1111111;
            h ^= std::hash<NTL::ZZ_pE>()(P.aff_y());
            return h;
        }
    };
}

