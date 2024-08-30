
#pragma once

#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pEX.h>
#include <optional>

size_t myroots(NTL::ZZ_pE *roots, NTL::ZZ_pEX const &f);
std::optional<NTL::ZZ_pE> sqrt(NTL::ZZ_pE const &alpha);

struct Fp2k
{
    NTL::ZZ_pEContext F;        // The field extension
    NTL::ZZ_pE Fp2_gen;         // The (image of the) generator of Fp2
    unsigned k;                 // The degree as an extension of Fp2
    size_t Fp2_gen_nonzero;
    NTL::ZZ_pEX mod_Fp2;
    NTL::mat_ZZ_p frob_action;    // For large fields this uses quite a bit of memory, but probably worth it
    Fp2k(unsigned degree);
    Fp2k() = delete;

    NTL::ZZ_pE frob(NTL::ZZ_pE alpha) const;

    // Compute the iota map when generating the fields, as a reduction of curve with CM
    // Could consider doing this only for Fp2 (currently this essentially precomputes the lifting for each field extension, using alot of memory)
    NTL::ZZ_pE starting_a;
    NTL::ZZ_pE starting_b;

    NTL::ZZ_pEX iota_x_num;
    NTL::ZZ_pEX iota_x_denom;
    NTL::ZZ_pEX iota_y_num;
    NTL::ZZ_pEX iota_y_denom;

    bool maximal; // whether to use iota or 2iota + 1

    size_t iota_degree;
};

std::optional<NTL::ZZ_pE> sqrt(Fp2k Fext, NTL::ZZ_pE const &alpha);
NTL::ZZ_pE lift(NTL::ZZ_pE const &alpha, Fp2k const &Fext);
NTL::ZZ_pE coerce(NTL::ZZ_pE const &alpha, Fp2k const &Fext);
NTL::ZZ_pEX lift(NTL::ZZ_pEX const &f, Fp2k const &Fext);
NTL::ZZ_pEX coerce(NTL::ZZ_pEX const &f, Fp2k const &Fext);
NTL::ZZ_pEX MinPoly(NTL::ZZ_pE const &alpha, Fp2k const &Fext);

