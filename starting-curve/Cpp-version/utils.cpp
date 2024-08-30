#include <cassert>
#include <unordered_map>
#include <vector>
#include <cassert>

#include <NTL/ZZ.h>
#include <NTL/ZZ_pEX.h>
#include <NTL/ZZ_pEXFactoring.h>
#include <NTL/ZZ_pXFactoring.h>

#include "Fp2k.hpp"
#include "ec.hpp"
#include "ecp.hpp"
#include "utils.hpp"

NTL::ZZ_pEX modPoly2(NTL::ZZ_pE const &j)
{
    // Construct the modular polynomial of level 2 evaluated at j
    NTL::ZZ_pEX Phi;
    NTL::SetCoeff(Phi, 3);

    auto j2 = j*j;
    auto j3 = j2*j;

    Phi[0] = j3 - 162000*j2 + 8748000000*j - 157464000000000;
    Phi[1] = 1488*j2 + 40773375*j + 8748000000;
    Phi[2] = -j2 + 1488*j - 162000;

    return Phi;
}

static size_t log2floor(NTL::ZZ const &n)
{
    return NTL::NumBits(n) - 1;
}


bool sutherland(NTL::ZZ_pE const &j)
{
    // Sutherlands algorithm for supersingularity testing
    // Extremely fast when j is NOT supersingular
    // Following eprint
    // CAUTION: Does not work for very small p < 100

    if (j == 1728 && NTL::ZZ_p::modulus() % 4 == 3){
        return true;
    }

    NTL::ZZ_pE j_old[3] {j,j,j}, j_new[3];
    NTL::ZZ_pE j_ext_old, j_ext_new;
    bool found_j_ext = false;

    NTL::ZZ_pEX Phi_2 = modPoly2(j);

    // First take a step in each direction
    if (myroots(j_new, Phi_2) < 3) {
        return false;
    }

    //check if there is a j-inv not defined over F_p
    for (int i = 0; i < 3; i++) {
        if (NTL::deg(NTL::rep(j_new[i])) == 1){
            j_ext_old = j_old[i];
            j_ext_new = j_new[i];
            found_j_ext = true;
        }
    }

    //if not, take one more step in each path
    if (!found_j_ext){
        for (int i = 0; i < 3; i++) {

            if (found_j_ext){
                continue;
            }

            Phi_2 = modPoly2(j_new[i]);

            // divide out by the previous step
            NTL::ZZ_pEX prev;
            NTL::SetCoeff(prev, 1);
            prev[0] = -j_old[i];
            Phi_2 /= prev;

            // Check if we're at the floor
            NTL::ZZ_pE roots[2];
            if (myroots(roots, Phi_2) < 2) {
                return false;
            }

            // break if we found a j-inv over F_p^2
            for (int ii = 0; ii < 2; ii++) {
                if (NTL::deg(NTL::rep(roots[ii])) == 1){
                    j_ext_old = j_new[i];
                    j_ext_new = roots[ii];
                    found_j_ext = true;
                }
            }    
        }
    }

    //must have found a j-inv over F_p^2 in the supersingular case at this point 
    if (!found_j_ext){
        return false;
    }

    // Then keep going along this path, either down the 2-isog vulcano, or just around the supersingular graph
    size_t bound = log2floor(NTL::ZZ_p::modulus())/2;
    for (size_t s = 0; s < bound + 1; s++) {  

        Phi_2 = modPoly2(j_ext_new);

        // divide out by the previous step
        NTL::ZZ_pEX prev;
        NTL::SetCoeff(prev, 1);
        prev[0] = -j_ext_old;
        Phi_2 /= prev;

        // Check if we're at the floor
        NTL::ZZ_pE roots[2];
        if (myroots(roots, Phi_2) < 2) {
            return false;
        }

        // update j-lists
        j_ext_old = j_ext_new;
        j_ext_new = roots[0];
    }

    return true;
}

bool sutherland_slow(NTL::ZZ_pE const &j)
{
    // Sutherlands algorithm for SS testing
    // Extremely fast when j is NOT supersingular

    NTL::ZZ_pE j_old[3] {j,j,j}, j_new[3];

    NTL::ZZ_pEX Phi_2 = modPoly2(j);

    // First take a step in each direction
    if (myroots(j_new, Phi_2) < 3) {
        return false;
    }

    // Then keep going, either down the 2-isog vulcano, or just around the SS graph
    size_t bound = log2floor(NTL::ZZ_p::modulus());
    for (size_t s = 0; s < bound; s++) {
        /*std::cout << "\n\nCurrent j-list:\n";
        std::cout << "j-old:\n";
        std::cout << j_old[0] << j_old[1] << j_old[2] << std::endl;
        std::cout << "j-new:\n";
        std::cout << j_new[0] << j_new[1] << j_new[2] << std::endl; */
        for (int i = 0; i < 3; i++) {

            Phi_2 = modPoly2(j_new[i]);

            // divide out by the previous step
            NTL::ZZ_pEX prev;
            NTL::SetCoeff(prev, 1);
            prev[0] = -j_old[i];
            Phi_2 /= prev;

            // Check if we're at the floor
            NTL::ZZ_pE roots[2];
            if (myroots(roots, Phi_2) < 2) {
                return false;
            }

            // update j-lists
            j_old[i] = j_new[i];
            j_new[i] = roots[0];
        }
    }

    return true;
}


int _BSGS(ecp const &Q, ecp const &base, int base_order) {

    assert (Q.curve() == base.curve());
    assert (&Q.field() == &base.field());
    NTL::ZZ_pEPush push(Q.field().F);

    size_t m = NTL::SqrRoot(base_order);

    std::unordered_map<ecp, size_t> baby_steps;

    ecp G = 0*base;

    baby_steps[G] = 0;


    for (size_t i = 1; i <= m+1; i++) {
        G += base;
        baby_steps[G] = i;
    }

    ecp min_G = -G;
    ecp gamma = Q;

    for (size_t i = 0; i < m; i++) {
        auto hit = baby_steps.find(gamma);

        if (hit != baby_steps.end())
            return hit->second + m*i;
        gamma += min_G;
    }

    return -1;
}

NTL::ZZ DLP(ecp const &Q, ecp const &base, int ell, int e) {
    ecp gamma = NTL::power(NTL::ZZ(ell), e-1)*base;

    NTL::ZZ ell_to_k(1);
    NTL::ZZ x(0);
    for (int k = 0; k < e; k++) {
        ecp H = NTL::power(NTL::ZZ(ell), e-1-k)*(Q - x*base);
        int d = _BSGS(H, gamma, ell);
        if (d == -1) {
            return NTL::ZZ(-1);
        }
        x += ell_to_k*d;
        ell_to_k *= ell;
    }
    return x;
}

std::unordered_map<int, int> factor(Integer const &N) {
    // Trial division for now
    int p = 1;

    std::unordered_map<int, int> factored;
    auto Ni = N;
    while (Ni > 1) {
        p = NTL::NextPrime(p+1);
        bool first = true;
        while ((Ni % p) == 0) {
            Ni /= p;
            if (first) {
                factored[p] = 1;
                first = false;
            } else {
                factored[p] += 1;
            }
        }
    }

    return factored;
}
