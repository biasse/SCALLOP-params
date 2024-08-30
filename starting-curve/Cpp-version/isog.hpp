#pragma once

#include <cassert>
#include <optional>
#include <vector>

#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pEX.h>
#include "Fp2k.hpp"
#include "ec.hpp"
#include "ecp.hpp"

class isog
{
private:
    std::shared_ptr<const ec> domain;
    std::shared_ptr<const ec> codomain;
    int _degree;
    NTL::ZZ_pEX psi, psi_der, psi_der_der, psi_der_der_der;
    NTL::ZZ_pE t; // Useful for degree 2 isogenies

public:
    isog(std::shared_ptr<const ec> E0, std::shared_ptr<const ec> E1, int degree, NTL::ZZ_pEX const &psi, NTL::ZZ_pEX const &psi_der, NTL::ZZ_pEX const &psi_der_der, NTL::ZZ_pEX const &psi_der_der_der) :
        domain(E0),
        codomain(E1),
        _degree(degree),
        psi(psi),
        psi_der(psi_der),
        psi_der_der(psi_der_der),
        psi_der_der_der(psi_der_der_der) {};

    isog(std::shared_ptr<const ec> E0, std::shared_ptr<const ec> E1, int degree, NTL::ZZ_pEX const &psi, NTL::ZZ_pE const &t) :
        domain(E0),
        codomain(E1),
        _degree(degree),
        psi(psi),
        t(t) {};

    ec const &get_domain() const { return *domain; }
    std::shared_ptr<const ec> const &get_domain_ptr() const { return domain; }
    ec const &get_codomain() const { return *codomain; }
    std::shared_ptr<const ec> const &get_codomain_ptr() const { return codomain; }

    ecp operator()(ecp const &P) const;
    ecp _even_evaluation(ecp const &P) const;
};

isog isogeny(ecp const &K, int degree);

NTL::ZZ_pEX kernel_polynomial(ecp const &K, int degree);

class isog_chain
{
    private:
        std::vector<isog> isog_parts;
    
    public:
        isog_chain(std::vector<std::pair<ecp,std::pair<int, int>>> kerGens);
        ecp operator()(ecp const &P) const {
            ecp Pi = P;
            for (isog phi : isog_parts) {
                Pi = phi(Pi);
            }
            return Pi;
        };

        ec const &get_domain() const { return isog_parts.front().get_domain(); }
        ec const &get_codomain() const { return isog_parts.back().get_codomain(); }
};