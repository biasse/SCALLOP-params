#pragma once

#include <map>

#include <NTL/ZZ_pE.h>

#include "quaternions.hpp"
#include "isog.hpp"
#include "utils.hpp"
#include "Fp2k.hpp"

int torsionToFieldDegree(Integer ell_e);

std::vector<std::pair<ecp,std::pair<int, int>>> idealToKernel(quat const &alpha, Integer const &N, ec const &E0, std::map<unsigned,Fp2k> &FieldExtensions, std::map<NTL::ZZ,std::pair<ecp,ecp>> &TorsionBases);
isog_chain idealToIsogeny(quat const &alpha, Integer const &N, ec const &E0, std::map<unsigned,Fp2k> &FieldExtensions, std::map<NTL::ZZ,std::pair<ecp,ecp>> &TorsionBases);
