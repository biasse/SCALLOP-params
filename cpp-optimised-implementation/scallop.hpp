#pragma once

#include <random>

#include "fp2.hpp"
#include "montgomery.hpp"
#include "isog.hpp"

//#include "params_512.hpp"
//#include "params_1024.hpp"
#include "params_1536.hpp"


void printVector(std::vector<int> const &vec) {
    std::cout << "es: " << std::endl;
    for (int elem : vec) {
        std::cout << elem << ", ";
    }
    std::cout << std::endl;
}

ProjA ActionIdealStep(xPoint &P, xPoint &Q, xPoint &Qm, ProjA A, NTL::ZZ Lpos, NTL::ZZ Lneg, std::vector<int> &es, std::chrono::duration<long, std::milli> &d1, std::chrono::duration<long, std::milli> &d2)
{
    auto start_d1 = std::chrono::steady_clock::now();
    NormalizeCoeff(A);

    NTL::ZZ L = Lpos*Lneg;

    std::cout << "K" << std::endl;
    xPoint K = PointOfOrderDividing(A, L);

    while (IsIdentity(K)) {
        K = PointOfOrderDividing(A, L);
    }

    //Computing omega
    std::cout << "Computing kernel generator.... " << std::endl;
    
    std::vector<IsogConsts> KerGens_1, KerGens_2, KerGens_3;
    ProjA Am = FourIsogChainPrecompute(P, A, strat, KerGens_1);

    ProjA Amm = FourIsogChainPrecompute(Q, A, strat, KerGens_2);

    auto ur = IsomorphismConstants(Amm, Am);
    xPoint Qmm = FourIsogChainEvaluate(Qm, KerGens_2);

    Qmm = IsomorphismEval(Qmm, ur);
    ProjA A_start = FourIsogChainPrecompute(Qmm, Am, strat, KerGens_3);

    auto ur_home = IsomorphismConstants(A_start, A);

    xPoint wK = FourIsogChainEvaluate(K, KerGens_1);
    wK = FourIsogChainEvaluate(wK, KerGens_3);
    wK = IsomorphismEval(wK, ur_home);

    auto PointDiff = xADDSUB(wK, K, A);

    xPoint Kp;
    xPoint Kn;
    if (Lpos > 1) {
        Kp = Ladder3pt(wK, K, PointDiff.first, (lamneg % Lpos), A);
        Kp = xMUL(Kp, Lneg, A);

        // Test that we got the correct sign
        xPoint wKp = FourIsogChainEvaluate(Kp, KerGens_1);
        wKp = FourIsogChainEvaluate(wKp, KerGens_3);
        wKp = IsomorphismEval(wKp, ur_home);
        assert (IsIdentity(xMUL(wKp, Lpos, A)));
        xPoint lamKp = xMUL(Kp, lampos, A);
        assert (IsIdentity(xMUL(lamKp, Lpos, A)));

        // Compute Kp and Kn accordingly
        if ((PointEqual(lamKp, wKp))) {
            Kn = Ladder3pt(wK, K, PointDiff.first, (lampos % Lneg), A);
            Kn = xMUL(Kn, Lpos, A);
        } else {
            std::cout << "SIGN FLIP THING" << std::endl;
            Kp = Ladder3pt(wK, K, PointDiff.second, (lamneg % Lpos), A);
            Kp = xMUL(Kp, Lneg, A);
            Kn = Ladder3pt(wK, K, PointDiff.second, (lampos % Lneg), A);
            Kn = xMUL(Kn, Lpos, A);
        }
    } else {        
        Kn = Ladder3pt(wK, K, PointDiff.first, (lampos % Lneg), A);

        // Test that we got the correct sign
        xPoint wKn = FourIsogChainEvaluate(Kn, KerGens_1);
        wKn = FourIsogChainEvaluate(wKn, KerGens_3);
        wKn = IsomorphismEval(wKn, ur_home);
        assert (IsIdentity(xMUL(wKn, Lneg, A)));
        xPoint lamKn = xMUL(Kn, lamneg, A);
        assert (IsIdentity(xMUL(lamKn, Lneg, A)));

        // Compute Kp and Kn accordingly
        if (!(PointEqual(lamKn, wKn))) {
            std::cout << "SIGN FLIP THING" << std::endl;
            Kn = Ladder3pt(wK, K, PointDiff.second, (lampos % Lneg), A);
        }

    }

    d1 = std::chrono::duration_cast<std::chrono::milliseconds>(d1 + std::chrono::steady_clock::now() - start_d1);
    auto start_d2 = std::chrono::steady_clock::now();
    ProjA Ai = A;

    //Compute isogeny, and update exponent vector
    std::cout << "Computing odd isogenies..." << std::endl;

    std::vector<xPoint> evalPts;
    if (Lneg == 1) {
        evalPts = {P, Q, Qm};
        Ai = xISOG_Chain(Kp, Ai, Lpos, evalPts, 1, ells, es);
    } else if (Lpos == 1) {
        evalPts = {P, Q, Qm};
        Ai = xISOG_Chain(Kn, Ai, Lneg, evalPts, -1, ells, es);
    } else {
        evalPts = {P, Q, Qm, Kn};
        Ai = xISOG_Chain(Kp, Ai, Lpos, evalPts, 1, ells, es);
        Kn = evalPts[3];
        evalPts.pop_back();
        Ai = xISOG_Chain(Kn, Ai, Lneg, evalPts, -1, ells, es);
    }

    P = evalPts[0];
    Q = evalPts[1];
    Qm = evalPts[2];

    d2 = std::chrono::duration_cast<std::chrono::milliseconds>(d2 + std::chrono::steady_clock::now() - start_d2);
    return Ai;
}

bool _all_zero(std::vector<int> const &vec) {
    for (int elem : vec) {
        if (elem != 0) {
            return false;
        }
    }
    return true;
}

ProjA GroupAction(xPoint &P, xPoint &Q, xPoint &Qm, ProjA const &A, std::vector<int> const &es_in)
{
    auto start = std::chrono::steady_clock::now();
    std::chrono::duration<long, std::milli> duration_compute_kernel = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start);
    std::chrono::duration<long, std::milli> duration_odd_isogenies = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start);
    std::vector<int> es(es_in);
    
    ProjA Ai = A;

    while (!(_all_zero(es))) {
        printVector(es);
        NTL::ZZ Lpos(1);
        NTL::ZZ Lneg(1);

        for (size_t i = 0 ; i < ells.size() ; i++) {
            if (es[i] > 0) {
                Lpos *= ells[i];
            }
            if (es[i] < 0) {
                Lneg *= ells[i];
            }
        }
        Ai = ActionIdealStep(P, Q, Qm, Ai, Lpos, Lneg, es, duration_compute_kernel, duration_odd_isogenies);

        auto duration = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - start);
        std::cout << ">>> We have used " << duration.count() << " seconds" << std::endl;
        std::cout << ">>> On finding ker gens: " << duration_compute_kernel.count() << " milliseconds" << std::endl;
        std::cout << ">>> On computing odd isos: " << duration_odd_isogenies.count() << " milliseconds" << std::endl;
    }
    return Ai;
}

std::vector<int> GenSecret(size_t length, int N) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dis(-N, N);

    std::vector<int> vec(length);
    for (int i = 0; i < length; ++i) {
        vec[i] = dis(gen);
    }
    return vec;
}