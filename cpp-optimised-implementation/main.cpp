#include <NTL/ZZ.h>

#include "montgomery.hpp"
#include "fp2.hpp"
#include "isog.hpp"
#include "scallop.hpp"


void test_fp2_arith()
{
    std::cout << "Testing Fp2 arithmetic" << std::endl;
    NTL::ZZ p(103);
    NTL::ZZ_p::init(p);

    fp2_elem a{NTL::ZZ_p(12), NTL::ZZ_p(10)};
    fp2_elem b{NTL::ZZ_p(42), NTL::ZZ_p(3)};

    fp2_elem coeff1{NTL::ZZ_p(49), NTL::ZZ_p(90)};
    fp2_elem coeff2{NTL::ZZ_p(62), NTL::ZZ_p(44)};

    fp2_elem a_sqr{NTL::ZZ_p(44), NTL::ZZ_p(34)};
    fp2_elem apb{NTL::ZZ_p(54), NTL::ZZ_p(13)};
    fp2_elem amb{NTL::ZZ_p(73), NTL::ZZ_p(7)};
    fp2_elem ab_inv{NTL::ZZ_p(43), NTL::ZZ_p(83)};
    fp2_elem a_sqroot {NTL::ZZ_p(66), NTL::ZZ_p(11)};
    fp2_elem a_pow {NTL::ZZ_p(65), NTL::ZZ_p(39)};
    assert (Fp2_sqr(a) == a_sqr);
    assert (Fp2_add(a, b) == apb);
    assert (Fp2_sub(a, b) == amb);
    assert (Fp2_div(a, b) == ab_inv);
    assert (Fp2_sqroot(a) == a_sqroot);
    assert (Fp2_pow(a, 123) == a_pow);

    auto roots = QuadraticRoot(coeff1, coeff2);
    assert ((roots.first == a) || (roots.second == a));
    assert ((roots.first == b) || (roots.second == b));

    std::cout << "  Success!\n\n\n\n" << std::endl;
}

void test_ec_params_1024()
{
    std::cout << "Testing EC arithmetic" << std::endl;
    //prime
    NTL::ZZ p;
    NTL::conv(p, "15922486913903522274738876334491173327471544233065982291422062886364771917539984633619947369979657884209891416505022976305949903782081999018601968658399212422459232014682162496366184150793214330561748936802196212819401893052051149876176865225396801900263676014808365383967944863381855436796098378947797615849584416454584458385591707302057385868132351");
    NTL::ZZ_p::init(p);

    NTL::ZZ_p A_0, A_1;
    NTL::conv(A_0, "9441462816738741402523832537038431675206053350665219869672459087171062099050012989213043047665999028842373979097626520609601732380406921502704227010948378677950926261216650604354277728513533256911748160796963384290747601854926112699155094011097433046522801436239578500871867010454362628828348179469507729212239950281818700688002010280955291469956886");
    NTL::conv(A_1, "15393241905357780709372343104273255578766701851965068876871056094205683936179533072517368856268240772747198173948856812668384523024263411280876291237727157314158232913760703417823150763000262604626224666643383569371303859527030199684944864642126660891277353825688615660864404467715468043674393179676397294278005826873824264231827260686176208366387282");

    ProjA A{{A_0, A_1}, Fp2_one()};

    NTL::ZZ_p P_0, P_1, Q_0, Q_1, Qm_0, Qm_1;
    NTL::conv(P_0, "14353701636938821171268290034267352062667169181353268795745398786626800083112523144224085987276394943940114921776742915472447262366134254351465161855824290101604810851963442952666271126102715267831160381957691661718732833320323215498147280327312431699906161752162724892491384348021562186350185200344926940962876984099715078485709541696511025013871874");
    NTL::conv(P_1, "13692625063031617304072917414100364081604020586987836950101409447375090066391794517475456676670992812502879343162870450112555441794037987957420674396340461768063238307539820033779965081766507198354708822499721313102032148883606970184198503348910022722136637844345546957865325968904005063794577157912954791182986726863396668863043030093031129072418794");
    NTL::conv(Q_0, "13197046818621393400637866126002431448715786102131685694567310625416128489590739496840294126332244465917256862647739935484253127529788253353399219289314626834748076843292749404382106180682806429417436894803444395677399845356301491965460073482388139378276527852378069019511676347253106914935557840835364161531881770516756489360964336428573389894471789");
    NTL::conv(Q_1, "9379799174354683183851730564075159641730480528368169906513795773509842845874924428667197315012818042126078716799949923128810646468842452222818715283846506810948318691352938206138182228131619622789184730037674829264109024343971815756783711885709450128357574153996255066012329094766118369891196808930710533609852112792785092108194068537959716112817579");
    NTL::conv(Qm_0, "14568109670088169988079472350244773753636577569053389831389405728117335367891502846421722912090086280557092621268330247722181222534380117146288475543152531051522014081130207748295646461639268616357337465457904604239102056486702447692316403075944036228192766673095321149492083800554281377403772282317488908724480733512255348690520890302960042843102813");
    NTL::conv(Qm_1, "12210357416695094027617085183496924203240927260555106204376096250016174889538158262294631091074412637094083252317930797490225984737429293238922162015222970871576684443900698090696191475877680924008255248569031655696699033790863567803457961257862754696400377571366667025076539001507146252529612485150546375860287857264084025273978580570643748218911853");

    xPoint P{{P_0, P_1}, Fp2_one()};
    xPoint Q{{Q_0, Q_1}, Fp2_one()};
    xPoint Qm{{Qm_0, Qm_1}, Fp2_one()};

    //std::vector<int> strat2{256, 128, 64, 32, 17, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 32, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 64, 32, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 32, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 128, 64, 32, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 32, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 64, 32, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 32, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1};
    
    std::vector<IsogConsts> KerGens_1, KerGens_2, KerGens_3;
    ProjA Am = FourIsogChainPrecompute(P, A, strat, KerGens_1);
    std::cout << "j-invariant away P: " << StringFp2(jInvariant(Am)) << std::endl;

    ProjA Amm = FourIsogChainPrecompute(Q, A, strat, KerGens_2);
    std::cout << "j-invariant away Q: " << StringFp2(jInvariant(Amm)) << std::endl;

    auto ur = IsomorphismConstants(Amm, Am);
    xPoint Qmm = FourIsogChainEvaluate(Qm, KerGens_2);

    Qmm = IsomorphismEval(Qmm, ur);
    ProjA A_start = FourIsogChainPrecompute(Qmm, Am, strat, KerGens_3);

    std::cout << "j-invariant orig: " << StringFp2(jInvariant(A)) << std::endl;
    std::cout << "j-invariant recovered: " << StringFp2(jInvariant(A_start)) << std::endl;

    auto ur_home = IsomorphismConstants(A_start, A);

    std::cout << "Finding K and evaluating" << std::endl;

    xPoint K = PointOfOrderDividing(A, NTL::ZZ(47));

    xPoint Kprev = K;
    xPoint Ki = xDBL(K, A);
    xPoint temp;
    
    for (size_t j = 2; j < 47; j++) {
        temp = Ki;
        assert (!IsIdentity(Ki));
        Ki = xADD(Ki, K, Kprev);
        Kprev = temp;
    }
    assert (IsIdentity(Ki));

    xPoint wK = FourIsogChainEvaluate(K, KerGens_1);

    wK = FourIsogChainEvaluate(wK, KerGens_3);
    wK = IsomorphismEval(wK, ur_home);

    assert (IsIdentity(xMUL(wK, NTL::ZZ(47), A)));

    NTL::ZZ lampos(18);
    NTL::ZZ lamneg(21);
    xPoint lamK = xMUL(K, lamneg, A);
    assert (IsIdentity(xMUL(K, NTL::ZZ(47), A)));
    assert (IsIdentity(xMUL(lamK, NTL::ZZ(47), A)));

    auto PointDiff = xADDSUB(wK, K, A);

    xPoint Kp = Ladder3pt(wK, K, PointDiff.first, lamneg, A);
    assert (IsIdentity(xMUL(Kp, NTL::ZZ(47), A)));

    // Test that we got the correct sign
    xPoint wKp = FourIsogChainEvaluate(Kp, KerGens_1);
    wKp = FourIsogChainEvaluate(wKp, KerGens_3);
    wKp = IsomorphismEval(wKp, ur_home);

    xPoint lamKp = xMUL(Kp, lampos, A);

    xPoint Kn;

    // Compute Kp and Kn accordingly
    if (!(PointEqual(lamKp, wKp))) {
        std::cout << "SIGN FLIP" << std::endl;
        Kp = Ladder3pt(wK, K, PointDiff.second, lamneg, A);
        assert (IsIdentity(xMUL(Kp, NTL::ZZ(47), A)));
    }
    NormalizePoint(Kp);

    std::vector<xPoint> evalPts{P, Q, Qm};
    ProjA A2 = xISOG(Kp, A, 47, evalPts);

    std::vector<IsogConsts> KerGens_4, KerGens_5;
    ProjA A2m = FourIsogChainPrecompute(evalPts[0], A2, strat, KerGens_4);
    std::cout << "j-invariant away P2: " << StringFp2(jInvariant(A2m)) << std::endl;

    ProjA A2mm = FourIsogChainPrecompute(evalPts[1], A2, strat, KerGens_5);
    std::cout << "j-invariant away Q2: " << StringFp2(jInvariant(A2mm)) << std::endl;
    assert (jInvariant(A2mm) == jInvariant(A2m));

    std::cout << "  Success!\n\n\n\n" << std::endl;
}

void time_four_isog_mult()
{
    std::cout << "Timing four isog vs double double" << std::endl;
    NTL::ZZ_p::init(p);

    NTL::ZZ_p A_0, A_1;
    NTL::conv(A_0, A_re);
    NTL::conv(A_1, A_im);

    ProjA A{{A_0, A_1}, Fp2_one()};

    NTL::ZZ_p P_0, P_1, Q_0, Q_1, Qm_0, Qm_1;
    NTL::conv(P_0, P_re);
    NTL::conv(P_1, P_im);
    NTL::conv(Q_0, Q_re);
    NTL::conv(Q_1, Q_im);
    NTL::conv(Qm_0, Qm_re);
    NTL::conv(Qm_1, Qm_im);

    xPoint P{{P_0, P_1}, Fp2_one()};
    xPoint Q{{Q_0, Q_1}, Fp2_one()};
    xPoint Qm{{Qm_0, Qm_1}, Fp2_one()};

    
    int rep = 1000;
    auto start = std::chrono::steady_clock::now();
    for (int i = 0; i < rep; i++) {
        xDBL(xDBL(P, A), A);
    }
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start);
    std::cout << "Multiplication by 2 took: " << duration.count() << " microseconds" << std::endl;
    xPoint Ker = xDBLe(P, e-2, A);
    assert (!(IsIdentity(xDBL(Ker, A))));
    assert (IsIdentity(xDBLe(Ker, 2, A)));

    IsogConsts Ks;
    ProjA A2 = FourIsogCurve(Ker, Ks);

    start = std::chrono::steady_clock::now();
    for (int i = 0; i < rep; i++) {
        FourIsogEval(Q, Ks);
    }
    duration = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start);
    std::cout << "Four isog eval took: " << duration.count() << " microseconds" << std::endl;
}

void test_scallop()
{
    std::cout << "Testing PEARL-SCALLOP" << std::endl;
    NTL::ZZ_p::init(p);

    NTL::ZZ_p A_0, A_1;
    NTL::conv(A_0, A_re);
    NTL::conv(A_1, A_im);

    ProjA A{{A_0, A_1}, Fp2_one()};

    NTL::ZZ_p P_0, P_1, Q_0, Q_1, Qm_0, Qm_1;
    NTL::conv(P_0, P_re);
    NTL::conv(P_1, P_im);
    NTL::conv(Q_0, Q_re);
    NTL::conv(Q_1, Q_im);
    NTL::conv(Qm_0, Qm_re);
    NTL::conv(Qm_1, Qm_im);

    xPoint P{{P_0, P_1}, Fp2_one()};
    xPoint Q{{Q_0, Q_1}, Fp2_one()};
    xPoint Qm{{Qm_0, Qm_1}, Fp2_one()};

    // Just a random L-infty norm ~(-20, 20) vector
    // std::vector<int> es = GenSecret(dim, N_ball);
    // auto PK_A = GroupAction(P, Q, Qm, A, es);

    // Key Exchange
    
    /*std::vector<int> es_A = GenSecret(dim, 3);
    std::vector<int> es_B = GenSecret(dim, 3);

    xPoint P_A = P;
    xPoint P_B = P;
    xPoint Q_A = Q;
    xPoint Q_B = Q;
    xPoint Qm_A = Qm;
    xPoint Qm_B = Qm;
    ProjA PK_A = GroupAction(P_A, Q_A, Qm_A, A, es_A);
    std::cout << "\n\n~~~~~~~Alice public key done~~~~~~~~\n\n\n" << std::endl;
    ProjA PK_B = GroupAction(P_B, Q_B, Qm_B, A, es_B);
    std::cout << "\n\n~~~~~~~Bob public key done~~~~~~~~\n\n\n" << std::endl;

    auto SS_A = GroupAction(P_B, Q_B, Qm_B, PK_B, es_A);
    std::cout << "\n\n~~~~~~~Alice shared key done~~~~~~~~\n\n\n" << std::endl;

    auto SS_B = GroupAction(P_A, Q_A, Qm_A, PK_A, es_B);
    std::cout << "\n\n~~~~~~~Bob shared key done~~~~~~~~\n\n\n" << std::endl;

    std::cout << "Shared keys equal??" << std::endl;
    std::cout << "j(E_AB) = " << StringFp2(jInvariant(SS_A)) << std::endl;
    std::cout << "j(E_AB) = " << StringFp2(jInvariant(SS_B)) << std::endl;
    assert (Fp2_equal(jInvariant(SS_A), jInvariant(SS_B))); */

    
    //Reduced vector 1024?
    std::vector<int> es{0,-3,5,-2,-2,-1,1,16,14,6,5,-17,16,27,8,-34,6,9,1,2,19,-24,21,35,-2,41,-11,-5,60,-11,80,6,20,13,15,8,22,2,-21,-12,7,-19,-68,-39,9,-68,-13,33,-1,20,-31,-104,-18,-23,6,30,11,11,-7,-11,9,9,8,33,-3,12,0,4,8,4,-4,-11,-1,0,11};
    auto PK_A = GroupAction(P, Q, Qm, A, es);

    std::cout << "  Success!\n\n\n\n" << std::endl;
}

int main()
{
    //test_fp2_arith();
    //test_ec_params_1024();
    //time_four_isog_mult();
    test_scallop();
    return 0;
}