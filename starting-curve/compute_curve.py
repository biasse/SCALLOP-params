from klpt import *
from id2iso import *
from param_loader import load

if __name__=="__main__":
    params = load('75-primes')
    p = params['p']
    F = params['F']
    f = params['f']
    D_com = params['D_com']
    D_chall = params['D_chall']
    T = params['T']
    B_2 = params['B_2']
    B_chall = params['B_chall']
    facToExt = params['facToExt']
    facToBasis = params['facToBasis']
    facToAction = params['facToAction']
    small_ns = params['small_ns']
    small_s = params['small_s']

    B = QuaternionAlgebra(-1, -p)
    i, j, k = B.gens()
    O0 = B.quaternion_order([1, i, (i+j)/2, (1+k)/2])
    I = O0.left_ideal((Integer(1)/Integer(2) + Integer(685007987407002328385657243320103317114507692539286019973979404696562077561845287517851438169205921643128283263047282449676994580335515837286045935409445563266189766896270322283618553453785413503026053008941)*j + Integer(339622259194838739769576394104774179369222497042180578849092505553706635091239698632368125210224223605737409434153645788235545862572169503978401359957725222863270532615665104120442395699141646241989356438805)/Integer(2)*k, Integer(1)/Integer(2)*i + Integer(1184726458466882292302821208958681168242462210075966352803615602016024949942781072437832740818306663546880857467395318464628923471991906839414666112879696222726971921532529517869735605685099989623317652873545)/Integer(2)*j + Integer(685007987407002328385657243320103317114507692539286019973979404696562077561845287517851438169205921643128283263047282449676994580335515837286045935409445563266189766896270322283618553453785413503026053008941)*k, Integer(762174358830860516036198801531727673805842353559073465826354053784865792517010385535100433014265443576309133450774482126432234667282038171696533736418710722795121227074097310995089000692120817932653504656175)*j, Integer(762174358830860516036198801531727673805842353559073465826354053784865792517010385535100433014265443576309133450774482126432234667282038171696533736418710722795121227074097310995089000692120817932653504656175)*k))
        
        
    O_start = I.right_order()
    delta = Integer(9606326337921206847938752593969970419136633589855810918422788159488757873099498049930964921611333196319529099574454436999312731782500515928142614952680748607954733962712409752928747365200187621269022899355581229637730676541138884421891907599597158372625349366357542367555115125368877715895311154324651689662847557618090043457299911794953302842862481)/Integer(762174358830860516036198801531727673805842353559073465826354053784865792517010385535100433014265443576309133450774482126432234667282038171696533736418710722795121227074097310995089000692120817932653504656175)*i + Integer(31532820412484549496770528269636797662904276836008793921008597228152532026126859825380363444457778970576206970784428987920681538210435390632681087323816696196259223294668683752212215939324)/Integer(152434871766172103207239760306345534761168470711814693165270810756973158503402077107020086602853088715261826690154896425286446933456407634339306747283742144559024245414819462199017800138424163586530700931235)*j - Integer(39830127181740430214048082814701718232352105340024973342555319105145974307097262053105000144370778586334286452924022230408579592869420623175179376180337798524873388783350137916426516953783)/Integer(762174358830860516036198801531727673805842353559073465826354053784865792517010385535100433014265443576309133450774482126432234667282038171696533736418710722795121227074097310995089000692120817932653504656175)*k
    assert delta in O_start
    print(delta/2 in O_start)
    omega = Integer(107262463439540776796592199985616134275233279428307433766342950569669432214610248132438916680001426895938208196632896833496294702258913582339408834412325768) - delta
    print(f"omega norm: {omega.reduced_norm()}")
    print(f"factored: {factor(omega.reduced_norm())}")
    assert omega in O_start
    print(omega/2 in O_start)
    I_omega_half1 = O_start*omega + O_start*2**518
    I_omega_half2 = O_start*omega.conjugate() + O_start*2**518

    print(factor(I_omega_half1.norm()))
    print(factor(I_omega_half2.norm()))
    E0 = EllipticCurve(F, [1,0])
    phi_start = IdealToIsogeny(O0, I, E0, facToBasis, facToAction)

    E_start = phi_start.codomain()

    K_1 = IdealToIsogenyGens(O0, pullback(I, I_omega_half1), facToBasis, facToAction)
    K_2 = IdealToIsogenyGens(O0, pullback(I, I_omega_half2), facToBasis, facToAction)


    print("Got the gens on E_0")
    print(K_1, K_2)

    K_1 = E0.lift_x(K_1[0][0].X)
    K_2 = E0.lift_x(K_2[0][0].X)
    K_1 = phi_start(K_1)
    K_2 = phi_start(K_2)

    print("All done!")
    print("Starting curve:")
    print(E_start)
    print("Effective orientation:")
    print(K_1) 
    print(K_2)