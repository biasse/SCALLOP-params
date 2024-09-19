def pullback(J, I):
    assert J.right_order() == I.left_order()
    return J*I + J.left_order()*I.norm()

def IdealGenerator(I):
    p = I.quaternion_algebra().ramified_primes()[0]
    bound = max(ceil(10*log(p,2)), 100)
    while True:
        alpha = sum([b * randint(-bound, bound) for b in I.basis()])
        if gcd(alpha.reduced_norm(), I.norm()**2) == I.norm():
            return alpha
        
def to_NTL(alpha):
    d = alpha.denominator()
    coeffs = list((alpha*d).coefficient_tuple())
    coeffs.append(d)
    s = ""
    for c in coeffs:
        s += f"NTL::conv<NTL::ZZ>(\"{c}\"), "
    s = s[:-2]
    print(s)

def prime_to_torsion(N, p):
    primeToTorsion = {}
    ks = {}
    for l, e in factor(N):
        for k in range(1,100):
            if (l**e).divides(p**(k) - (-1)**k):
                primeToTorsion[l**e] = k
                if not k in ks.keys():
                    ks[k] = []
                ks[k].append((l,e))
                break
                
    return ks, primeToTorsion

def ReducedBasis(I):
    def _matrix_to_gens(M, B):
        return [sum(c * g for c, g in zip(row, B)) for row in M]
    
    def gram_matrix(basis):
        M = []
        for a in basis:
            M.append([QQ(2) * a.pair(b) for b in basis])
        return Matrix(QQ, M)

    B = I.basis()
    G = gram_matrix(I.basis())
    U = G.LLL_gram().transpose()
    reduced_basis_elements = _matrix_to_gens(U, B)
    return reduced_basis_elements

def connectingIdeal(O1, O2):
    I = O1*O2
    I *= I.norm().denominator()
    return I

def isIsomorphic(O1, O2):
    I = connectingIdeal(O1, O2)
    for alpha in ReducedBasis(I):
        if alpha.reduced_norm() == I.norm():
            return True
    return False

#512
#p = 167286323857221689112346016933207258999493176647479781908348180838625562682489086996433613891517156716513168813389511283523303877305870043625319556074092350622078682027926127121250390130545025689333161800187261717666012808114506310888652381411304126074312458239

# 1536
p = 23771452727449498898004557547225354988120021364701398066313318892838230235651289792849929330578203612228877729959916570183072595360843544888617115522542149031138111127294471888062615222219308670025415410876727090874823170475481539518868548152615795234974491096957906735993259268864790983747139027153544490397993988855511337575848418797409297084369809840481360670399419714736596634148006572830334348569025487101216814643806207

B = QuaternionAlgebra(-1, -p)
i, j, k = B.gens()
O0 = B.quaternion_order([1, i, (i+j)/2, (1+k)/2])

#512
#I = O0.left_ideal((Integer(1)/Integer(2) + Integer(8028170362878520627000397624849124892182804416632981394797083245489770480253778314941516911176802745559682787465824306790119665950954458423109837617464626788671289332024103172895239588435040315)*j + Integer(5517406830432160569531460434783887168508512136923704264395014092668084428703376490606065995099857645960283681480115728322350240403716673861139024369914489493617824498398267835375843424885744371)/Integer(2)*k, Integer(1)/Integer(2)*i + Integer(12559263776757650835333586527979651668129931543141805391191597038825242054633959978712546746135088882968796071885052394436683145813123728240242016740549303638376323965584247860375004143396377679)/Integer(2)*j + Integer(8028170362878520627000397624849124892182804416632981394797083245489770480253778314941516911176802745559682787465824306790119665950954458423109837617464626788671289332024103172895239588435040315)*k, Integer(9038335303594905702432523481381769418319221840032754827793305565746663241668668234659306370617473264464539876682584061379516693108420201050690520555231896565997074231991257847875423784141061025)*j, Integer(9038335303594905702432523481381769418319221840032754827793305565746663241668668234659306370617473264464539876682584061379516693108420201050690520555231896565997074231991257847875423784141061025)*k))  
#e = 256
#delta = Integer(159552265354474160974185212535768133044337794207791175162387613622343972203672892620294385823603408615196894574815513736838239809843521949481306153521731911016128471720392737067923933657432024617399001841668116552111895154362008623304629452031945255102495907218)/Integer(9038335303594905702432523481381769418319221840032754827793305565746663241668668234659306370617473264464539876682584061379516693108420201050690520555231896565997074231991257847875423784141061025)*i + Integer(4902147998203463609140318443449451568899667218603278364071763525977316879466888129190741568857630859366516877631404224168670484233439)/Integer(549808088545802169163220059998624582926252705370476475637848044457656659938712373541945240502154425889132548817234649507275543645826026872116847498349322773005331519118581555533905335425)*j + Integer(1783766779885671839023099111918400074849355529222311908852243800862128393632296029189893648014519111591965727974458431623743887874422)/Integer(2265448505717903648678647719629633568480208600021494350448048777709962475776333040824260396756471632053613761796496696167081706043625485738879153547554205089925814009013648544638789612645)*k
#omega = Integer(2310820388073888658175157102723307102542652012212569935585880870910116436225) - delta
    
#1536
e = 762
I = O0.left_ideal((1/2 + 148517744027306236980775731157810778948711920318943987807540783716230652454126991998817343738733811348148028698509926219317522751190441122264308788716483931431089914614568165641455111020608050592475658*j + 242984351222788300942193050348357225225870083836269207547058002966214500516954773806743967376129739548314897225693549140777334730285172725935880441211266193130076180906476532964360808915185387420641049/2*k, 1/2*i + 149693651364335684515151031299529758203083734966198966244227867833479591615789768264757395647222339923105619847706762233064649002805474938505778219003657215579613028579153303678225016696525221880174401/2*j + 148517744027306236980775731157810778948711920318943987807540783716230652454126991998817343738733811348148028698509926219317522751190441122264308788716483931431089914614568165641455111020608050592475658*k, 196339001293561992728672040823943491714476909401234086895642935399847046066372271035750681511676039735710258536700155686920991866545323832220829330107461704354844604742814918321292912805855304650407725*j, 196339001293561992728672040823943491714476909401234086895642935399847046066372271035750681511676039735710258536700155686920991866545323832220829330107461704354844604742814918321292912805855304650407725*k))
delta = 14427982307211528997591287295588710162146728674124622434118957127836601304769368466350591244746945199823588548991363387194649728996822181428250059536920366422372008367243557589187016744957373254413949462748533591537217340629643172951166628342814073380741948755501024813793495563495792122922278303094317310361348642732139906970542461670846471681786737481329496753902314796473890792809882896404064236867911428232295540692783849/196339001293561992728672040823943491714476909401234086895642935399847046066372271035750681511676039735710258536700155686920991866545323832220829330107461704354844604742814918321292912805855304650407725*i + 2254876261945793443139191746292308322213422382177367390294242508815183654231117779707515164552068236937773894415177101886243178080045515938514448813615381457252990221726161615992642516699290500504458274104731464005232/28048428756223141818381720117706213102068129914462012413663276485692435152338895862250097358810862819387179790957165098131570266649331976031547047158208814907834943534687845474470416115122186378629675*j + 10693684360183525167427939253204907990223517408624315170265841172302140115949058443940299111135329250502341254569657034800580906771795725633084386094694018403795213753935992499445227812171785787545048030611925378709609/196339001293561992728672040823943491714476909401234086895642935399847046066372271035750681511676039735710258536700155686920991866545323832220829330107461704354844604742814918321292912805855304650407725*k
omega = 19086769594710305957479658851316502070038121206500346144568806255427709054642899899552063491995584993331593399986913140812754447230564684801247389937557608054790951503629749489741455198170946951717020050611071188339934920597600736 - delta

O_start = I.right_order()
print(f"omega norm: {omega.reduced_norm()}")
print(f"factored: {factor(omega.reduced_norm())}")
assert delta in O_start
assert delta/2 not in O_start
assert omega in O_start
assert omega/2 not in O_start

I_omega_half1 = O_start*omega + O_start*2**e
I_omega_half2 = O_start*omega.conjugate() + O_start*2**e
I_o1 = pullback(I, I_omega_half1)
I_o2 = pullback(I, I_omega_half2)

alpha = IdealGenerator(I)

print("this should be true:")
print(isIsomorphic(I_omega_half1.right_order(), I_omega_half2.right_order()))

print("\n----- I_start -----")
print(f"I_start norm: {factor(I.norm())}")
to_NTL(alpha)

print("\n----- I_o1 -----")
print(f"I_o1 norm: {factor(I_o1.norm())}")
print(I_o1)
to_NTL(IdealGenerator(I_o1))

print("\n----- I_o2 -----")
print(f"I_o2 norm: {factor(I_o2.norm())}")
print(I_o2)
to_NTL(IdealGenerator(I_o2))

ks, tors = prime_to_torsion(I.norm(), p)
klist = ks.keys()

s = "std::array<unsigned, " + str(len(klist)) + "> ks {"

for k in sorted(klist):
    s += f"{k}, "

s = s[:-2]
s += "};"

print(s)

for k in sorted(ks.keys()):
    row = ks[k]
    s = f'extension k = {k}:  '
    for l, e in row:
        if e > 1:
            s += f'{l}^{e}, '
        else:
            s += f'{l}, '
    print(s)

print(factor(I.norm()))
print(factor(p+1))