from sage.all import *
from utilities import BiDLP
from ec import TorsionBasis
from xonly import xPoint

####################################
#                                  #
#     cost.py from DeuringFTP      #
#                                  #
####################################

def smoothPart(n, B):
    r"""
    Compute the maximal odd divisor of n which is B-smooth.
    """
    n = ZZ(n)
    n = n.prime_to_m_part(2)
    return prod(l**e for l,e in n.factor(limit=B) if l < B)

def cost_model(p, newconst=None):
    r"""
    Given a prime characteristic p, return a function in two arguments
    l^e and k which roughly describes the cost of the Deuring algorithm
    when using l^e-torsion in GF(p^(2k)).
    """
    logp = p.log(2).n()
    loglogp = logp.log(2).n()
    log3 = log(3,2).n()

    # experimentally guessed using Sage 9.7 on a 2016 laptop
    if logp < 32:
        cutoff = 79
        def karatsuba(k):
            return k**log3 * (1/30 + logp/1e3 + logp**2/5e5 + (logp>=30)/50 + (logp>=31)/30)
    elif logp <= 64:
        cutoff = 44
        def karatsuba(k):
            return k**log3 * (1/30 + logp/190)
    else:
        cutoff = 55
        def karatsuba(k):
            return k**log3 * (1/30 + loglogp/50)

    def quasilinear(k):
        fun0 = lambda k: k*log(k,2).n() * (1/10 + logp/200)
        off = karatsuba(cutoff) - fun0(cutoff)
        return off + fun0(k)

    def field_model(kk):
        if kk < cutoff:
            oneop = karatsuba(RR(kk))
        elif kk > 100:
            oneop = kk**3 #Max limit for field size
        else:
            oneop = quasilinear(RR(kk))
        return RR(oneop)

    if newconst:
        c1, c2, c3, c4 = newconst
    else:
        c1, c2, c3, c4 = (0.01, 2.65, 1.25, 0.022)  # Empirically guesstimated

    def model(le, k):
        le = ZZ(le)
        (l,e), = le.factor()
        logl = RR(l).log(2)
        oneop = field_model(k)
        # Roughly some additions + isogeny + evaluating. For additions contain e**2 due to isogeny computations not using optimal strategies (reduces to e log e).
        return RR(e*c2*e*oneop*logl + c3*e*l*(k+c4*logl**2) + c1*oneop*l*logl) #(the size of the biggest isogeny to evaluate the point in can be uppper bounded by the current \ell by sorting them)

    return model

def adj_cost(cost, ell):
    r"""
    Adjusted cost to greedily solve the approximate subset sum
    """
    return round(cost/max(log(ell, 2), 1), 2)

def choose_torsion(p, q, lowbound, newconst=None):
    r"""
    Given p,q,N as in the Deuring correspondence, greedily find a set
    of pairs (l^e,k) for the specific characteristic which minimize
    the cost of the Deuring algorithm according to the cost_model().
    """
    facToExt = dict()

    # establish a baseline: take every l
    print(lowbound)
    le = ZZ.one()
    while lcm(le for le in facToExt.keys()) < lowbound*2*q:
        le += 1
        if p.divides(le):
            continue
        if le % 2 == 0:
            continue
        if not is_pseudoprime(le.radical()):
            continue
        k = Mod(p, le).multiplicative_order()*2
        #if k%2 == 0 and pow(p, k//2, le) - le == -1: # Use twist in this case (cant just divide k by 2, since (ZZ/2^eZZ)^* is not cyclic...)
        #    k //= 2
        facToExt[le] = k

    model = cost_model(p, newconst=newconst)

    # now, keep increasing k, looking for small-ish l defined over small k
    k = ZZ.zero()
    while True:
        # sort pairs (l^e,k) by cost estimate
        tups = sorted(facToExt.items(), key = lambda tup: adj_cost(model(*tup), tup[0].prime_factors()[0]))
        # compute T to check what's the worst (l^e,k) pair we currently have to use
        it = 0
        T = ZZ.one()
        while it < len(tups) and T//T.gcd(2*q) < lowbound:
            tup = tups[it]
            cost = adj_cost(model(*tup), tup[0].prime_factors()[0])
            T = lcm(T, tup[0])
            it += 1

        facToExt = dict(tups[:it])
        T //= T.gcd(2 * q)   # point divisions
        assert T >= lowbound

        k += 1
        # figure out up to which prime l it's worth searching for this k
        maxlebits = 0
        while adj_cost(model(next_prime(1 << maxlebits), k), next_prime(1 << maxlebits)) < cost:
            maxlebits += 1
        maxle = 0
        for i in reversed(range(maxlebits+1)):
            if adj_cost(model(next_prime(maxle | (1 << i)), k), next_prime(maxle | (1 << i))) < cost:
                maxle |= 1 << i
        maxle = min(maxle, 100000)
        print(maxle, k)
        if maxle < 1 or k > 100:   # no l^e at all -> done (and we're never gonna do k > 100)
            break

        # trial-divide the order in degree k to find new l^e
        on_curve = smoothPart(p**k - (-1)**k, maxle)
        #on_twist = smoothPart(p**k + (-1)**k, maxle)
        #for fac in (on_curve, on_twist):
        for l,e in on_curve.factor():
            for i in range(1,e+1):
                if l**i in facToExt:
                    facToExt[l**i] = min(k, facToExt[l**i])
                else:
                    facToExt[l**i] = k

    facToExt = dict(sorted(facToExt.items(), key = lambda tup: tup[0].prime_factors()))
    assert T >= lowbound

    def shave_excess(T, facToExt, B):
        tups = sorted([(le, k) for le, k in facToExt.items() if T % le == 0], key = lambda tup: adj_cost(model(*tup), tup[0].prime_factors()[0]))
        for le, k in tups:
            l, e = factor(le)[0]
            if T//l > B:
                return shave_excess(T//l, facToExt, B)
        return T
    
    T = shave_excess(T, facToExt, lowbound)
    return T, facToExt



####################################
#                                  #
#    Action matrices on basis      #
#                                  #
####################################

def endo_i(P):
    E0 = P.curve()
    F = E0.base_field()
    x, y = P.xy()
    return E0(-x, F(sqrtm1)*y)

def endo_j(P):
    E0 = P.curve()
    pi = E0.base_field().frobenius_endomorphism()
    x,y = P.xy()
    return E0(pi(x), pi(y))

def EvalEndomorphism(alpha, P, ord):
    d = lcm(c.denominator() for c in alpha)
    E = P.curve()
    assert P*ord == 0
    assert d == 1 or d == 2
    if gcd(d, ord) == 2:
        alpha = d*alpha
        Fbig, _ = E.base_field().extension(4,'A').objgen()
        Ebig = E.base_extend(Fbig)
        P = Ebig(P).division_points(2)[0]
    iP = endo_i(P)
    jP = endo_j(P)
    kP = endo_i(jP)
    coeffs = [coeff % ord for coeff in alpha]
    return E(sum(c*Q for c, Q in zip(coeffs, [P, iP, jP, kP])))

def ActionMatrix(alpha, basis, ord):
    P, Q = basis
    alphaP = EvalEndomorphism(alpha, P, ord)
    alphaQ = EvalEndomorphism(alpha, Q, ord)
    a, c = BiDLP(alphaP, P, Q, ord)
    b, d = BiDLP(alphaQ, P, Q, ord)
    #return Matrix(Integers(P.order()), [[a, b],[c,d]])
    return [[a, b],[c,d]]

proof.all(False)

# 512
#p = Integer(167286323857221689112346016933207258999493176647479781908348180838625562682489086996433613891517156716513168813389511283523303877305870043625319556074092350622078682027926127121250390130545025689333161800187261717666012808114506310888652381411304126074312458239)

# 1536
p = Integer(23771452727449498898004557547225354988120021364701398066313318892838230235651289792849929330578203612228877729959916570183072595360843544888617115522542149031138111127294471888062615222219308670025415410876727090874823170475481539518868548152615795234974491096957906735993259268864790983747139027153544490397993988855511337575848418797409297084369809840481360670399419714736596634148006572830334348569025487101216814643806207)

# 2048
#p = Integer(12383345627484975026282352322381621431157228601006950549534308373576108702370134246030714664863826135269122108134969632266058955780067877643619141232561095816669506871207293086566253817648100790943214906574600705093385136118325398200319911685063169182801861057092342096220018419014168812818333860508447888355884884008514670069222088418657620851840717986794126514854813357200934360714804108003400452998106893691991574083281389243586574099783589886418938352208079307208545021269224558144670071437980620093757860952890370611721221573382616667668884022910073760554575618908629958538684405288635125884481680886147749954225954200241067601546030410619884061938102340933311965901014735615296421602529381736952736707884469749887352221690180403199)

F = GF((p,2), name='z2', modulus=var('x')**2 + 1)
sqrtm1 = F.gens()[0]

#filename = 'Precomputed/75-primes.py'
#filename = 'Precomputed/75-primes-512.py'
filename = 'Precomputed/75-primes-1536.py'
#filename = 'Precomputed/75-primes-2048.py'

try:
    open(filename, "r")
    print('precomputed file already exists!')
except:
    assert is_pseudoprime(p)

    print("Choosing torsion....")
    T, facToExt = choose_torsion(p, 1, (2**(0.48*RR(log(p,2)))))

    T = ZZ(T) # Only for the big one
    print(f"Gonna find torsion = {factor(T)}")
    print(facToExt)
    f = (p+1).valuation(2)

    sec_param = 128
    if p == 507227047723007: #Toy paramter
        sec_param = 25
    if p == 136319:
        sec_param = 5
    if p == 8513034219037441780170691209753296498696014329521974009944792576819199999999:
        sec_param = 48 #Just for now to keep things easy
    
    f2 = min(f, 128)
    if f2 < sec_param:
        f3 = Integer(p + 1).valuation(3)
        D_chall = 2**f2*3**f3
    else:
        f3 = 0
        D_chall = 2**f2
    D_com = T.prime_to_m_part(3**f3)

    #assert D_com > p
    assert D_chall >= 2**sec_param

    assert p % 4 == 3
    B = QuaternionAlgebra(-1, -p)
    i, j, k = B.gens()
    O0 = B.quaternion_order([1, i, (i+j)/2, (1+k)/2])

    E0 = EllipticCurve(F, [1,0])
    E0.set_order((p+1)**2, num_checks=0)

    facToBasis = {}
    facToAction = {}

    B_2 = TorsionBasis(E0, 2**f, xOnly = 1)
    P2, Q2 = B_2
    P2.set_order(2**f)
    Q2.set_order(2**f)
    M_0 = ActionMatrix(B(1), B_2, 2**f)
    M_1 = ActionMatrix(i, B_2, 2**f)
    M_2 = ActionMatrix((i+j)/2, B_2, 2**f)
    M_3 = ActionMatrix((1+k)/2, B_2, 2**f)
    M_theta = ActionMatrix(j + (1+k)/2, B_2, 2**f)
    # M_theta = 2*M_2 - M_1 + M_3 #theta = j + (1+k)/2
    BasisAction = (M_0, M_1, M_2, M_3, M_theta)
    facToAction[2**f] = BasisAction

    print("Chall-basis starting...")
    B_chall = TorsionBasis(E0, D_chall, xOnly = 1)
    Pc, Qc = B_chall
    Pc.set_order(D_chall)
    Qc.set_order(D_chall)
    M_0 = ActionMatrix(B(1), B_chall, D_chall)
    M_1 = ActionMatrix(i, B_chall, D_chall)
    M_2 = ActionMatrix((i+j)/2, B_chall, D_chall)
    M_3 = ActionMatrix((1+k)/2, B_chall, D_chall)
    M_theta = ActionMatrix(j + (1+k)/2, B_chall, D_chall)
    # M_theta = 2*M_2 - M_1 + M_3 #theta = j + (1+k)/2
    BasisAction = (M_0, M_1, M_2, M_3, M_theta)
    facToAction[D_chall] = BasisAction

    degToField = {1 : F}
    degToModulus = {1 : str(var('x')**2 + 1)}
    for l, ee in factor(T):
        for e in range(1, ee+1):
            le = l**e
            extdeg = facToExt[le]
            twist = not le.divides(p**extdeg - (-1)**extdeg)
            if extdeg in degToField.keys():
                Fbig = degToField[extdeg]
            else:
                Fbig = F.extension(extdeg,'z'+str(2*extdeg))
                degToField[extdeg] = Fbig
                degToModulus[extdeg] = str(Fbig.modulus())
            Ebig = E0.base_extend(Fbig)
            order = p**extdeg - (-1)**extdeg
            Ebig.set_order(order**2, num_checks=0)
            print(l, e, Fbig)
            if twist:
                print("Registering new field")
                Fbigger = GF((p,4*extdeg), name='temp')
                if not Fbigger.has_coerce_map_from(F):
                    phi1 = Hom(Fbig, Fbigger).an_element()
                    Fbigger.register_coercion(Hom(F, Fbigger)(phi1(Fbig(F.gens()[0]))))
                    if extdeg > 1:
                        Fbigger.register_coercion(phi1)
                assert Fbigger(Fbig(sqrtm1)) == Fbigger(sqrtm1) #THIS ASSERTION IS KEY, IF THIS IS NOT TRUE, THE FIELD EXTENSIONS ARE NOT CONSISTENT
                Ebigger = E0.base_extend(Fbigger)
                P, Q = TorsionBasis(Ebigger, l**e, xOnly = 1)
                PmQ = P-Q
                xP, xQ, xPmQ = xPoint(P.xy()[0], Ebig), xPoint(Q.xy()[0], Ebig), xPoint(PmQ.xy()[0], Ebig)
            else:
                P, Q, xP, xQ, xPmQ = TorsionBasis(Ebig, l**e, xOnly = 2)
            facToBasis[l**e] = [xP, xQ, xPmQ]     
            # Compute basis action
            le = l**e
            Basis = P, Q   
            M_0 = ActionMatrix(B(1), Basis, le)
            M_1 = ActionMatrix(i, Basis, le)
            M_2 = ActionMatrix((i+j)/2, Basis, le)
            M_3 = ActionMatrix((1+k)/2, Basis, le)
            M_theta = ActionMatrix(j + (1+k)/2, Basis, le)
            # M_theta = 2*M_2 - M_1 + M_3 #theta = j + (1+k)/2
            BasisAction = (M_0, M_1, M_2, M_3, M_theta)
            facToAction[l**e] = BasisAction

    printable_facToBasis = {}
    for key in facToBasis.keys():
        printable_facToBasis[key] = [[str(c) for c in P.X] for P in facToBasis[key]]

    with open(filename, 'w') as f:
        f.write(f'{p}\n')
        f.write(f'{(p+1).valuation(2)}\n') #f
        f.write(f'{D_com}\n')
        f.write(f'{D_chall}\n')
        f.write(f'{T}\n')
        f.write(f'{[[str(c) for c in P.xy()] for P in B_2]}\n')
        f.write(f'{[[str(c) for c in P.xy()] for P in B_chall]}\n')
        f.write(f'{facToExt}\n')
        f.write(f'{degToModulus}\n')
        f.write(f'{printable_facToBasis}\n')
        f.write(f'{facToAction}\n')