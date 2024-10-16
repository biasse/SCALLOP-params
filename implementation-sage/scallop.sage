import cypari2
pari = cypari2.Pari()
from ast import literal_eval


def discrete_log_pari(a, base, order):
    r"""
    Wrapper around pari discrete log.
    Works like a.log(b), but allows
    us to use the optional argument
    order.
    """
    x = pari.fflog(a, base, order)
    return ZZ(x)


def BiDLP(R, P, Q, D): 
    r"""
    Given points R, and a D-torsion basis P, Q
    returns a, b s.t. R = [a]P + [b]Q
    """
    ePQ = P.weil_pairing(Q, D, algorithm="pari")
    eRQ = R.weil_pairing(Q, D, algorithm="pari")
    eRP = R.weil_pairing(-P, D, algorithm="pari")

    a = discrete_log_pari(eRQ, ePQ, D)
    b = discrete_log_pari(eRP, ePQ, D)
    assert R == a*P + b*Q
    return a, b


def torsionBasis(E, D):
    F = E.base_field()
    p = F.characteristic()

    cof = (p+1) // D
    facD = factor(D)
    Drad = radical(D) 
    
    while True:
        Psmalls = []
        x = F.random_element()
        try:
            P = E.lift_x(x)*cof
        except:
            continue
        Psmall = P*(D/Drad)
        fullorder = True
        for (l, e) in facD:
            Psmall_i = Psmall*(Drad/l)
            Psmalls.append(Psmall_i)
            if not Psmall_i:
                fullorder = False
                break
        if fullorder:
            P.set_order(D)
            break

    while True:
        x = F.random_element()
        try:
            Q = E.lift_x(x)*cof
        except:
            continue
        Qsmall = Q*(D/Drad)
        basis = True
        for i, (l, e) in enumerate(facD):
            Qsmalli = Qsmall*(Drad/l)
            if not Qsmalli or Psmalls[i].weil_pairing(Qsmalli, l) == 1:
                basis = False
                break
        if basis:
            P.set_order(D)
            Q.set_order(D)
            return P, Q
        
def completeBasis(E, D, P):
    F = E.base_field()
    p = F.characteristic()

    cof = (p+1) // D
    facD = factor(D)
    Drad = radical(D) 

    Psmalls = []
    Psmall = P*(D/Drad)
    for (l, e) in facD:
        Psmall_i = Psmall*(Drad/l)
        Psmalls.append(Psmall_i)

    while True:
        x = F.random_element()
        try:
            Q = E.lift_x(x)*cof
        except:
            continue
        Qsmall = Q*(D/Drad)
        basis = True
        for i, (l, e) in enumerate(facD):
            Qsmalli = Qsmall*(Drad/l)
            if not Qsmalli or Psmalls[i].weil_pairing(Qsmalli, l) == 1:
                basis = False
                break
        if basis:
            Q.set_order(D)
            return Q


def actionMatrix(E, K1, K2, P, Q, order):

    print("     >computing isogenies")
    phi_1 = E.isogeny(K1, algorithm='factored')
    phi_2 = E.isogeny(K2, algorithm='factored')
    isom = phi_1.codomain().isomorphism_to(phi_2.codomain())
    phi_1 = isom * phi_1

    print("     >eval isogenies")

    phi_1P = phi_1(P)
    phi_1Q = phi_1(Q)
    phi_2P = phi_2(P)
    phi_2Q = phi_2(Q)

    deg_inv = pow(phi_2.degree(), -1, order)
    
    print("     >solving DLPs")
    x_1, x_2 = BiDLP(phi_1P, phi_2P, phi_2Q, order)
    x_1 = (deg_inv*x_1) % order
    x_2 = (deg_inv*x_2) % order
    
    x_3, x_4 = BiDLP(phi_1Q, phi_2P, phi_2Q, order)
    x_3 = (deg_inv*x_3) % order
    x_4 = (deg_inv*x_4) % order

    trace = (214524926879081553593184399971232268550466558856614867532685901139338864429220496264877833360002853791876416393265793666992589404517827164678817668824651536) % order #The isomorphism might fuck up the sign, but we can use the trace of omega to repair it
    
    norm = (736335108039604595805923406147184530889923370574768772191969612422073040099331944991573923112581267542507986451953227192970402893063850485730703075899286013451337291468249027691733891486704001513279827771740183629161065194874727962517148100775228363421083691764065477590823919364012917984605619526140822066036736) % order

    print("Trace matrix:")
    print((x_1 + x_4)%order)
    print((-(x_1 + x_4))%order)
    print("Trace real")
    print(trace) #Why are these not the same lol, what am I doing wrong

    print("Norm matrix:")
    print((x_1*x_4 - x_3*x_2)%order)
    print("Norm real")
    print(norm) #Why are these not the same lol, what am I doing wrong

    """
    x_4 = (trace - x_1) % order
    x_3 = ((x_1*x_4 - norm) * pow(x_2, -1, order)) % order
    testPoint = x_3*phi_2Q + x_4*phi_2P
    #if testPoint != phi_2.degree()*phi_1P:
    #    x_1 = -x_1 % order
    #    x_2 = -x_2 % order"""

    return [[x_1, x_3], [x_2, x_4]]


def getKernel(P, Q, L, Lpos, M):   

    ker_gens = []

    lam1s = []
    lam2s = []

    for ell, _ in factor(L):
        Mat = Matrix(GF(ell), M)

        #find coeffs and kernel generator
        if Lpos%ell == 0:
            pos = 0
        else:
            pos = 1

        print(f"ell: {ell}")
        print(f"lams: {Mat.eigenvalues()}")

        lams = sorted([int(lam) for lam in Mat.eigenvalues()])
        lam1s.append(lams[0])
        lam2s.append(lams[1])
        lam = lams[pos] #positive direction = smallest absolute value for instance?
        #print(f"ell: {ell}, eigenvalues: {lams}") <- useful for debugging, these should always be the same...
        
        a, b = (Mat - Matrix(GF(ell), [[lam, 0], [0, lam]])).right_kernel().basis()[0]
        if a != 0:
            ker_gen = P + b*(a**(-1))*Q
        else:
            ker_gen = Q

        ker_gen *= (L/ell)
        assert ker_gen
        ker_gens.append(ker_gen)

    K = sum(ker_gens)
    K.set_order(L)

    print("CRT lams:")
    print(crt(lam1s, [ell for ell, _ in factor(L)]))
    print(crt(lam2s, [ell for ell, _ in factor(L)]))
    
    return K

def getKernelPrecomputedEigen(P, Q, L, Lpos, omega, omega_bar):
    wP = omega(P)
    wQ = omega(Q)

    x_1, x_2 = BiDLP(wP, P, Q, L)    
    x_3, x_4 = BiDLP(wQ, P, Q, L)

    lam1s = []
    lam2s = []
    for ell, _ in factor(L):
        Mat = Matrix(GF(ell), [[x_1, x_3], [x_2, x_4]])

        print(f"ell: {ell}")
        print(f"lams: {Mat.eigenvalues()}")

        lams = sorted([int(lam) for lam in Mat.eigenvalues()])
        lam1s.append(lams[0])
        lam2s.append(lams[1])

    print("CRT lams:")
    lampos = crt(lam1s, [ell for ell, _ in factor(L)])
    lamneg = crt(lam2s, [ell for ell, _ in factor(L)])
    print(lampos)
    print(lamneg)
    print(f"L = {L}")

    Lneg = L//Lpos

    Ppos = omega(P) - (lamneg % L)*P
    Qpos = omega(Q) - (lamneg % L)*Q
    print("these should be true false false true")
    print(omega(Ppos) == lampos*Ppos)
    print(omega(Ppos) == lamneg*Ppos)

    Pneg = omega(P) - (lampos % L)*P
    Qneg = omega(Q) - (lampos % L)*Q
    print(omega(Pneg) == lampos*Pneg)
    print(omega(Pneg) == lamneg*Pneg)

    primary = [Ppos, Pneg]
    backup = [Qpos, Qneg]

    ker_gens = []

    for ell, _ in factor(L):
        dir = 0
        if Lpos%ell != 0:
            dir = 1

        Ki = (L/ell)*primary[dir]
        if not Ki:
            Ki = (L/ell)*backup[dir]
        assert Ki
        ker_gens.append(Ki)

    K = sum(ker_gens)
    K.set_order(L)

    return K



def ActionIdeal(E, P, Q, Lpos, Lneg):
    # Takes in an effectively oriented curve (E, K1, K2) and an ideal norm L
    # Computes the action of an ideal of norm L

    L = Lpos*Lneg

    ell_e = P.order()

    print("Computing torsion basis")
    B1, B2 = torsionBasis(E, L)

    #print("Computing action matrix")
    #M = actionMatrix(E, P, Q, B1, B2, L)

    print("Computing omega")
    phi_P = E.isogeny(P, algorithm='factored')
    phi_Q = E.isogeny(Q, algorithm='factored')

    Qm = completeBasis(E, ell_e, Q)
    phi_Q_dual = phi_Q.codomain().isogeny(phi_Q(Qm), algorithm='factored')
    omega = phi_Q_dual.codomain().isomorphism_to(E) * phi_Q_dual * phi_P.codomain().isomorphism_to(phi_Q.codomain()) * phi_P

    Pm = completeBasis(E, ell_e, P)
    phi_P_dual = phi_P.codomain().isogeny(phi_P(Pm), algorithm='factored')
    omega_bar = phi_P_dual.codomain().isomorphism_to(E) * phi_P_dual * phi_Q.codomain().isomorphism_to(phi_P.codomain()) * phi_Q

    print("Finding kernel generator of ideal")
    #ker_gen = getKernel(B1, B2, L, Lpos, M)
    ker_gen = getKernelPrecomputedEigen(B1, B2, L, Lpos, omega, omega_bar)

    print("Computing isogeny corresponding to ideal")
    phi = E.isogeny(ker_gen, algorithm='factored')

    print("Evaluating isogeny")
    E_l = phi.codomain()

    P_l = phi(P)
    Q_l = phi(Q)

    P_l.set_order(ell_e)
    Q_l.set_order(ell_e)

    return E_l, P_l, Q_l


def GroupAction(E, K1, K2, vec, ells):    
    
    import time
    t_start = time.time()

    E_i = E
    K1_i = K1
    K2_i = K2

    inv = False
    while not all([e <= 0 for e in vec]):
        print("\n\nVector currently looks like:")
        print(vec)
        print(f"Have used {time.time() - t_start} seconds")
        Lpos = 1
        Lneg = 1
        for i in range(len(vec)):
            if vec[i] > 0:
                Lpos *= ells[i]
                vec[i] -= 1
            if vec[i] < 0:
                Lneg *= ells[i]
                vec[i] += 1
        E_i, K1_i, K2_i = ActionIdeal(E_i, K1_i, K2_i, Lpos, Lneg)

    print("\n\n\n~~~~~DONE~~~~~\n")
    print(f"Took a total of {time.time() - t_start} seconds")

    return E_i, K1_i, K2_i


if __name__ == "__main__":
    proof.all(False)

    #param_lvl = "1024"
    #param_lvl = "512"
    #param_lvl = "1536"
    param_lvl = "1536-140"

    with open("params_" + param_lvl + ".txt", "r") as file:
        p = Integer(literal_eval(file.readline()))
        F = GF((p,2), name='z2', modulus=var('x')**2 + 1)
        z2 = F.gens()[0]
        A = F(literal_eval(file.readline()))
        E = EllipticCurve(F, [0, A, 0, 1, 0])
        E.set_order((p+1)**2)
        P = E([F(c) for c in literal_eval(file.readline())])
        Q = E([F(c) for c in literal_eval(file.readline())])
        ells = literal_eval(file.readline())
        e = literal_eval(file.readline())

    P.set_order(2**e)
    Q.set_order(2**e)
    Qm = completeBasis(E, 2**e, Q)

    print(E)
    print((P).xy()[0])
    print((Q).xy()[0])
    print((Qm).xy()[0])

    #es = [randint(-5, 5) for _ in range(75)]
    es = [3 for _ in range(75)]
    #es = [0,-3,5,-2,-2,-1,1,16,14,6,5,-17,16,27,8,-34,6,9,1,2,19,-24,21,35,-2,41,-11,-5,60,-11,80,6,20,13,15,8,22,2,-21,-12,7,-19,-68,-39,9,-68,-13,33,2,-3,-6,-139,-1,-4,26,10,-6,1,18,-13,-31,13,14,-6,32,-14,6,0,-3,8,2,6,-4,0,11]
    GroupAction(E, P, Q, es, ells)

    """
    alice = [randint(-3, 3) for _ in range(75)]
    bob = [randint(-3, 3) for _ in range(75)]

    E_A, P_A, Q_A = GroupAction(E, P, Q, alice.copy(), ells)
    print("\n\n\n ALICE PUBLIC KEY DONE")
    E_B, P_B, Q_B = GroupAction(E, P, Q, bob.copy(), ells)
    print("\n\n\n BOB PUBLIC KEY DONE")

    E_BA, _, _ = GroupAction(E_B, P_B, Q_B, alice.copy(), ells)
    print("\n\n\n ALICE SHARED SECRET DONE")
    E_AB, _, _ = GroupAction(E_A, P_A, Q_A, bob.copy(), ells)
    print("\n\n\n BOB SHARED SECRET DONE")

    print(E_AB.j_invariant())
    print(E_BA.j_invariant())

    assert E_AB.j_invariant() == E_BA.j_invariant()
    """
