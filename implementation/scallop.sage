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

def customDual(E, phi, dual_gen, deg):
    
    K = phi(dual_gen)
    K.set_order(2**(deg))

    phi_hat = phi.codomain().isogeny(K, algorithm='factored')

    return phi_hat.codomain().isomorphism_to(E) * phi_hat


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

def ActionIdeal(E, K1, K2, Lpos, Lneg):
    # Takes in an effectively oriented curve (E, K1, K2) and an ideal norm L
    # Computes the action of an ideal of norm L

    L = Lpos*Lneg

    print("Computing torsion basis")
    P, Q = torsionBasis(E, L)

    print("Computing action matrix")
    M = actionMatrix(E, K1, K2, P, Q, L)

    E_i = E
    P_i = P
    Q_i = Q

    K1_i = K1
    K2_i = K2

    cofac = L

    for ell, _ in factor(L):
        print(f"Computing isogeny of norm {ell}")
        Mat = Matrix(GF(ell), M)

        #find coeffs and kernel generator
        if Lpos%ell == 1:
            pos = 0
        else:
            pos = 1

        lams = sorted([int(lam) for lam in Mat.eigenvalues()])
        lam = lams[pos] #positive direction = smallest absolute value for instance?

        print(lams)
        
        a, b = (Mat - Matrix(GF(ell), [[lam, 0], [0, lam]])).right_kernel().basis()[0]
        if a != 0:
            ker_gen = P_i + b*(a**(-1))*Q_i
            Q_i *= ell
            P_i *= ell
        else:
            ker_gen = Q_i
            P_i *= ell
            Q_i *= ell

        ker_gen *= (cofac/ell)
        assert ker_gen

        if ell > 300:
            phi_ell = E_i.isogeny(ker_gen, algorithm='velu_sqrt')
        else:
            phi_ell = E_i.isogeny(ker_gen)

        P_i = phi_ell(P_i)
        Q_i = phi_ell(Q_i)
        K1_i = phi_ell(K1_i)
        K2_i = phi_ell(K2_i)

        E_i = phi_ell.codomain()
        cofac /= ell

        P_i.set_order(cofac)
        Q_i.set_order(cofac)

        K1_i.set_order(2**518)
        K2_i.set_order(2**518)

    assert cofac == 1

    return E_i, K1_i, K2_i

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
    
    vec = [-e for e in vec]
    inv = True

    print("\n\n\n~~~~~DONE~~~~~\n")
    print(f"Took a total of {time.time() - t_start} seconds")

    return E_i, K1_i, K2_i

if __name__ == "__main__":
    proof.all(False)

    with open("STARTING_CURVE.txt", "r") as file:
        p = Integer(literal_eval(file.readline()))
        F = GF((p,2), name='z2', modulus=var('x')**2 + 1)
        z2 = F.gens()[0]
        A = F(literal_eval(file.readline()))
        E = EllipticCurve(F, [0, A, 0, 1, 0])
        E.set_order((p+1)**2)
        P = E([F(c) for c in literal_eval(file.readline())])
        Q = E([F(c) for c in literal_eval(file.readline())])

    P.set_order(2**518)
    Q.set_order(2**518)

    es = [randint(-20, 20) for _ in range(75)]

    with open("split_primes.txt", "r") as file:
        ells = file.readline()

    ells = [int(ell) for ell in ells.split(" ")]

    GroupAction(E, P, Q, es, ells)