#code from SIBC;
#implements the algorithm from the SIKE spec


def compute_strategy(ell, n):

    assert(ell == 2)

    p = 5633    #cost for doubling       
    q = 5461   #cost for 2-isogeny evaluation

    S = {1:[]}
    C = {1:0 }
    for i in range(2, n+1):
        b, cost = min(((b, C[i-b] + C[b] + b*p + (i-b)*q) for b in range(1,i)), key=lambda t: t[1])
        S[i] = [b] + S[i-b] + S[b]
        C[i] = cost
    return S[n], C[n]

if __name__=="__main__":
    #print(compute_strategy(2, 518))
    print(compute_strategy(2, 381))