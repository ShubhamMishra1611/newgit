"""
R1 check of the DSSYK chord transfer matrix (2402.00635 eq. (2)):  H|n> = |n+1> + [n]_q |n-1>.
(i) spectrum of the truncated matrix -> support [-2/sqrt(1-q), 2/sqrt(1-q)]  (q-Hermite: E = 2cos(theta)/sqrt(1-q))
(ii) moments <0|H^{2k}|0> vs Touchard-Riordan chord counting (sum over chord diagrams of q^{#crossings})
(iii) internal consistency of 2402.00635: eq.(4) L_A = 2 sin(pi alpha), eq.(5) H = L_A/sqrt(1-q), 2 pi alpha = pi - 2 theta
      => H = 2 cos(theta)/sqrt(1-q), whereas the printed eq.(3) has cos(theta)/sqrt(1-q).
"""
import itertools
import numpy as np
import sympy as sp

def Hmat(q, n):
    H = np.zeros((n, n))
    for k in range(n - 1):
        H[k + 1, k] = 1.0                      # |k> -> |k+1>
        H[k, k + 1] = (1 - q**(k + 1))/(1 - q)  # |k+1> -> [k+1]|k>
    return H

for q in [0.2, 0.5, 0.8]:
    n = 1500  # symmetrised (similar) matrix: off-diagonals sqrt([k]_q); avoids non-normal round-off
    Hs = np.diag(np.sqrt((1 - q**np.arange(1, n))/(1 - q)), 1); Hs = Hs + Hs.T
    ev = np.linalg.eigvalsh(Hs)
    print(f"q={q}: max eig (n=1500) = {ev.max():.5f};  2/sqrt(1-q) = {2/np.sqrt(1-q):.5f};  1/sqrt(1-q) = {1/np.sqrt(1-q):.5f}")

# moments, symbolic in q
qs = sp.symbols('q')
def moment(k):
    # <0|H^k|0> with the non-symmetric matrix (exact)
    n = k + 2
    M = sp.zeros(n, n)
    for j in range(n - 1):
        M[j + 1, j] = 1
        M[j, j + 1] = sum(qs**i for i in range(j + 1))
    v = sp.zeros(n, 1); v[0] = 1
    w = (M**k)*v
    return sp.expand(w[0])

def chord_count(k):
    # sum over perfect matchings of 2k points on a line of q^{#crossings}
    pts = list(range(2*k))
    tot = 0
    def matchings(p):
        if not p:
            yield []
            return
        a = p[0]
        for i in range(1, len(p)):
            b = p[i]
            rest = p[1:i] + p[i+1:]
            for m in matchings(rest):
                yield [(a, b)] + m
    for m in matchings(pts):
        c = sum(1 for (a, b), (c_, d) in itertools.combinations(m, 2) if (a < c_ < b < d) or (c_ < a < d < b))
        tot += qs**c
    return sp.expand(tot)

for k in [1, 2, 3, 4]:
    print(f"<0|H^{2*k}|0> = {moment(2*k)}   chord count = {chord_count(k)}")
th = sp.symbols('theta')
print("eq.(4)+(5) with 2 pi alpha = pi - 2 theta:  L_A/sqrt(1-q) =", sp.simplify(2*sp.sin(sp.pi*(sp.pi - 2*th)/(2*sp.pi))), "/ sqrt(1-q)")
