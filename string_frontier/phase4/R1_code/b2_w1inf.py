"""R1: 2105.14346 eqs (2.5)-(2.7). Check that (2.5) with the rescaling (2.6) gives (2.7) exactly
(exact rational arithmetic, kappa symbolic), and Jacobi for the wedge algebra (2.7),(2.9)."""
from fractions import Fraction as F
from math import factorial
import itertools, sympy as sp
kap = sp.symbols('kappa')
def fac(x):
    return sp.factorial(x) if x >= 0 else sp.zoo
half = [F(k, 2) for k in range(2, 9)]          # p = 1, 3/2, ..., 4
def modes(p):
    return [p - 1 - j for j in range(int(2*p - 1))]   # 1-p <= m <= p-1, step 1
bad = 0; n_checked = 0
for p, q in itertools.product(half, half):
    k, l = 4 - 2*p, 4 - 2*q
    for m in modes(p):
        for n in modes(q):
            P = p + q - 2
            if not (1 - P <= m + n <= P - 1):    # target outside wedge -> (2.7) coefficient must vanish
                c27 = m*(q - 1) - n*(p - 1)
                if c27 != 0: bad += 1
                continue
            a, b = (2 - k)/2, (2 - l)/2
            c25 = -kap/2*(n*(2 - k) - m*(2 - l))*fac(int(a - m + b - n - 1))/(fac(int(a - m))*fac(int(b - n))) \
                  * fac(int(a + m + b + n - 1))/(fac(int(a + m))*fac(int(b + n)))
            # w^p_m = (1/kap)(p-m-1)!(p+m-1)! H^{4-2p}_m ;  H^{4-2P}_{m+n} = kap w^P_{m+n}/((P-m-n-1)!(P+m+n-1)!)
            lhs = (1/kap)**2*fac(int(p - m - 1))*fac(int(p + m - 1))*fac(int(q - n - 1))*fac(int(q + n - 1))*c25 \
                  * kap/(fac(int(P - m - n - 1))*fac(int(P + m + n - 1)))
            n_checked += 1
            if sp.simplify(lhs - (m*(q - 1) - n*(p - 1))) != 0: bad += 1
print(f"(2.5)+(2.6) -> (2.7): checked {n_checked} structure constants, mismatches: {bad}")
def br(x, y):
    (p, m), (q, n) = x, y
    c = m*(q - 1) - n*(p - 1)
    return {} if c == 0 else {(p + q - 2, m + n): c}
def brlin(X, y):
    out = {}
    for g, c in X.items():
        for g2, c2 in br(g, y).items():
            out[g2] = out.get(g2, 0) + c*c2
    return out
gens = [(p, m) for p in half[:5] for m in modes(p)]
fails = 0; cnt = 0
for x, y, z in itertools.combinations(gens, 3):
    tot = {}
    for A, B, C in [(x, y, z), (y, z, x), (z, x, y)]:
        for g, c in brlin(br(A, B), C).items():   # [[A,B],C]
            tot[g] = tot.get(g, 0) + c
    cnt += 1
    if any(v != 0 for v in tot.values()): fails += 1
print(f"Jacobi on {cnt} triples of wedge generators (p<=3): failures {fails}")
