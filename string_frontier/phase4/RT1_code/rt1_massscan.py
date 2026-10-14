"""Claim 4: open c=0 family with external mass off trajectory.
 (i)  a=1, m^2=-1+delta : D_max(delta) for levels n<=N
 (ii) exact delta* (level-4 spin-0 = 0 at D=26) via sympy
 (iii) level-1 decoupling edge m^2=-(1+3a)/4: D_max vs (2a+3)^2+1
 (iv) grid over {m^2 <= -a} region (a in [-0.7,1], m^2 in [-(1+3a)/4, -a]) : max D_max
"""
import sys; sys.path.insert(0, '.')
from rt1core import *
from fractions import Fraction as Fr
import sympy as sp
mode = sys.argv[1]; N = int(sys.argv[2])

def Dmax(a, m2, N, tol=Fr(1, 10**8)):
    best = None
    for n in range(1, N + 1):
        f, k = open_facs(n, a, m2, 0)
        if k == 0: continue
        d, bad = Dcrit_factors(f, k, tol=tol)
        if d is not None and (best is None or d < best[0]):
            best = (d, n, bad)
    return best

if mode == 'line':
    for dl in [Fr(0), Fr(1, 100), Fr(1, 20), Fr(1, 10), Fr(1, 5), Fr(1, 4), Fr(3, 10), Fr(7, 20), Fr(2, 5), Fr(4022, 10000),
               Fr(4023, 10000), Fr(9, 20), Fr(1, 2), Fr(3, 5), Fr(3, 4), Fr(1), Fr(5, 4), Fr(3, 2)]:
        b = Dmax(1, -1 + dl, N)
        print(f"a=1 delta={str(dl):>10} m2={float(-1+dl):+.4f}  D_max={float(b[0]):.6f}  level {b[1]} spins {b[2]}", flush=True)
elif mode == 'exact':
    a, D, dl, x = sp.symbols('a D delta x')
    nu = (D - 3) / 2
    def rho0(n, av, m2):
        kap = (n - av - 4 * m2) / 2
        R = sp.expand(sp.prod([av + i + kap * (x - 1) for i in range(1, n + 1)]))
        P = sp.Poly(R, x)
        c0 = 0
        for (k,), rk in P.terms():
            if k % 2 == 0:
                m = k // 2
                c0 += rk * sp.factorial(k) / 2**k * nu / (sp.factorial(m) * sp.rf(nu, k - m + 1))
        return sp.factor(sp.simplify(c0))
    for n in [2, 3, 4, 5]:
        r = rho0(n, 1, -1 + dl)
        print(f"level {n} spin0 at a=1, m2=-1+delta:", r)
        num = sp.numer(sp.together(r))
        sols = [sp.N(s, 20) for s in sp.Poly(num.subs(D, 26), dl).real_roots() if 0 < s < 1.5]
        print("   roots in delta at D=26:", sols)
    # level-1 edge
    m2e = -(1 + 3 * a) / 4
    r2 = rho0(2, a, m2e)
    print("edge level2 spin0:", r2, " -> D bound:", sp.solve(sp.numer(sp.together(r2)), D))
elif mode == 'edge':
    for av in [Fr(-1, 2), Fr(-1, 4), Fr(0), Fr(1, 4), Fr(1, 2), Fr(3, 4), Fr(9, 10), Fr(1), Fr(5, 4), Fr(3, 2), Fr(2)]:
        m2 = -(1 + 3 * av) / 4
        b = Dmax(av, m2, N)
        print(f"edge a={str(av):>5} m2={float(m2):+.4f} D_max={float(b[0]):.6f} level {b[1]} spins {b[2]}  (2a+3)^2+1={float((2*av+3)**2+1):.6f}", flush=True)
elif mode == 'grid':
    best = None
    rows = 0
    for ia in range(-11, 17):           # a = ia/16 in [-0.6875, 1]
        av = Fr(ia, 16)
        lo = -(1 + 3 * av) / 4; hi = -av
        if lo > hi: continue
        K = 32
        for j in range(K + 1):
            m2 = lo + (hi - lo) * Fr(j, K)
            b = Dmax(av, m2, N, tol=Fr(1, 10**6))
            rows += 1
            if b and (best is None or b[0] > best[0]):
                best = (b[0], av, m2, b[1])
            if j in (0, K) or j % 8 == 0:
                print(f"a={float(av):+.4f} m2={float(m2):+.5f} D_max={float(b[0]):.4f} (level {b[1]})", flush=True)
    print(f"grid points={rows}; max D_max = {float(best[0]):.6f} at a={best[1]}, m2={best[2]} (level {best[3]})")
