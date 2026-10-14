"""Claim 1: D_max(a) for the open c=0 family on trajectory (m^2=-a), rho rule, exact arithmetic.
mode 'grid': for a on a grid, Dcrit(n;a) for n<=NLOW by exact bisection -> D_max, binding (n,spin);
             then verify all levels NLOW<n<=NHIGH are rho-positive at D_max (exact, at the bisection lower bound).
"""
import sys, time; sys.path.insert(0, '.')
from rt1core import *
from fractions import Fraction as Fr
import sympy as sp

A = sp.symbols('a')
def formula(av):
    av = sp.Rational(av)
    a1, a2, a3 = 0.388368678214866756, 0.515763405625102549, 0.522154773496927503
    D = sp.symbols('D')
    if av <= a1:
        return float((26*av+30)/(3-av)), 'L3'
    if av <= a2:
        N5 = (D**2*(av**4 - 4*av**3 - 14*av**2 + 36*av + 45) + D*(90*av**4 + 120*av**3 - 800*av**2 - 2000*av - 1250)
              + 1304*av**4 + 8224*av**3 + 19464*av**2 + 20464*av + 8080)
        return float(min(r for r in sp.Poly(N5, D).nroots() if r.is_real and r > 3)), 'L5'
    if av <= a3:
        N4 = (D**2*(av**4 - 4*av**3 - 4*av**2 + 16*av) + D*(54*av**4 + 36*av**3 - 228*av**2 - 288*av - 64)
              + 296*av**4 + 1336*av**3 + 2368*av**2 + 2000*av + 704)
        return float(min(r for r in sp.Poly(N4, D).nroots() if r.is_real and r > 3)), 'L4'
    return float((8*av**2+14*av+4)/(av*(2-av))), 'L2'

NLOW = int(sys.argv[1]); NHIGH = int(sys.argv[2]); avals = [Fr(s) for s in sys.argv[3].split(',')]
for a in avals:
    t0 = time.time()
    res = []
    for n in range(1, NLOW + 1):
        f, k = open_facs(n, a, -a, 0)
        if k == 0:
            continue
        d, bad = Dcrit_factors(f, k, tol=Fr(1, 10**10))
        if d is not None:
            res.append((d, n, bad))
    res.sort()
    Dmax, nb, sb = res[0]
    second = res[1]
    # verify higher levels at Dmax (lower bound of bisection, so low levels are positive there)
    worst = None
    for n in range(NLOW + 1, NHIGH + 1):
        f, k = open_facs(n, a, -a, 0)
        if k == 0: continue
        ok, bad = ok_factors(f, k, Dmax)
        if not ok:
            worst = (n, bad); break
    fD, fl = formula(a)
    print(f"a={str(a):>8} ({float(a):+.6f}) D_max={float(Dmax):.10f} [level {nb} spins {sb}] formula={fD:.10f} ({fl}) "
          f"| next-binding: level {second[1]} D={float(second[0]):.6f} | levels {NLOW+1}..{NHIGH} at D_max: "
          f"{'ALL POSITIVE' if worst is None else 'FAIL '+str(worst)}  ({time.time()-t0:.1f}s)", flush=True)
