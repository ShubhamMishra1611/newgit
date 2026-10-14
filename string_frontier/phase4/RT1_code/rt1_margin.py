"""Margin of higher levels: min_{6<=n<=N} Dcrit(n;a) - D_max(a) (c=0 on trajectory), exact bisection tol 1e-7."""
import sys; sys.path.insert(0, '.')
from rt1core import *
from fractions import Fraction as Fr
N = int(sys.argv[1]); avals = [Fr(s) for s in sys.argv[2].split(',')]
for a in avals:
    low = []
    for n in range(1, 6):
        f, k = open_facs(n, a, -a, 0)
        if k == 0: continue
        d, b = Dcrit_factors(f, k, tol=Fr(1, 10**8))
        if d is not None: low.append((d, n))
    Dm = min(low)
    hi = []
    for n in range(6, N + 1):
        f, k = open_facs(n, a, -a, 0)
        if k == 0: continue
        d, b = Dcrit_factors(f, k, tol=Fr(1, 10**7))
        if d is not None: hi.append((d, n, b))
    h = min(hi)
    print(f"a={float(a):+.5f} D_max={float(Dm[0]):.6f} (L{Dm[1]})  min_(6<=n<={N}) Dcrit = {float(h[0]):.6f} at n={h[1]} spins {h[2]}  margin={float(h[0]-Dm[0]):.4f}", flush=True)
