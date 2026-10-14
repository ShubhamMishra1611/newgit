"""Claim 1 lower endpoint: is the c=0 on-trajectory region really empty for a <= -21/29 (incl. a < -1)?
Exact rational arithmetic, rho rule. For each a: D_max(a) = min_n Dcrit(n;a) over n<=N."""
import sys; sys.path.insert(0, '.')
from rt1core import *
from fractions import Fraction as Fr
N = int(sys.argv[1]) if len(sys.argv) > 1 else 40
avals = [Fr(-k, 50) for k in range(36, 301)]  # -0.72 ... -6
avals += [Fr(-21, 29) - Fr(1, 10**4), Fr(-21, 29) + Fr(1, 10**4)]
for a in sorted(set(avals)):
    best = None
    for n in range(1, N + 1):
        f, k = open_facs(n, a, -a, 0)
        # quick test just above D=3
        ok, bad = ok_factors(f, k, Fr(3) + Fr(1, 10**6)) if k != 0 else (True, [])
        if not ok:
            best = (Fr(3), n, bad)
            break
        d, bad = Dcrit_factors(f, k, tol=Fr(1, 10**6)) if k != 0 else (None, [])
        if d is not None and (best is None or d < best[0]):
            best = (d, n, bad)
    print(f"a={float(a):+.5f}  D_max(N<={N}) = {float(best[0]) if best else 'inf':>10}  set by level {best[1] if best else '-'} spins {best[2] if best else ''}", flush=True)
