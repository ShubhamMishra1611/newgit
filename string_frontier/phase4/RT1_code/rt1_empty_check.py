"""Fine scan: for a in [-6, -21/29) is there ANY D>3 with all levels n<=N rho-positive?
By monotonicity in D it suffices to test D = 3 + 1e-6 (and D=3+1e-9 for a few)."""
import sys; sys.path.insert(0, '.')
from rt1core import *
from fractions import Fraction as Fr
N = int(sys.argv[1]); step = Fr(1, int(sys.argv[2]))
D0 = Fr(3) + Fr(1, 10**6)
a = Fr(-6); cnt = 0; survivors = []
while a < Fr(-21, 29):
    killed = None
    for n in range(1, N + 1):
        f, k = open_facs(n, a, -a, 0)
        if k == 0:
            continue
        ok, bad = ok_factors(f, k, D0)
        if not ok:
            killed = n; break
    if killed is None:
        survivors.append(a)
    cnt += 1
    a += step
print(f"N={N} step={step}: scanned {cnt} a-values in [-6,-21/29); survivors (all levels positive at D=3+1e-6): {[float(s) for s in survivors]}")
