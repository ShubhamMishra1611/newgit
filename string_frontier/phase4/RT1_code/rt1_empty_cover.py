"""For a in [-6,-21/29): record which levels n<=N fail at D=3+1e-6; report, for each level, the a-intervals
(on the grid) where it fails, and check that consecutive grid points always share a common failing level
(so no window can open between grid points unless two failing intervals just touch)."""
import sys; sys.path.insert(0, '.')
from rt1core import *
from fractions import Fraction as Fr
N = int(sys.argv[1]); step = Fr(1, int(sys.argv[2]))
D0 = Fr(3) + Fr(1, 10**6)
a = Fr(-6); rows = []
while a < Fr(-21, 29):
    fails = set()
    for n in range(1, N + 1):
        f, k = open_facs(n, a, -a, 0)
        if k == 0:
            continue
        ok, bad = ok_factors(f, k, D0)
        if not ok:
            fails.add(n)
    rows.append((a, fails))
    a += step
gaps = [(float(rows[i][0]), float(rows[i+1][0])) for i in range(len(rows)-1) if not (rows[i][1] & rows[i+1][1])]
print(f"N={N}, step={step}, points={len(rows)}; min #failing levels per point = {min(len(r[1]) for r in rows)}")
print("adjacent grid pairs with NO common failing level:", gaps)
