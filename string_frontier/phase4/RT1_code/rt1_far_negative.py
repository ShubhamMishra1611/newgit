"""Coarse check of emptiness for very negative intercepts a (c=0 and c=1), levels <= N, at D = 3+1e-6."""
import sys; sys.path.insert(0, '.')
from rt1core import *
from fractions import Fraction as Fr
N = int(sys.argv[1])
D0 = Fr(3) + Fr(1, 10**6)
for c in (0, 1):
    for a in [Fr(-7), Fr(-8), Fr(-17, 2), Fr(-10), Fr(-25, 2), Fr(-15), Fr(-20), Fr(-61, 3)]:
        killer = None
        for n in range(1, N + 1):
            if n - c < 1: continue
            f, k = open_facs(n, a, -a, c)
            if k == 0: continue
            ok, bad = ok_factors(f, k, D0)
            if not ok: killer = (n, bad[:4]); break
        print(f"c={c} a={str(a):>6}: first failing level at D=3+1e-6: {killer}", flush=True)
