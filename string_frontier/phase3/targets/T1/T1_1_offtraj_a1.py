"""T1 (i'): a=1, external mass m^2 = -1 + delta. Fine scan of D_max(delta) on [0, 0.6] (levels <= 30)
and bisection for the first delta* > 0 where D_max(1, delta) exceeds 26."""
from math import factorial
from flint import fmpq
from pwcore import poch_poly, level_crit, D_of_nu
out = []
def log(*s):
    s = " ".join(str(z) for z in s); print(s, flush=True); out.append(s)
N = 30
def Dmax(d):
    best = None
    for n in range(1, N + 1):
        nu, j = level_crit(poch_poly(2, n) * fmpq(1, factorial(n)), (n + 3 - 4 * d) / 2)
        if nu is None: continue
        D = float(D_of_nu(nu).mid())
        if best is None or D < best[0]: best = (D, n, j)
    return best
first = None
for k in range(0, 601):
    d = fmpq(k, 1000)
    b = Dmax(d)
    if k % 25 == 0: log(f"  delta={float(d):.3f}: D_max={b[0]:.6f} (n,j)=({b[1]},{b[2]})")
    if first is None and k > 0 and b[0] > 26 + 1e-12: first = d
log(f"first grid point (step 1e-3) with D_max(1,delta) > 26: delta = {first} ({float(first):.3f}), m^2 = {float(-1+first):.3f}")
lo, hi = first - fmpq(1, 1000), first
for _ in range(30):
    mid = (lo + hi) / 2
    if Dmax(mid)[0] > 26: hi = mid
    else: lo = mid
b = Dmax(hi)
log(f"bisection: delta* = {float(lo):.9f} .. {float(hi):.9f}; just above: D_max={b[0]:.9f} set by (n,j)=({b[1]},{b[2]})")
with open(__file__.replace('.py', '_output.txt'), 'w') as f:
    f.write("\n".join(out) + "\n")
