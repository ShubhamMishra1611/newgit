"""T1 (i): margin check -- on a grid, min over levels 6..12 of D_crit(n;a) minus D_max(a) (levels<=5).
A strictly positive margin with small grid step rules out hidden windows where a level >= 6 binds."""
from math import factorial
from flint import fmpq
from pwcore import poch_poly, level_crit, D_of_nu
out = []
def log(*s):
    s = " ".join(str(z) for z in s); print(s, flush=True); out.append(s)
def Dc(a, n):
    nu, j = level_crit(poch_poly(a + 1, n) * fmpq(1, factorial(n)), (n + 3 * a) / 2)
    return float('inf') if nu is None else float(D_of_nu(nu).mid())
worst = (1e9, None)
den = 1000
for k in range(-723, 1001):
    a = fmpq(k, den)
    dmax = min(Dc(a, n) for n in range(1, 6))
    dhi = min(Dc(a, n) for n in range(6, 13))
    if dmax <= 3: continue
    g = dhi - dmax
    if g < worst[0]: worst = (g, float(a), dmax, dhi)
log(f"grid step 1/{den}: min over a of [min_(6<=n<=12) D_crit(n;a) - D_max(a)] = {worst[0]:.5f} at a={worst[1]:.3f} (D_max={worst[2]:.4f}, levels>=6 give {worst[3]:.4f})")
with open(__file__.replace('.py', '_output.txt'), 'w') as f:
    f.write("\n".join(out) + "\n")
