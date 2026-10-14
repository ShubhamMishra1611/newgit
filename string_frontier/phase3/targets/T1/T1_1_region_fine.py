"""
T1 (i): the (a, D) unitarity region of the intercept families (exact arithmetic, rigorous root isolation).
  family c=0:  A = G(-al_s)G(-al_t)/G(-al_s-al_t), al(x)=a+x, R_n = (al_t+1)_n/n!,     n>=1
  family c=1:  A = G(-al_s)G(-al_t)/G(1-al_s-al_t),               R_n = (al_t+1)_{n-1}/n!, n>=1 (massive levels)
  externals m^2 = -a + delta (delta=0 by default: externals = level-0 state of the leading trajectory).
For each a on a rational grid, levels n <= NMAX: D_crit(n;a) (all spins >=0 for 3 < D <= D_crit(n;a)),
D_max(a) = min_n D_crit(n;a), the (level, spin) that sets it, and min over levels n >= 5 (to test
whether any high level cuts the region further).  Positivity uses rho_j = c_j/kappa^j (correct sign
below the two-particle threshold).
Usage: python3 T1_1_region.py c NMAX amin amax denom [delta]
"""
import sys, time
from math import factorial
from flint import fmpq
from pwcore import poch_poly, level_crit, D_of_nu

c = int(sys.argv[1]); NMAX = int(sys.argv[2])
amin, amax, den = fmpq(sys.argv[3]), fmpq(sys.argv[4]), int(sys.argv[5])
delta = fmpq(sys.argv[6]) if len(sys.argv) > 6 else fmpq(0)

def R(n, a):
    return poch_poly(a + 1, n - c) * fmpq(1, factorial(n))
def kap(n, a):
    return (n + 3 * a - 4 * delta) / 2

def closed_form(a):
    """D_max from the exact level-1..3 formulas (c=0, delta=0), for comparison."""
    a = float(a); cands = []
    if c == 0:
        if a > 1: return 0.0
        if 0 < a < 2: cands.append((8*a*a + 14*a + 4) / (a*(2 - a)))          # L2 j=0
        if -1 < a < 1: cands.append((26*a + 30) / (3 - a))                     # L3 j=0
        if -3*a*a + 6*a + 1 > 0: cands.append((30*a*a + 48*a + 26) / (-3*a*a + 6*a + 1))  # L3 j=1
    else:
        if a > 0: return 0.0
        if -1 < a < 1: cands.append((8*a + 10) / (1 - a))                      # L3 j=0
    return min(cands) if cands else float('inf')

out = []
def log(*s):
    s = " ".join(str(z) for z in s); print(s, flush=True); out.append(s)
log(f"family c={c}, delta={delta}, levels {1}..{NMAX}, a in [{amin},{amax}] step 1/{den}")
log(f"{'a':>9} | {'D_max(a)':>14} | set by (n,j) | {'L1-3 formula':>12} | {'min_{n>=5} D_crit':>16} at n | D_crit(NMAX)")
a = amin
t0 = time.time()
while a <= amax:
    best = None; best5 = None; last = None
    for n in range(1, NMAX + 1):
        if n - c < 0: continue
        nu, j = level_crit(R(n, a), kap(n, a))
        if nu is None:
            continue
        Dn = float(D_of_nu(nu).mid())
        if n == NMAX: last = Dn
        if best is None or Dn < best[0]:
            best = (Dn, n, j, D_of_nu(nu))
        if n >= 5 and (best5 is None or Dn < best5[0]):
            best5 = (Dn, n)
    cf = closed_form(a)
    Ds = 'inf' if best is None else (f"{best[0]:.9f}" if best[0] > 3 else "<=3 (none)")
    log(f"{str(a):>9} | {Ds:>14} | ({best[1] if best else '-'},{best[2] if best else '-'}) | {cf:12.6f} | "
        f"{(best5[0] if best5 else float('inf')):16.6f} at {best5[1] if best5 else '-'} | {last}")
    a += fmpq(1, den)
log(f"time {time.time()-t0:.1f}s")
tag = f"c{c}_N{NMAX}" + (f"_delta{str(delta).replace('/','over')}" if delta != 0 else "")
with open(__file__.replace('.py', f'_{tag}_output.txt'), 'w') as f:
    f.write("\n".join(out) + "\n")
