"""
T1 (i'): c=0 intercept family with externals no heavier than the level-0 resonance:
m^2 = -a + delta with delta <= 0.  Fine grid; D_max(a,delta) from levels <= NMAX.
Tests whether max D over {a <= 1+4 delta, delta <= 0} is 26, attained only at (a,delta)=(1,0).
"""
import sys, time
from math import factorial
from flint import fmpq
from pwcore import poch_poly, level_crit, D_of_nu
NMAX = int(sys.argv[1])
out = []
def log(*s):
    s = " ".join(str(z) for z in s); print(s, flush=True); out.append(s)
def Dmax(a, delta):
    best = None
    for n in range(1, NMAX + 1):
        R = poch_poly(a + 1, n) * fmpq(1, factorial(n)); kap = (n + 3 * a - 4 * delta) / 2
        nu, j = level_crit(R, kap)
        if nu is None: continue
        Dn = float(D_of_nu(nu).mid())
        if best is None or Dn < best[0]: best = (Dn, n, j)
    return best
t0 = time.time()
glob = (0, None)
log(f"levels <= {NMAX}; rows delta, columns a; entries D_max ('--' = no D>3 allowed)")
avals = [fmpq(k, 16) for k in range(-8, 17)]
log(" delta\\a " + "".join(f"{float(a):7.3f}" for a in avals))
for kd in range(0, -17, -1):
    d = fmpq(kd, 32)
    row = []
    for a in avals:
        b = Dmax(a, d)
        if b is None:
            row.append("   inf "); continue
        row.append("    -- " if b[0] <= 3 else f"{b[0]:7.3f}")
        if b[0] > glob[0]: glob = (b[0], (a, d, b[1], b[2]))
    log(f"{float(d):7.4f} " + "".join(row))
log(f"global max over grid: D = {glob[0]:.6f} at a={glob[1][0]}, delta={glob[1][1]} (set by level {glob[1][2]}, spin {glob[1][3]})")
log(f"time {time.time()-t0:.1f}s")
with open(__file__.replace('.py', f'_N{NMAX}_output.txt'), 'w') as f:
    f.write("\n".join(out) + "\n")
