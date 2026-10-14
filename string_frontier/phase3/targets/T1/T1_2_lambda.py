"""
T1 (ii)(c): the Cheung-Hillman-Remmen lambda-families (2408.03362).
  closed M (eq.(7)), residues (eq.(13)):  R(n,t) = [ ((1+lam)/2 + lam t)_{n-1} / ((1+3lam)/2)_{n-1} ]^2,
  planar A (eq.(18)), residues (eq.(16)): Rhat(n,t) = ((1+lam)/2 + lam t)_{n-1} / ((1+3lam)/2)_{n-1},
  poles at s = mu(n) = (n+lam-1)/lam (eq.(5)), massless externals, n >= 1.
  lam = 1: M = VS (stripped), A = type-I stripped (eq.(20)).
For each rational lam: D_crit(lam) = min over levels n <= NMAX of the rigorous per-level critical D,
and the (level, spin) that sets it.
Usage: python3 T1_2_lambda.py NMAX
"""
import sys, time
from flint import fmpq, fmpq_poly
from pwcore import level_crit, D_of_nu, level_ok_at

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 30
out = []
def log(*s):
    s = " ".join(str(z) for z in s); print(s, flush=True); out.append(s)

def Rhat(n, lam):
    # prod_{k=0}^{n-2} (t + ((1+lam)/2 + k)/lam)  (positive constant dropped)
    p = fmpq_poly([1])
    for k in range(n - 1):
        p = p * fmpq_poly([((1 + lam) / 2 + k) / lam, 1])
    return p

def Dcrit(lam, closed):
    best = None
    for n in range(1, NMAX + 1):
        R = Rhat(n, lam) ** 2 if closed else Rhat(n, lam)
        kap = (n + lam - 1) / lam / 2
        nu, j = level_crit(R, kap)
        if nu is None: continue
        Dn = float(D_of_nu(nu).mid())
        if best is None or Dn < best[0]: best = (Dn, n, j)
    return best

lams = [fmpq(1, 16), fmpq(1, 8), fmpq(1, 4), fmpq(3, 8), fmpq(1, 2), fmpq(5, 8), fmpq(3, 4), fmpq(7, 8), fmpq(15, 16),
        fmpq(1), fmpq(17, 16), fmpq(9, 8), fmpq(5, 4), fmpq(3, 2), fmpq(7, 4), fmpq(2), fmpq(3), fmpq(4), fmpq(8), fmpq(16)]
t0 = time.time()
for closed in (True, False):
    log(f"\n=== {'closed M (eq.7)' if closed else 'planar A (eq.18)'}: D_crit(lambda), levels <= {NMAX} ===")
    for lam in lams:
        b = Dcrit(lam, closed)
        log(f"  lambda={str(lam):>6} ({float(lam):7.4f}): D_crit = {('inf' if b is None else f'{b[0]:.6f}'):>12}  set by (n,j)=({b[1] if b else '-'},{b[2] if b else '-'})")
# CHR statements: M positive in D=4 for all lam>=0 ; A positive in D=4 for all lam>=0
log("\nexact checks of 2408.03362 statements at D=4 (levels <= NMAX):")
for lam in lams:
    okM = all(level_ok_at(Rhat(n, lam) ** 2, (n + lam - 1) / lam / 2, 4)[0] for n in range(1, NMAX + 1))
    okA = all(level_ok_at(Rhat(n, lam), (n + lam - 1) / lam / 2, 4)[0] for n in range(1, NMAX + 1))
    log(f"  lambda={str(lam):>6}: M positive at D=4: {okM};  A positive at D=4: {okA}")
log(f"time {time.time()-t0:.1f}s")
with open(__file__.replace('.py', f'_N{NMAX}_output.txt'), 'w') as f:
    f.write("\n".join(out) + "\n")
