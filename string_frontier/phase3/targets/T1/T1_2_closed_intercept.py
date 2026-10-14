"""
T1 (ii)(d): closed-string 'intercept' families of Virasoro-Shapiro ratio type
   M_{a,S} = Gamma(-al_s)Gamma(-al_t)Gamma(-al_u)/[Gamma(1+al_s)Gamma(1+al_t)Gamma(1+al_u)],  al(x) = a + x,
with S := al_s + al_t + al_u = 3a + 4m^2 an INTEGER (needed for polynomial residues).
Residue at al_s = n (s = n - a):  (-1)^{1+S} [ (al_t + 1)_{n-1-S} ]^2 / (n!)^2   (derived in the report; equals
(al_t+1)_m (al_u+1)_m up to sign, so t<->u symmetric).  External mass m^2 = (S - 3a)/4, kappa = (s-4m^2)/2.
   S = -1: a = 1 is the closed bosonic tachyon amplitude (2201.11575 eq. (3.20));
   S =  0: a = 0 is the stripped Virasoro-Shapiro amplitude.
D_crit(a) = min over levels n <= NMAX of the rigorous per-level critical D.
Usage: python3 T1_2_closed_intercept.py NMAX
"""
import sys, time
from math import factorial
from flint import fmpq
from pwcore import poch_poly, level_crit, D_of_nu
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 30
out = []
def log(*s):
    s = " ".join(str(z) for z in s); print(s, flush=True); out.append(s)
def Dcrit(a, S):
    best = None
    for n in range(1, NMAX + 1):
        m = n - 1 - S
        if m < 0: continue
        R = poch_poly(a + 1, m) ** 2 * fmpq(1, factorial(n) ** 2)
        kap = (n - a - (S - 3 * a)) / 2
        nu, j = level_crit(R, kap)
        if nu is None: continue
        Dn = float(D_of_nu(nu).mid())
        if best is None or Dn < best[0]: best = (Dn, n, j)
    return best
t0 = time.time()
for S in (-1, 0):
    log(f"\n=== S = {S}: externals m^2 = ({S} - 3a)/4; levels <= {NMAX} ===")
    for k in range(-8, 13):
        a = fmpq(k, 8)
        b = Dcrit(a, S)
        log(f"  a={str(a):>5} (m^2={float((S-3*a)/4):+.4f}, level-1 mass^2 = {float(1-a):+.3f}): D_crit = "
            f"{('inf' if b is None else (f'{b[0]:.6f}' if b[0] > 3 else '<=3'))}  set by (n,j)=({b[1] if b else '-'},{b[2] if b else '-'})")
log(f"time {time.time()-t0:.1f}s")
with open(__file__.replace('.py', f'_N{NMAX}_output.txt'), 'w') as f:
    f.write("\n".join(out) + "\n")
