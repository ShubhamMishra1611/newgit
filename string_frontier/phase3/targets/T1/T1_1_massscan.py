"""
T1 (i'): two-parameter family (a, m^2): c=0 intercept amplitude G(-al_s)G(-al_t)/G(-al_s-al_t),
al(x)=a+x, with external scalars of arbitrary mass m^2 = -a + delta (delta=0: on the trajectory).
Computes D_max(a, delta) = min over levels n<=NMAX of the rigorous per-level critical D.
Modes:
  line_a1   : a = 1, delta in a grid                       (vary only the external mass)
  edge      : delta = (a-1)/4 (level-1 scalar decouples: rho_0(level 1) = 0), a in a grid
  grid      : 2D grid
Usage: python3 T1_1_massscan.py MODE NMAX
"""
import sys, time
from math import factorial
from flint import fmpq
from pwcore import poch_poly, level_crit, D_of_nu

mode, NMAX = sys.argv[1], int(sys.argv[2])
out = []
def log(*s):
    s = " ".join(str(z) for z in s); print(s, flush=True); out.append(s)

def Dmax(a, delta):
    best = None
    for n in range(1, NMAX + 1):
        R = poch_poly(a + 1, n) * fmpq(1, factorial(n))
        kap = (n + 3 * a - 4 * delta) / 2
        nu, j = level_crit(R, kap)
        if nu is None: continue
        Dn = float(D_of_nu(nu).mid())
        if best is None or Dn < best[0]:
            best = (Dn, n, j)
    return best

t0 = time.time()
if mode == 'line_a1':
    log(f"a=1, m^2 = -1 + delta; levels <= {NMAX}")
    for k in list(range(0, 41)) + list(range(44, 121, 4)):
        d = fmpq(k, 80)
        b = Dmax(fmpq(1), d)
        log(f"  delta={str(d):>6} (m^2={float(-1+d):+.4f}): D_max = {b[0]:.6f}  set by (n,j)=({b[1]},{b[2]})" if b else f"  delta={d}: no bound")
elif mode == 'edge':
    log(f"edge of the level-1 constraint: delta=(a-1)/4, i.e. m^2 = -(1+3a)/4; levels <= {NMAX}")
    for k in range(-4, 25):
        a = fmpq(k, 8)
        d = (a - 1) / 4
        b = Dmax(a, d)
        log(f"  a={str(a):>5}: m^2={float(-a+d):+.4f}: D_max = {b[0]:.6f} set by (n,j)=({b[1]},{b[2]});  level-2 formula 2(2a^2+6a+5) = {float(2*(2*a*a+6*a+5)):.4f}" if b else f"  a={a}: no bound")
elif mode == 'grid':
    log(f"2D grid: rows a, columns m^2; entry = D_max (levels <= {NMAX}); '  --  ' = excluded for all D>3")
    m2s = [fmpq(k, 4) for k in range(-6, 3)]
    log("   a \\ m^2 " + "".join(f"{float(m):>9.2f}" for m in m2s))
    for k in range(-2, 9):
        a = fmpq(k, 4)
        row = []
        for m2 in m2s:
            d = m2 + a
            b = Dmax(a, d)
            row.append("    inf  " if b is None else ("    --   " if b[0] <= 3 else f"{b[0]:9.3f}"))
        log(f"  {float(a):6.2f}   " + "".join(row))
log(f"time {time.time()-t0:.1f}s")
with open(__file__.replace('.py', f'_{mode}_N{NMAX}_output.txt'), 'w') as f:
    f.write("\n".join(out) + "\n")
