"""
T1 task (0): baselines, exact arithmetic (python-flint), levels up to NMAX.
 (a) open bosonic Veneziano Gamma(-1-s)Gamma(-1-t)/Gamma(-2-s-t), tachyon externals m^2=-1.
     Our level n (alpha(s)=n, s=n-1) = level n-1 in the labelling of 2201.11575 / 2608.25183.
     R_n(t) = (t+2)_n / n!, kappa = (s-4m^2)/2 = (n+3)/2.
 (b) type-I stripped function Gamma(-s)Gamma(-t)/Gamma(1-s-t) (massless):
     R_n(t) = (t+1)_{n-1}/n!, s=n, kappa = n/2.   [2201.11575 eq. (1.9), Table 1]
For every level: rigorous critical dimension D_crit(n) (all spins nonneg. for 3<D<=D_crit(n)),
the responsible spin, and an exact check at the critical D (26 resp. 10) listing exact zeros.
"""
import sys, time
from math import factorial
from flint import fmpq
from pwcore import poch_poly, level_crit, level_ok_at, gegenbauer_coeffs, D_of_nu

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 100
out = []
def log(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)

def run(name, Rfun, kfun, nmin, Dcrit_expected):
    log(f"\n=== {name}: levels {nmin}..{NMAX} ===")
    t0 = time.time()
    worst = None
    rows = []
    allok = True
    zeros_all = []
    for n in range(nmin, NMAX + 1):
        R, kap = Rfun(n), kfun(n)
        nu, j = level_crit(R, kap)
        D = None if nu is None else D_of_nu(nu)
        ok, neg, zer = level_ok_at(R, kap, Dcrit_expected)
        okp, negp, _ = level_ok_at(R, kap, fmpq(Dcrit_expected) + fmpq(1, 1000))
        allok &= ok
        if zer: zeros_all.append((n, zer))
        rows.append((n, D, j, ok, negp))
        if D is not None and (worst is None or float(D.mid()) < float(worst[1].mid())):
            worst = (n, D, j)
        if n <= 12 or n % 10 == 0:
            log(f"  n={n:3d}: D_crit(n) = {('inf' if D is None else D.str(12))}  (spin {j});"
                f"  exact all rho_j>=0 at D={Dcrit_expected}: {ok};  negative spins at D={Dcrit_expected}+1/1000: {negp[:6]}")
    log(f"  exact zeros at D={Dcrit_expected} (level, spins): {zeros_all}")
    log(f"  all levels {nmin}..{NMAX} nonnegative at D={Dcrit_expected}: {allok}")
    log(f"  minimum over levels of D_crit(n): level {worst[0]}, spin {worst[2]}, D = {worst[1].str(15)}")
    hi = [(n, float(D.mid())) for (n, D, j, ok, negp) in rows if n >= 4 and D is not None]
    log(f"  min over n>=4 of D_crit(n) = {min(d for _, d in hi):.6f} at n={min(hi, key=lambda z: z[1])[0]}")
    log(f"  D_crit(n) for n = NMAX-4..NMAX: {[(n, round(d,4)) for n, d in hi[-5:]]}")
    log(f"  time {time.time()-t0:.1f}s")
    return rows

rows_b = run("open bosonic (a=1), tachyons", lambda n: poch_poly(2, n) * fmpq(1, factorial(n)),
             lambda n: fmpq(n + 3, 2), 1, 26)
# normalization cross-checks against 2608.25183 eq.(5),(S18) and 2201.11575 Table 1
c = gegenbauer_coeffs(poch_poly(2, 2) * fmpq(1, 2), fmpq(5, 2), 30)
log(f"\n  check 2608.25183 eq.(5): our level-2 spin-0 coefficient at D=30 = {c[0]}, (26-D)/(8(D-1)) = {fmpq(26-30, 8*29)}")
rows_I = run("type-I stripped (massless)", lambda n: poch_poly(1, n - 1) * fmpq(1, factorial(n)),
             lambda n: fmpq(n, 2), 1, 10)
c = gegenbauer_coeffs(poch_poly(1, 4) * fmpq(1, 120), fmpq(5, 2), 7)
D = 7
log(f"  check 2201.11575 Table 1 B_(5,0) at D=7: ours {c[0]}, table {fmpq(9*D*D-250*D+1616, 1920*(D*D-1))}")
with open(__file__.replace('.py', '_output.txt'), 'w') as f:
    f.write("\n".join(out) + "\n")
