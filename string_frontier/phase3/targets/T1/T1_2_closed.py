"""
T1 (ii): closed-string partial-wave positivity, exact arithmetic.
 (a) Virasoro-Shapiro, stripped (2408.03362 eq.(1),(3)): R_n(t) = prod_{k=1}^{n-1}(t+k)^2/(n!)^2, s=n, m=0.
 (b) four-dilaton type-II amplitude of 2102.08951 eq.(3.21) (m=1 units):
       M = -8 pi G [(st)^2+(tu)^2+(su)^2] Gamma(-s)Gamma(-t)Gamma(-u)/[Gamma(1+s)Gamma(1+t)Gamma(1+u)],
     residue at s=n: 8 pi G [(nt)^2+(tu)^2+(nu)^2] prod_{k=1}^{n-1}(t+k)^2/(n!)^2,  u=-n-t.
     Footnote 14 there: "experimentally" positive for 2 <= D <= 24.148...
 (c) bosonic closed tachyon amplitude (2201.11575 eq.(3.20)-(3.22)): residue = (open residue)^2,
     R_n = [(t+2)_n/n!]^2 at s=n-1, m^2=-1 (our level n).
Usage: python3 T1_2_closed.py NMAX
"""
import sys, time
from math import factorial
from flint import fmpq, fmpq_poly
from pwcore import poch_poly, level_crit, level_ok_at, D_of_nu, sign_polys, first_neg_root, T

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 60
out = []
def log(*s):
    s = " ".join(str(z) for z in s); print(s, flush=True); out.append(s)

def scan(name, Rf, kf, nmin):
    log(f"\n=== {name}: levels {nmin}..{NMAX} ===")
    t0 = time.time(); rows = []
    for n in range(nmin, NMAX + 1):
        nu, j = level_crit(Rf(n), kf(n))
        D = None if nu is None else D_of_nu(nu)
        rows.append((n, D, j))
        if n <= 10 or n % 10 == 0:
            log(f"  n={n:3d}: D_crit(n) = {'inf' if D is None else D.str(15)}  (spin {j})")
    fin = [(n, D, j) for (n, D, j) in rows if D is not None]
    n0, D0, j0 = min(fin, key=lambda z: float(z[1].mid()))
    log(f"  => D_crit = min_n D_crit(n) = {D0.str(20)} at level {n0}, spin {j0}   ({time.time()-t0:.1f}s)")
    rest = [(n, float(D.mid())) for (n, D, j) in fin if n != n0]
    n1, d1 = min(rest, key=lambda z: z[1])
    log(f"     next-most-constraining level: n={n1}, D_crit={d1:.6f};  D_crit(NMAX)={float(fin[-1][1].mid()):.4f}")
    return n0, j0, D0

# (a) VS stripped
VS = lambda n: (poch_poly(1, n - 1) ** 2) * fmpq(1, factorial(n) ** 2)
nA, jA, DA = scan("(a) Virasoro-Shapiro stripped (type II scalar function)", VS, lambda n: fmpq(n, 2), 1)
# exact form of the critical root
P = sign_polys(VS(nA), fmpq(nA, 2))[jA]
log(f"  critical polynomial in nu=(D-3)/2 (level {nA}, spin {jA}): {P}  -> factor: {P.factor()}")
ok = all(level_ok_at(VS(n), fmpq(n, 2), 24)[0] for n in range(1, NMAX + 1))
log(f"  exact: all levels <= {NMAX} nonnegative at D=24: {ok};  at D=24+1/1000: level {nA} ok? {level_ok_at(VS(nA), fmpq(nA,2), fmpq(24001,1000))[0]}")

# (b) dilaton amplitude (3.21)
def DIL(n):
    s = fmpq(n); t = T; u = -s - t
    pref = (s * t) ** 2 + (t * u) ** 2 + (s * u) ** 2
    return pref * VS(n)
nB, jB, DB = scan("(b) 2102.08951 eq.(3.21) four-dilaton amplitude", DIL, lambda n: fmpq(n, 2), 1)
P = sign_polys(DIL(nB), fmpq(nB, 2))[jB]
log(f"  critical polynomial in nu (level {nB}, spin {jB}): {P}")
log(f"  factorization: {P.factor()}")
log(f"  => critical D = 3 + 2*nu_root = {DB.str(25)}")
for Dt in (fmpq(24), fmpq(24148, 1000), fmpq(24149, 1000), fmpq(25)):
    allok = all(level_ok_at(DIL(n), fmpq(n, 2), Dt)[0] for n in range(1, NMAX + 1))
    log(f"  exact check: all spins of all levels <= {NMAX} nonnegative at D={Dt} ({float(Dt)}): {allok}")

# (c) bosonic closed (KLT square of open bosonic)
BC = lambda n: (poch_poly(2, n) ** 2) * fmpq(1, factorial(n) ** 2)
nC, jC, DC = scan("(c) closed bosonic tachyons (square of open residue)", BC, lambda n: fmpq(n + 3, 2), 1)
P = sign_polys(BC(nC), fmpq(nC + 3, 2))[jC]
log(f"  critical polynomial in nu (level {nC}, spin {jC}): {P.factor()}")
with open(__file__.replace('.py', '_output.txt'), 'w') as f:
    f.write("\n".join(out) + "\n")
