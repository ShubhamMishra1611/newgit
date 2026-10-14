"""
T1 (i): high-level check of the c=0 / c=1 intercept families at selected rational a:
 per-level D_crit(n;a) up to NMAX (rigorous root isolation), and an EXACT check that every
 spin at every level n<=NMAX is nonnegative at the rational point D = D_max(a) (or just below it
 when D_max is irrational).  Also tracks the minimum along fixed daughter trajectories Delta=n-j.
Usage: python3 T1_1_highlevel_check.py c NMAX a1 a2 ...
"""
import sys, time
from math import factorial
from flint import fmpq
from pwcore import poch_poly, level_crit, level_ok_at, D_of_nu, sign_polys, first_neg_root

c = int(sys.argv[1]); NMAX = int(sys.argv[2]); avals = [fmpq(s) for s in sys.argv[3:]]
out = []
def log(*s):
    s = " ".join(str(z) for z in s); print(s, flush=True); out.append(s)
for a in avals:
    t0 = time.time()
    R = lambda n: poch_poly(a + 1, n - c) * fmpq(1, factorial(n))
    kap = lambda n: (n + 3 * a) / 2
    crit = []
    traj = {}   # Delta -> min D_crit over n for spin j = n - c - Delta (top spin is n-c)
    for n in range(1, NMAX + 1):
        if n - c < 0: continue
        P = sign_polys(R(n), kap(n))
        best = None
        for j, Pj in P.items():
            st, r = first_neg_root(Pj)
            Dj = 3.0 if st == 'neg0' else (float('inf') if st == 'ok' else float(D_of_nu(r).mid()))
            Dl = (n - c) - j
            if n >= 6 and Dj < traj.get(Dl, (float('inf'),))[0]:
                traj[Dl] = (Dj, n)
            if best is None or Dj < best[0]:
                best = (Dj, j)
        crit.append((n, best[0], best[1]))
    nset, Dmax, jset = min(crit, key=lambda z: z[1])
    # exact check at a rational D just below/at Dmax
    Dq = fmpq(int(Dmax * 10**6), 10**6)
    if abs(Dmax - round(Dmax * 10**6) / 10**6) < 1e-12:
        Dq = fmpq(round(Dmax * 10**6), 10**6)
    allok = all(level_ok_at(R(n), kap(n), Dq)[0] for n in range(max(1, c), NMAX + 1))
    high = [(n, d) for (n, d, j) in crit if n >= 6]
    log(f"\nc={c}, a={a}: D_max = {Dmax:.9f} set by (level {nset}, spin {jset});  exact check: all spins of levels <= {NMAX} "
        f"nonnegative at D={Dq} ({float(Dq):.6f}): {allok}")
    log(f"   min_(n>=6) D_crit(n) = {min(d for _, d in high):.6f} at n={min(high, key=lambda z: z[1])[0]};"
        f"  D_crit(n) at n=10,20,..: {[(n, round(d, 3)) for (n, d) in high if n % 10 == 0]}")
    log(f"   daughter trajectories Delta=n_top-j (n>=6): min D_crit = " +
        ", ".join(f"D{dl}:{traj[dl][0]:.3f}@n{traj[dl][1]}" for dl in sorted(traj)[:12]))
    log(f"   time {time.time()-t0:.1f}s")
with open(__file__.replace('.py', f'_c{c}_N{NMAX}_output.txt'), 'a') as f:
    f.write("\n".join(out) + "\n")
