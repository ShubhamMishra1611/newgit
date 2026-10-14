"""RT3 check of T3 claims 1-2 (completeness + face identity).

(i)  Identity  #{n in Z_{>=0}^J : X(n) <= T} = 2^-J sum_{U subset [J]} #(Z^U cap E_T)
     checked by brute force for J = 2..5 and random T.
(ii) Brute force vs my MITM (rt3_bp_core.Shell) for small J.
(iii) Completeness of T3's J=4 (seed 11, R=1500) window counts: an INDEPENDENT
     algorithm (loop over n1, vectorise (n2,n3), solve n4 analytically) lists every
     orthant point with |Lambda| <= eps_max and the counts N[0,eps], N[-eps,0] are
     compared with T3's printed table (T3_A1_bp_discretuum.out).
"""
import sys, os, math, itertools, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rt3_bp_core import Shell, smooth_count

rng = np.random.default_rng(2026)

print("(i) identity check")
for J in (2, 3, 4, 5):
    for trial in range(3):
        c = (1 + rng.random(J)) ** 2
        T = float(rng.uniform(20, 60))
        nm = int(math.isqrt(int(T))) + 2
        # LHS
        lhs = 0
        full = {}
        for n in itertools.product(range(-nm, nm + 1), repeat=J):
            X = sum(ci * k * k for ci, k in zip(c, n))
            if X <= T:
                if all(k >= 0 for k in n):
                    lhs += 1
        # RHS: for each subset U count full-lattice points supported in U (zero outside)
        rhs = 0
        for r in range(J + 1):
            for U in itertools.combinations(range(J), r):
                cnt = 0
                for m in itertools.product(range(-nm, nm + 1), repeat=r):
                    X = sum(c[i] * k * k for i, k in zip(U, m))
                    if X <= T:
                        cnt += 1
                rhs += cnt
        print(f"  J={J} T={T:6.2f}: orthant count {lhs}, 2^-J * sum_U #(Z^U cap E) = {rhs}/{2**J} = {rhs/2**J}"
              f"  -> {'OK' if rhs == lhs * 2**J else 'MISMATCH'};  smooth face formula {smooth_count(c, T):.1f}")

print("\n(ii) MITM vs brute force (J=3..5) on random windows")
for J in (3, 4, 5):
    c = (1 + rng.random(J)) ** 2
    S = 150.123
    sh = Shell(c, S + 30)
    nm = int(math.isqrt(int((S + 30) / c.min()))) + 1
    allX = []
    for n in itertools.product(range(0, nm + 1), repeat=J):
        X = sum(ci * k * k for ci, k in zip(c, n))
        if X <= S + 30:
            allX.append(X)
    allX = np.sort(np.array(allX))
    for lo, hi in ((S - 5, S + 5), (S, S + 0.4), (S - 30, S + 30)):
        bf = int(np.sum((allX >= lo) & (allX < hi)))
        print(f"  J={J} [{lo:.2f},{hi:.2f}): brute {bf}, mitm {sh.count(lo, hi)}")

print("\n(iii) completeness of T3's J=4 table (seed 11, R=1500), independent 3-loop+analytic algorithm")
t0 = time.time()
q = 1.0 + np.random.default_rng(11).random(4)   # T3's charge recipe (re-implemented)
c = q * q
S = 1500.0 ** 2 + 0.123456789
# T3's eps grid: geomspace(3/rho_sm, 0.02*S, 22); recompute rho_sm with my own formula
from rt3_bp_core import smooth_rho
rho = smooth_rho(c, S)
eps_grid = np.geomspace(3.0 / rho, 0.02 * S, 22)
epsmax = eps_grid[8]          # first 9 windows (eps <= 0.0226) are listed point by point here
lo_all, hi_all = S - 2 * epsmax, S + 2 * epsmax
# (n2,n3) table
n2m = int(math.isqrt(int(hi_all / c[1]))) + 1
n3m = int(math.isqrt(int(hi_all / c[2]))) + 1
N2, N3 = np.meshgrid(np.arange(n2m + 1), np.arange(n3m + 1), indexing="ij")
S23 = (c[1] * N2 * N2 + c[2] * N3 * N3).ravel()
S23 = np.sort(S23[S23 <= hi_all])
vals = []
n1m = int(math.isqrt(int(hi_all / c[0]))) + 1
for n1 in range(n1m + 1):
    base = c[0] * n1 * n1
    if base > hi_all:
        break
    s = base + S23[: np.searchsorted(S23, hi_all - base, side="right")]
    t = (hi_all - s) / c[3]
    n4 = np.floor(np.sqrt(np.maximum(t, 0))).astype(np.int64)
    x = s + c[3] * n4 * n4
    # guard against float floor error: also test n4+1 and n4-1
    for d in (-1, 0, 1):
        nn = n4 + d
        ok = nn >= 0
        xx = s[ok] + c[3] * nn[ok] * nn[ok]
        sel = (xx >= lo_all) & (xx < hi_all)
        vals.append(xx[sel])
vals = np.sort(np.concatenate(vals))  # the -1,0,+1 guard gives distinct n4 per s, so no duplicates (np.unique would merge float64 ties)
Lam = (vals - S) / 2
T3_pos = [2, 10, 24, 80, 257, 838, 2462, 7320, 22550]
T3_neg = [5, 15, 39, 87, 247, 732, 2345, 7267, 22579]
print(f"  listed {len(vals)} orthant points with |Lambda| < {epsmax:.4g}  [{time.time()-t0:.1f}s]")
print("   eps          N[0,eps] RT3   T3   |  N[-eps,0] RT3   T3")
allok = True
for k in range(9):
    e = eps_grid[k]
    Np = int(np.sum((Lam >= 0) & (Lam < e)))
    Nm = int(np.sum((Lam < 0) & (Lam >= -e)))
    ok = (Np == T3_pos[k]) and (Nm == T3_neg[k])
    allok &= ok
    print(f"   {e:10.3e}  {Np:9d}  {T3_pos[k]:9d}   |  {Nm:9d}  {T3_neg[k]:9d}   {'OK' if ok else 'DIFF'}")
print("  all identical:" , allok)
# local density in that window vs smooth (faces) vs BP
print(f"  density over |Lambda|<{epsmax:.3g}: {len(Lam)/(2*epsmax):.5e}; smooth {rho:.5e}")
