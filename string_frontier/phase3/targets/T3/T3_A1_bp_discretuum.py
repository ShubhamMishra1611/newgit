"""T3 (A1): Bousso-Polchinski discretuum near Lambda = 0 -- density, level statistics,
minimum |Lambda| scaling, sensitivity to the charges.

Run:  python3 -I T3_A1_bp_discretuum.py > T3_A1_bp_discretuum.out   (~10-20 min, <4 GB RAM)

Model (BP hep-th/0004134 eq. (2.21)):  Lambda = lambda_bare + (1/2) sum n_i^2 q_i^2, lambda_bare<0.
Notation: R^2 = 2|lambda_bare|, Q = sum q_i^2 n_i^2, Lambda = (Q - R^2)/2.  Units: arbitrary
(BP's (sqrt2 kappa4)^4 = 1 not needed; everything below is a statement about the quadratic form).
We count DISTINCT Lambda values, i.e. orthant points n_i >= 0 (BP's sign degeneracy D = 2^J
removed).  "Flux cutoff" = the largest flux reached on the shell, n_max = R/q_min.

Predictions compared against:
  * BP (2.26) with D = 2^J:  rho_BP = omega_{J-1} R^{J-2} / (2^J prod q_i)   [per unit Lambda]
  * our smooth orthant formula (bp_lib.smooth_orthant_*), which adds the coordinate-face
    corrections:  N_orth(X) ~ 2^{-J} sum_U V_|U| X^{|U|/2}/prod_U q.
  * Poisson (Berry-Tabor) level statistics for a generic flat-torus spectrum.
"""
import sys, os, math, time, json
import numpy as np
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bp_lib import (BPModel, smooth_orthant_count, smooth_orthant_density_Lambda,
                    bp_leading_density)

T0 = time.time()
OUT = {}


def charges(J, seed, kind="random", delta=None):
    rng = np.random.default_rng(seed)
    u = rng.random(J)
    if kind == "random":
        q = 1.0 + u            # q_i in [1,2): incommensurate O(1) charges
        return q * q
    if kind == "nearcomm":   # q_i^2 = 1 + delta*u_i
        return 1.0 + delta * u
    raise ValueError


def hdr(s):
    print("\n" + "=" * 90 + "\n" + s + "\n" + "=" * 90, flush=True)


# --------------------------------------------------------------------------------------
hdr("PART (a)+(b): counting function N(eps) and level statistics near Lambda = 0")
CONFIGS = [  # (J, R, seed)
    (4, 1500.0, 11),
    (6, 330.0, 12),
    (8, 90.0, 13),
    (10, 46.0, 14),
    (12, 30.0, 15),
]
LEVELS_TARGET = 600_000
summary_ab = []
for J, R, seed in CONFIGS:
    t = time.time()
    q2 = charges(J, seed)
    R2 = float(R * R + 0.123456789)      # generic shell radius
    rho_bp = bp_leading_density(q2, R2)
    rho_sm = smooth_orthant_density_Lambda(q2, R2, 0.0)
    eps_max = 0.02 * R2                  # largest window used for N(eps)
    m = BPModel(q2, R2, 2 * eps_max + 1e-6)
    print(f"\nJ={J}  R={R}  q_i^2 in [{q2.min():.3f},{q2.max():.3f}]  n_max=R/q_min={R/np.sqrt(q2.min()):.0f}"
          f"  half sizes {len(m.A):,}/{len(m.B):,}")
    print(f"  density at Lambda=0: BP(2.26),D=2^J: {rho_bp:.4e};  smooth orthant (with faces): {rho_sm:.4e}"
          f"  -> mean spacing {1/rho_sm:.3e}")
    # (a) N(eps) for Lambda in [0, eps] and [-eps, 0]
    eps_grid = np.geomspace(3.0 / rho_sm, eps_max, 22)
    print("   eps          N[0,eps]    N_pred     ratio   (N-Np)/sqrt(Np) |  N[-eps,0]  ratio  | local slope dlnN/dlneps")
    Np_prev = None
    rows = []
    for k, eps in enumerate(eps_grid):
        N = m.count_Q_window(R2, R2 + 2 * eps)
        Nm = m.count_Q_window(R2 - 2 * eps, R2)
        Np = float(smooth_orthant_count(q2, R2 + 2 * eps) - smooth_orthant_count(q2, R2))
        Npm = float(smooth_orthant_count(q2, R2) - smooth_orthant_count(q2, R2 - 2 * eps))
        rows.append((eps, N, Np, Nm, Npm))
    rows = np.array(rows)
    for k, (eps, N, Np, Nm, Npm) in enumerate(rows):
        if k > 0 and rows[k - 1, 1] > 0 and N > 0:
            slope = math.log(N / rows[k - 1, 1]) / math.log(eps / rows[k - 1, 0])
        else:
            slope = float("nan")
        print(f"   {eps:10.3e}  {int(N):11d}  {Np:11.1f}  {N/Np:7.4f}  {(N-Np)/math.sqrt(Np):+7.2f}      |  {int(Nm):10d} {Nm/Npm:7.4f} |  {slope:6.3f}")
    # power-law fit on the range with N >= 10 and eps <= 1e-3 R^2
    sel = (rows[:, 1] >= 10) & (rows[:, 0] <= 1e-3 * R2)
    p, c = np.polyfit(np.log(rows[sel, 0]), np.log(rows[sel, 1]), 1)
    print(f"  fitted N ~ eps^p on [N>=10, eps<=1e-3 R^2]: p = {p:.5f}  ({sel.sum()} points)")

    # (b) level statistics: list ~LEVELS_TARGET levels centred on Lambda = 0
    W = 0.5 * LEVELS_TARGET / rho_sm
    Q, _, _ = m.list_Q_window(R2 - 2 * W, R2 + 2 * W, max_points=5 * LEVELS_TARGET)
    Lam = (Q - np.longdouble(R2)) / 2
    # unfold with the smooth density expanded to 2nd order about Lambda=0 (window << R^2, so
    # this is exact to O(W^2/R^4)); rescaled so the empirical mean spacing is 1.
    h = 1e-6 * R2
    drho = (smooth_orthant_density_Lambda(q2, R2, h) - smooth_orthant_density_Lambda(q2, R2, -h)) / (2 * h)
    Lf = Lam.astype(np.float64)
    x = rho_sm * Lf + 0.5 * drho * Lf * Lf
    x = x - x[0]
    x = x * (len(x) - 1) / x[-1]
    s = np.diff(x)
    s = s / s.mean()
    M = len(s)
    ks = stats.kstest(s, "expon")
    # Wigner surmise (GOE) for reference
    frac01 = np.mean(s < 0.1); frac001 = np.mean(s < 0.01)
    print(f"  level statistics on {M:,} spacings (Lambda in [{float(Lam[0]):.3e},{float(Lam[-1]):.3e}]):")
    print(f"     var(s) = {s.var():.4f} (Poisson 1, GOE 0.286);  <s^2> = {np.mean(s*s):.4f} (Poisson 2, GOE 4/pi=1.273)")
    print(f"     P(s<0.1) = {frac01:.4f} (Poisson {1-math.exp(-0.1):.4f}, GOE {1-math.exp(-math.pi*0.01/4):.4f});"
          f"  P(s<0.01) = {frac001:.5f} (Poisson {1-math.exp(-0.01):.5f})")
    print(f"     KS vs Exp(1): D = {ks.statistic:.5f}, p = {ks.pvalue:.3g};  max s = {s.max():.2f}"
          f" (Poisson expectation ~ ln M + 0.577 = {math.log(M)+0.5772:.2f})")
    # empirical density vs prediction over the listed window
    emp_rho = (len(x) - 1) / (float(Lam[-1]) - float(Lam[0]))
    print(f"     empirical mean density over window / smooth prediction = {emp_rho/rho_sm:.5f};"
          f"  / BP(2.26) = {emp_rho/rho_bp:.4f}")
    # number variance Sigma^2(L)
    Ls = [0.5, 1, 2, 5, 10, 30, 100, 300, 1000, 3000, 10000, 30000]
    rng = np.random.default_rng(7)
    sig2 = []
    for Lw in Ls:
        if Lw > 0.2 * x[-1]:
            break
        st = rng.uniform(0, x[-1] - Lw, 20000)
        cnt = np.searchsorted(x, st + Lw) - np.searchsorted(x, st)
        sig2.append((Lw, cnt.var()))
    print("     number variance Sigma^2(L) [Poisson: L]:  " + "  ".join(f"L={a:g}:{b:.3g}" for a, b in sig2))
    # histogram of s for the record
    hist, edges = np.histogram(s, bins=np.linspace(0, 6, 61), density=True)
    OUT[f"J{J}"] = dict(R=R, R2=R2, rho_bp=rho_bp, rho_smooth=rho_sm, emp_rho=emp_rho,
                        fit_p=p, var_s=float(s.var()), P_s_lt_0p1=float(frac01), KS_D=float(ks.statistic),
                        KS_p=float(ks.pvalue), sigma2=sig2, hist=hist.tolist(), edges=edges.tolist(),
                        Neps=rows.tolist())
    summary_ab.append((J, R, M, p, s.var(), frac01, ks.statistic, ks.pvalue, emp_rho / rho_sm, sig2[-1]))
    print(f"  [{time.time()-t:.1f}s]", flush=True)
    del m, Q, Lam, x, s

hdr("PART (c): minimum |Lambda| vs J and flux cutoff (exactly re-verified with Fractions)")
# For each (J,R): NSHIFT generic shell radii R^2 + u, u ~ U[0,1).  Record exact min|Lambda|,
# and the density estimated from the 100 smallest |Lambda| (order statistics of a Poisson process).
NSHIFT = 6
MINCONF = {
    4: [200.0, 600.0, 1800.0, 5400.0],
    6: [40.0, 80.0, 160.0, 320.0],
    8: [20.0, 35.0, 60.0, 100.0],
    10: [14.0, 22.0, 34.0, 48.0],
    12: [10.0, 15.0, 22.0, 30.0],
    16: [8.0, 11.0, 15.0, 19.0],
    20: [7.0, 9.0, 11.5, 13.5],
    24: [6.0, 8.0, 10.0, 11.5],
}
rows_c = []
print(" J     R    n_max  |  <min|Lam|>(exact)   1/(2 rho_sm)   ratio | rho_emp(k<=100)/rho_sm | min over shifts")
for J, Rlist in MINCONF.items():
    q2 = charges(J, 100 + J)
    for R in Rlist:
        t = time.time()
        rng = np.random.default_rng(1000 + J)
        mins, dens = [], []
        R2base = R * R
        rho = smooth_orthant_density_Lambda(q2, R2base, 0.0)
        Wq = 2 * 400.0 / rho + 1e-9           # window wide enough to hold >=100 levels each side
        m = BPModel(q2, R2base + 1.0, Wq)       # enumerate once up to the largest shifted radius
        for sh in range(NSHIFT):
            R2 = float(R2base + rng.random())
            c = m.closest(R2, k=40)
            mins.append(c[0][0] / 2.0)        # |Lambda| = |Q-R^2|/2  (exact)
            Qw, _, _ = m.list_Q_window(R2 - Wq / 2, R2 + Wq / 2, max_points=10_000_000)
            lam = np.sort(np.abs((Qw - np.longdouble(R2)) / 2).astype(np.float64))
            k = min(100, len(lam))
            if k >= 20:
                # MLE of Poisson density from k smallest |Lambda| (two-sided): rho = k/(2 lam_k)
                dens.append(k / (2 * lam[k - 1]))
        mins = np.array(mins)
        rho_emp = np.mean(dens) if dens else float("nan")
        rows_c.append((J, R, R / math.sqrt(q2.min()), mins.mean(), 1 / (2 * rho), mins.mean() * 2 * rho,
                       rho_emp / rho, mins.min()))
        print(f"{J:3d} {R:7.1f} {R/math.sqrt(q2.min()):6.0f}  |  {mins.mean():12.4e}     {1/(2*rho):12.4e}  {mins.mean()*2*rho:6.3f} |"
              f"       {rho_emp/rho:7.4f}         | {mins.min():.3e}   [{time.time()-t:.1f}s, half {len(m.A):,}]", flush=True)
        del m
rows_c = np.array(rows_c)
print("\n fitted exponent of <min|Lambda|> vs R at fixed J (prediction from rho ~ R^{J-2}: -(J-2)):")
for J in MINCONF:
    s = rows_c[rows_c[:, 0] == J]
    p = np.polyfit(np.log(s[:, 1]), np.log(s[:, 3]), 1)[0]
    pr = np.polyfit(np.log(s[:, 1]), np.log(s[:, 4]), 1)[0]
    print(f"   J={J:2d}: fitted slope {p:+.3f}  (smooth-density prediction over same R range {pr:+.3f}; leading -(J-2) = {-(J-2)})")
OUT["min_scaling"] = rows_c.tolist()
print(" mean of ratio <min|Lam|> * 2 rho over all (J,R):", f"{rows_c[:,5].mean():.3f} +- {rows_c[:,5].std()/math.sqrt(len(rows_c)):.3f}"
      " (Poisson: 1; each entry averages only 6 shifts, so expect ~+-0.4 scatter per entry)")

hdr("PART (c'): fixed flux cutoff R = 11 (n_max <= 11), varying J")
for J in (4, 6, 8, 10, 12, 16, 20, 24):
    q2 = charges(J, 200 + J)
    R2 = 11.0 ** 2 + 0.3183098861837907
    rho = smooth_orthant_density_Lambda(q2, R2, 0.0)
    m = BPModel(q2, R2, 2 * 50.0 / rho + 1e-9)
    c = m.closest(R2, k=40)
    tot = smooth_orthant_count(q2, R2)
    print(f"  J={J:2d}: #orthant points inside shell ~ {tot:.3e};  rho(0) = {rho:.3e};  1/(2rho) = {1/(2*rho):.3e};"
          f"  exact min|Lambda| = {c[0][0]/2:.3e};  verified n = {c[0][2]}")
    del m

hdr("PART (d): sensitivity to the charges")
print("(d1) near-commensurate charges q_i^2 = 1 + delta*u_i, J = 8, R^2 = 400 + 0.3711 (generic)")
J = 8
R2 = 400.0 + 0.3711
print("  delta     delta*R^2 | rho_sm(0)  | N[|Lam|<1/rho_sm*50] /100 | max gap / mean spacing (|Lam|<1) | P(s<0.1) var(s) | min|Lam| exact")
for delta in (0.3, 3e-2, 3e-3, 3e-4, 3e-5, 0.0):
    q2 = charges(J, 300, "nearcomm", delta) if delta > 0 else np.ones(J)
    m = BPModel(q2, R2, 4.0)
    rho = smooth_orthant_density_Lambda(q2, R2, 0.0) if delta > 0 else smooth_orthant_density_Lambda(np.ones(J) * (1 + 1e-12), R2, 0.0)
    Qw, _, _ = m.list_Q_window(R2 - 2.0, R2 + 2.0, max_points=50_000_000)
    lam = ((Qw - np.longdouble(R2)) / 2).astype(np.float64)
    near = np.sum(np.abs(lam) < 50 / rho)
    # distinct values (exactly commensurate case has huge degeneracy)
    lu = np.unique(np.round(lam, 12))
    gaps = np.diff(lu)
    meansp = 1.0 / rho
    sp = gaps / gaps.mean() if len(gaps) > 1 else np.array([np.nan])
    c = m.closest(R2, k=10)
    print(f"  {delta:7.1e}  {delta*R2:9.3g} | {rho:10.3e} | {near/100:8.3f}                   | {gaps.max()/meansp:12.3g}"
          f"                    | {np.mean(sp<0.1):.3f}  {np.nanvar(sp):7.3f} | {c[0][0]/2:.3e}  (#distinct in |Lam|<1: {len(lu):,} vs pred {2*rho:.3g})")
    del m

print("\n(d2) repeated charges: J = 12 fluxes but only 3 distinct charges (each used 4 times)")
print("     prediction: distinct values are sum_a q_a^2 m_a, m_a in Z>=0 (Lagrange: every m_a is a sum of 4 squares)")
print("     -> distinct-value density per unit Lambda = X^2 / prod_a q_a^2 at X=R^2 (effective J_eff = 6), vs generic J=12.")
rng = np.random.default_rng(400)
qa = (1 + rng.random(3)) ** 2
q2 = np.repeat(qa, 4)
for R in (8.0, 11.0, 15.0):
    R2 = R * R + 0.271828
    m = BPModel(q2, R2, 4.0)
    Qw, _, _ = m.list_Q_window(R2 - 1.0, R2 + 1.0, max_points=80_000_000)
    lam = ((Qw - np.longdouble(R2)) / 2).astype(np.float64)
    lu = np.unique(np.round(lam, 10))
    dens_distinct = len(lu) / 1.0
    pred_eff = R2 ** 2 / np.prod(qa)
    rho12 = smooth_orthant_density_Lambda(q2 * (1 + 1e-9 * np.arange(12)), R2, 0.0)
    gaps = np.diff(lu)
    print(f"  R={R:5.1f}: orthant pts in |Lam|<0.5: {len(lam):,}; distinct Lambda values: {len(lu):,} (pred J_eff=6: {pred_eff:.4g});"
          f" mean degeneracy {len(lam)/len(lu):.1f}; generic-J=12 density would be {rho12:.3e};"
          f" min|Lam| = {np.min(np.abs(lu)):.3e} (1/(2 rho_eff) = {1/(2*pred_eff):.3e})")
    del m

hdr("PART (e): BP's own numbers -- what J and q are needed for Lambda_obs, under Poisson statistics")
# P(at least one vacuum with |Lambda| < eps) = 1 - exp(-2 rho eps), rho = BP density with D = 2^J
for J, logq in ((100, -1.6), (100, -1.43), (100, -1.73), (300, -0.87)):
    q = 10 ** logq
    lam_bare = 1.0
    R2 = 2 * lam_bare
    log_rho = (math.log10(2 * math.pi ** (J / 2) / math.gamma(J / 2)) + (J / 2 - 1) * math.log10(R2)
               - J * math.log10(2) - J * logq)
    print(f"  J={J}, q=10^{logq}: log10 rho(0) [distinct values per unit lambda, D=2^J] = {log_rho:.1f}"
          f";  expected # vacua with |lambda|<1e-120: 10^{log_rho-120+math.log10(2):.1f}")

json.dump(OUT, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "T3_A1_bp_discretuum_data.json"), "w"))
print(f"\nTotal runtime {time.time()-T0:.0f}s")
