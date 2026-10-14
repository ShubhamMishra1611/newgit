"""T3 (A1c): distribution and scaling of min|Lambda| in the BP model, with many shell radii.

Run:  python3 -I T3_A1c_min_statistics.py > T3_A1c_min_statistics.out   (~3-5 min)

For each (J, R) we draw NSH generic shell radii R^2 + u (u ~ U[0,1)), find the exact minimum
|Lambda| = |Q - R^2|/2 over all orthant lattice points (MITM search, then exact Fraction
re-evaluation), and normalise by the smooth density: z = 2 rho_sm(0) min|Lambda|.
Poisson (flat, gap-free) prediction: z ~ Exp(1) exactly, independent of J and R.
We test that with a KS test on the pooled z's, and fit the R-exponent of <min|Lambda|>.
"""
import sys, os, math, time
import numpy as np
from scipy import stats
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bp_lib import BPModel, smooth_orthant_density_Lambda

NSH = 40
CONF = {4: [300.0, 900.0, 2700.0], 6: [50.0, 100.0, 200.0], 8: [25.0, 40.0, 64.0],
        10: [16.0, 24.0, 36.0], 12: [12.0, 16.0, 22.0], 16: [8.0, 10.5, 14.0],
        20: [7.0, 8.5, 10.5], 24: [6.0, 7.5, 9.0]}
allz = []
print(" J     R   n_max | <min|Lam|>      1/(2rho)     <z>=<2 rho min>  +- s.e. | KS p (z~Exp(1)) | min|Lam| range")
fits = {}
for J, Rs in CONF.items():
    rng = np.random.default_rng(5000 + J)
    q = 1 + rng.random(J)
    q2 = q * q
    means = []
    for R in Rs:
        t = time.time()
        rho = smooth_orthant_density_Lambda(q2, R * R, 0.0)
        m = BPModel(q2, R * R + 1.0, 2 * 30.0 / rho + 1e-9)
        z = []
        mins = []
        for k in range(NSH):
            R2 = float(R * R + rng.random())
            rho_k = smooth_orthant_density_Lambda(q2, R2, 0.0)
            c = m.closest(R2, k=30)
            mn = c[0][0] / 2
            mins.append(mn)
            z.append(2 * rho_k * mn)
        z = np.array(z); mins = np.array(mins)
        allz.extend(z.tolist())
        ks = stats.kstest(z, "expon")
        means.append((R, mins.mean(), 1 / (2 * rho)))
        print(f"{J:3d} {R:7.1f} {R/q.min():5.0f} | {mins.mean():.4e}   {1/(2*rho):.4e}    {z.mean():.3f} +- {z.std()/math.sqrt(NSH):.3f}"
              f"     |   {ks.pvalue:.3f}         | [{mins.min():.2e}, {mins.max():.2e}]  ({time.time()-t:.1f}s)", flush=True)
        del m
    a = np.array(means)
    p_emp = np.polyfit(np.log(a[:, 0]), np.log(a[:, 1]), 1)[0]
    p_pred = np.polyfit(np.log(a[:, 0]), np.log(a[:, 2]), 1)[0]
    fits[J] = (p_emp, p_pred)
allz = np.array(allz)
ks = stats.kstest(allz, "expon")
print(f"\nPOOLED over all (J,R): N={len(allz)}, <z> = {allz.mean():.4f} +- {allz.std()/math.sqrt(len(allz)):.4f} (Exp(1): 1),"
      f" var z = {allz.var():.3f} (1), KS D = {ks.statistic:.4f}, p = {ks.pvalue:.3f}")
print("R-exponent of <min|Lambda|> at fixed J  (measured | smooth-density prediction over same R | asymptotic -(J-2)):")
for J, (pe, pp) in fits.items():
    print(f"   J={J:2d}: {pe:+.2f} | {pp:+.2f} | {-(J-2)}")
