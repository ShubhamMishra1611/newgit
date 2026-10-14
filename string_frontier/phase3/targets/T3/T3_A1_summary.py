"""T3 (A1) post-processing: flat-prior test on DISJOINT shells Lambda in [eps_k, eps_{k+1})
(nested windows in the main output are correlated).  Reads T3_A1_bp_discretuum_data.json.
Run: python3 -I T3_A1_summary.py > T3_A1_summary.out
"""
import json, os, math
import numpy as np
from scipy import stats

D = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "T3_A1_bp_discretuum_data.json")))
print(" J | disjoint shells (Lam>0 and Lam<0) | chi2/dof of (n-np)/sqrt(np) | p-value | max|pull| | max|n/np-1| for np>=1e4")
for J in (4, 6, 8, 10, 12):
    d = D[f"J{J}"]
    rows = np.array(d["Neps"])  # eps, N+, Np+, N-, Np-
    pulls = []
    rel = []
    for col_n, col_p in ((1, 2), (3, 4)):
        n = np.diff(np.concatenate([[0], rows[:, col_n]]))
        p = np.diff(np.concatenate([[0], rows[:, col_p]]))
        ok = p > 5
        pulls.extend(((n[ok] - p[ok]) / np.sqrt(p[ok])).tolist())
        big = p >= 1e4
        rel.extend(np.abs(n[big] / p[big] - 1).tolist())
    pulls = np.array(pulls)
    chi2 = float(np.sum(pulls ** 2))
    dof = len(pulls)
    print(f"{J:3d}|  {dof:3d}                              |  {chi2/dof:6.2f}                     | {stats.chi2.sf(chi2, dof):.3f}   |"
          f"  {np.abs(pulls).max():.2f}    |  {max(rel):.2e}")
print("\nNB the two outermost shells (eps > 1e-3 R^2) test the curvature of the smooth density, not flatness;"
      " the smooth formula includes it.")
