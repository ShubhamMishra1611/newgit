"""RT3: is the mild chi^2 excess of the flatness test systematic? extra seeds at J=4 (R=800) and J=10 (R=40)."""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rt3_bp2_flatness_power import flatness
from scipy import stats
for J, R in ((4, 800.0), (10, 40.0), (8, 70.0)):
    chis = []
    for seed in range(500, 506):
        rng = np.random.default_rng(seed + 10 * J)
        c = (1 + rng.random(J)) ** 2
        r = flatness(c, R * R + 0.3183, nlev=50_000, verbose=False)
        chis.append(r['chi2_dof'])
        print(f"  J={J} seed={seed+10*J}: p={r['p']:.4f} chi2/dof={r['chi2_dof']:.2f} (p={r['chi2_p']:.3f}) maxpull={r['maxpull']:.2f}"
              f" var_s={r['var_s']:.3f} KS p={r['ks_p']:.3f}", flush=True)
    tot = np.sum(np.array(chis) * 44); dof = 44 * len(chis)
    print(f"  J={J}: pooled chi2/dof = {tot/dof:.3f} over {dof} dof, p = {stats.chi2.sf(tot, dof):.3g}")
