"""RT3 diagnostic: where do the chi^2 excesses of the flatness test come from? (per-shell pulls)"""
import sys, os, math
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rt3_bp_core import Shell, smooth_count, smooth_rho
for J, R, seed in ((4, 800.0, 101), (10, 40.0, 104)):
    rng = np.random.default_rng(seed); c = (1 + rng.random(J)) ** 2; S = R * R + 0.3183
    rho = smooth_rho(c, S); epsmax = 1e-3 * S
    sh = Shell(c, S + 2 * epsmax + 1e-9)
    eps = np.geomspace(3.0 / rho, epsmax, 22)
    print(f"J={J} R={R}")
    prevp = prevm = 0; pp0 = pm0 = 0.0
    for e in eps:
        Np = sh.count(S, S + 2 * e); Nm = sh.count(S - 2 * e, S)
        Pp = float(smooth_count(c, S + 2 * e) - smooth_count(c, S)); Pm = float(smooth_count(c, S) - smooth_count(c, S - 2 * e))
        np_, nm_ = Np - prevp, Nm - prevm; pp_, pm_ = Pp - pp0, Pm - pm0
        print(f"  eps={e:9.3e} shell+: n={np_:10d} pred={pp_:12.1f} pull={(np_-pp_)/math.sqrt(pp_):+6.2f} rel={np_/pp_-1:+.2e} | shell-: pull={(nm_-pm_)/math.sqrt(pm_):+6.2f} rel={nm_/pm_-1:+.2e}")
        prevp, prevm, pp0, pm0 = Np, Nm, Pp, Pm
