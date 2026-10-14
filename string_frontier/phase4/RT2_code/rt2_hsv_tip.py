"""2D scan of the HSV curvature-corrected KPV potential to locate the shark-fin tip (validation item of claim 6)."""
import numpy as np, sys
sys.path.insert(0, '/home/user/newgit/string_frontier/phase4/RT2_code')
from rt2_kpv import v_hsv, extrema
best = (None, 1e9)
for x in np.arange(14.0, 21.01, 0.1):
    for q in np.arange(0.0780, 0.0860, 0.00025):
        m, _ = extrema(lambda t: v_hsv(t, q, x))
        if any(t < np.pi/2 for t, _ in m):
            if x < best[1]: best = (q, x)
    if best[0] is not None: break
print('HSV shark-fin tip (2D scan, dg_sM=0.1, dp/M=0.00025): p/M=%.5f g_sM=%.2f  [HSV: 0.0826, 20.363]' % best)
