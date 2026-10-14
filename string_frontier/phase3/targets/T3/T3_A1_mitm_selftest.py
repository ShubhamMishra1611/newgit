"""T3 (A1) self-test: MITM orthant enumeration (bp_lib) vs brute force for J=3,4,5; smooth face-corrected count vs exact."""
import sys, itertools, math, numpy as np
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bp_lib import *
rng = np.random.default_rng(0)
for J in (3,4,5):
    q2 = (1+rng.random(J))**2
    R2 = 150.0
    m = BPModel(q2, R2, 30.0)
    # brute force
    nmax = [int(math.sqrt((R2+30)/c))+1 for c in q2]
    vals = []
    for n in itertools.product(*[range(k+1) for k in nmax]):
        Q = sum(c*k*k for c,k in zip(q2,n))
        if Q <= R2+30: vals.append(Q)
    vals = np.sort(np.array(vals))
    for lo,hi in ((R2-5,R2+5),(R2,R2+0.3),(R2-20,R2+20)):
        bf = int(np.sum((vals>=lo)&(vals<hi)))
        print(J, lo, hi, 'mitm', m.count_Q_window(lo,hi), 'brute', bf)
    Q, ia, ib = m.list_Q_window(R2-5, R2+5)
    bfw = vals[(vals>=R2-5)&(vals<R2+5)]
    print('  list equal:', np.allclose(Q, bfw), len(Q), len(bfw))
    c = m.closest(R2, k=5)
    j = np.argmin(np.abs(vals-R2)); print('  closest exact', c[0][:2], 'brute', vals[j]-R2)
    print('  smooth count N(R2+30)', smooth_orthant_count(q2, R2+30), 'actual', len(vals))
