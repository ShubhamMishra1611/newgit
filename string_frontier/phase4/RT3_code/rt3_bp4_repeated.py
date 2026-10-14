"""RT3 check of T3 claim 5b (repeated charges => J_eff).  J=12 fluxes, 3 distinct charges c_a used 4 times each.
Independent count: distinct Lambda values = #{m in Z^3_{>=0} : sum c_a m_a in [S-1, S+1)} (Lagrange four-square:
every m_a >= 0 is a sum of 4 squares), compared with my MITM orthant enumeration (distinct values by exact
integer-vector key) and with the J_eff=6 density X^2/prod c_a."""
import sys, os, math
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rt3_bp_core import Shell
rng = np.random.default_rng(400)       # T3's recipe for the 3 charges (re-implemented)
ca = (1 + rng.random(3)) ** 2
c = np.repeat(ca, 4)
for R in (8.0, 11.0):
    S = R * R + 0.271828
    # independent: simplex lattice count
    m1 = np.arange(int((S + 1) / ca[0]) + 1)
    cnt = 0
    for a in m1:
        for b in range(int((S + 1 - ca[0] * a) / ca[1]) + 1):
            base = ca[0] * a + ca[1] * b
            lo = math.ceil(max(S - 1 - base, 0) / ca[2] - 1e-12); hi = math.floor((S + 1 - base) / ca[2] - 1e-12)
            # count m3 with base + ca[2] m3 in [S-1, S+1)
            m3 = np.arange(max(lo - 1, 0), hi + 2)
            v = base + ca[2] * m3
            cnt += int(np.sum((v >= S - 1) & (v < S + 1)))
    # MITM orthant enumeration: distinct values via the (m_a) key
    sh = Shell(c, S + 1.0)
    X, ia, ib = sh.listing(S - 1.0, S + 1.0, with_n=True)
    na = sh.na[ia].astype(np.int64); nb = sh.nb[ib].astype(np.int64)
    n = np.hstack([na, nb])
    m = np.stack([(n[:, 4*a_:4*a_+4] ** 2).sum(1) for a_ in range(3)], 1)   # m_a = sum of the 4 squares with charge c_a
    keys = np.unique(m, axis=0)
    pred = S ** 2 / np.prod(ca)
    print(f"R={R}: orthant points in |Lam|<0.5: {len(X):,}; distinct (by m-key): {len(keys):,}; simplex count (independent): {cnt:,};"
          f" J_eff=6 prediction X^2/prod c = {pred:.1f}; mean degeneracy {len(X)/len(keys):.0f}")
print("T3: R=8: 271,373 pts, 315 distinct (pred 282.8), deg 861.5;  R=11: 3,473,478 pts, 1,059 distinct (pred 1007), deg 3280")
