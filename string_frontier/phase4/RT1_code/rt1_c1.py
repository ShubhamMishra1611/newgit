"""Claim 5: c=1 family on trajectory (m^2=-a): residue (al_t+1)_{n-1}.
mode sel  : D_max at selected a from levels<=NLOW, check levels up to NHIGH at D_max
mode empty: for a in [-6, -7/11) check whether any D>3 survives levels<=N (test at D=3+1e-6)
mode upper: for a in (0, 2] check level-2 spin-0 kills (rho_0=-a/4)"""
import sys; sys.path.insert(0, '.')
from rt1core import *
from fractions import Fraction as Fr
mode = sys.argv[1]
if mode == 'sel':
    NLOW, NHIGH = int(sys.argv[2]), int(sys.argv[3])
    for a in [Fr(s) for s in sys.argv[4].split(',')]:
        res = []
        for n in range(2, NLOW + 1):
            f, k = open_facs(n, a, -a, 1)
            if k == 0: continue
            d, bad = Dcrit_factors(f, k, tol=Fr(1, 10**10))
            if d is not None: res.append((d, n, bad))
        res.sort(); Dm = res[0]
        worst = None
        for n in range(NLOW + 1, NHIGH + 1):
            f, k = open_facs(n, a, -a, 1)
            if k == 0: continue
            ok, bad = ok_factors(f, k, Dm[0])
            if not ok: worst = (n, bad); break
        form = (8 * a + 10) / (1 - a)
        print(f"c=1 a={str(a):>7} D_max={float(Dm[0]):.10f} [level {Dm[1]} spins {Dm[2]}] (8a+10)/(1-a)={float(form):.10f} "
              f"next {[(r[1], round(float(r[0]),4)) for r in res[1:3]]} | levels {NLOW+1}..{NHIGH}: {'ALL POSITIVE' if worst is None else 'FAIL '+str(worst)}", flush=True)
elif mode == 'empty':
    N = int(sys.argv[2]); step = Fr(1, int(sys.argv[3]))
    D0 = Fr(3) + Fr(1, 10**6)
    a = Fr(-6); surv = []; rows = []
    while a < Fr(-7, 11):
        fails = set()
        for n in range(2, N + 1):
            f, k = open_facs(n, a, -a, 1)
            if k == 0: continue
            ok, bad = ok_factors(f, k, D0)
            if not ok: fails.add(n)
        if not fails: surv.append(a)
        rows.append((a, fails)); a += step
    gaps = [(float(rows[i][0]), float(rows[i+1][0])) for i in range(len(rows)-1) if not (rows[i][1] & rows[i+1][1])]
    print(f"c=1 empty-check N={N} step={step}: points={len(rows)} survivors={[float(s) for s in surv]} gaps={gaps}")
