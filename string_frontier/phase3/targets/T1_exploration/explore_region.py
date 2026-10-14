import numpy as np
from numpy.polynomial import polynomial as P
from scipy.special import eval_gegenbauer, roots_gegenbauer, gammaln
import math
# Gegenbauer projection via Gauss-Gegenbauer quadrature (exact for polynomial integrands)
def pw(n, a, D):
    lam = (D-3)/2
    xs, ws = roots_gegenbauer(n+2, lam)
    t = -(n+3*a)/2*(1-xs)
    R = np.ones_like(xs)
    for k in range(n):
        R *= (a + t + 1 + k)
    R /= math.factorial(n)
    cs = []
    for j in range(n+1):
        Cj = eval_gegenbauer(j, lam, xs)
        norm = np.sum(ws*Cj*Cj)
        cs.append(np.sum(ws*R*Cj)/norm)
    return np.array(cs)
def ok(a, D, N):
    for n in range(1, N+1):
        c = pw(n, a, D)
        scale = np.max(np.abs(c))
        if np.any(c < -1e-10*scale):
            return False, n
    return True, None
for N in (2, 6, 12, 20):
    best = (0, None)
    rows = []
    for a in np.linspace(-0.6, 1.0, 33):
        lo, hi = 3.5, 60.0
        if not ok(a, lo, N)[0]:
            rows.append((round(a,3), None)); continue
        for _ in range(40):
            mid = (lo+hi)/2
            if ok(a, mid, N)[0]: lo = mid
            else: hi = mid
        rows.append((round(a,3), round(lo,3)))
        if lo > best[0]: best = (lo, a)
    print(f"N={N}: Dmax(a) =", rows[::4], " argmax:", best)
# just above a = 1
for a in (1.0, 1.001, 1.01):
    print(a, [ok(a, D, 12) for D in (4.0, 10.0, 26.0)])
