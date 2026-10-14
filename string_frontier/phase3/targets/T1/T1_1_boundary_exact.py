"""
T1 (i): exact piecewise boundary D_max(a) of the c=0 intercept family on the on-trajectory slice.
Step 1: dense scan a in (-21/29, 1], step 1/5000, levels n <= 8 (rigorous per-level roots) -> which
        (level, spin) is binding where.  (Levels n >= 6 never bind: see T1_1_highlevel_check output.)
Step 2: exact transition points between consecutive binding pieces via resultants in D.
"""
import sympy as sp
from math import factorial
from flint import fmpq
from pwcore import poch_poly, level_crit, D_of_nu
out = []
def log(*s):
    s = " ".join(str(z) for z in s); print(s, flush=True); out.append(s)

NL = 8
pieces = []
k0, k1, den = -3620, 5000, 5000    # a from -0.724 to 1
prev = None
for k in range(k0, k1 + 1):
    a = fmpq(k, den)
    best = None
    for n in range(1, NL + 1):
        nu, j = level_crit(poch_poly(a + 1, n) * fmpq(1, factorial(n)), (n + 3 * a) / 2)
        if nu is None: continue
        D = float(D_of_nu(nu).mid())
        if best is None or D < best[0]: best = (D, n, j)
    key = (best[1], best[2]) if best[0] > 3 else ('none', '-')
    if key != prev:
        pieces.append([key, k, k, best[0], best[0]])
        prev = key
    else:
        pieces[-1][2] = k; pieces[-1][4] = best[0]
log(f"binding (level, spin) pieces on a-grid step 1/{den}, levels <= {NL}:")
for key, ka, kb, Da, Db in pieces:
    log(f"  {key}: a in [{ka/den:.4f}, {kb/den:.4f}]   D_max from {Da:.6f} to {Db:.6f}")

# exact symbolic boundary polynomials (numerators of rho_0 at levels 2,3,4,5, from T1_1_lowlevels_symbolic / T1_1_transitions)
a, D = sp.symbols('a D')
N2 = D*a**2 - 2*D*a + 8*a**2 + 14*a + 4
N3 = D*a - 3*D + 26*a + 30                     # times -(a-1)(a+1)
N4 = (D**2*a**4 - 4*D**2*a**3 - 4*D**2*a**2 + 16*D**2*a + 54*D*a**4 + 36*D*a**3 - 228*D*a**2 - 288*D*a - 64*D
      + 296*a**4 + 1336*a**3 + 2368*a**2 + 2000*a + 704)
N5 = (D**2*a**4 - 4*D**2*a**3 - 14*D**2*a**2 + 36*D**2*a + 45*D**2 + 90*D*a**4 + 120*D*a**3 - 800*D*a**2
      - 2000*D*a - 1250*D + 1304*a**4 + 8224*a**3 + 19464*a**2 + 20464*a + 8080)   # times -(a-1)
polys = {(2, 0): N2, (3, 0): N3, (4, 0): N4, (5, 0): N5}
log("\nexact boundary pieces: D_max(a) = smallest root D>3 of the listed polynomial (rho_0 numerator):")
for k, p in polys.items():
    log(f"  level {k[0]} spin 0: {sp.factor(p)} = 0")
seq = [pc[0] for pc in pieces if pc[0][0] != 'none']
log("\ntransition points (common roots in (a,D)):")
for k1_, k2_ in zip(seq, seq[1:]):
    if k1_ not in polys or k2_ not in polys:
        log(f"  {k1_}->{k2_}: polynomial not tabulated"); continue
    res = sp.factor(sp.resultant(polys[k1_], polys[k2_], D))
    # locate the root near the numerical transition
    kk = [pc for pc in pieces if pc[0] == k2_][0][1] / den
    fac = [f for f, m in sp.factor_list(res)[1] if sp.Poly(f, a).degree() > 0]
    best = None
    for f in fac:
        for r in sp.Poly(f, a).nroots(n=25):
            if abs(sp.im(r)) < 1e-20 and (best is None or abs(sp.re(r) - kk) < abs(best[0] - kk)):
                best = (sp.re(r), f)
    ar = best[0]
    Dvals = [sp.re(s) for s in sp.Poly(polys[k1_].subs(a, ar), D).nroots(n=25) if abs(sp.im(s)) < 1e-15 and sp.re(s) > 3]
    log(f"  {k1_} -> {k2_}: a* = {sp.N(ar, 18)} (root of {best[1]}),  D(a*) = {sp.N(min(Dvals), 15) if Dvals else '?'}")
with open(__file__.replace('.py', '_output.txt'), 'w') as f:
    f.write("\n".join(out) + "\n")
