"""Claim 3 (parity lemma) and the closed-string explanation for claim 6.
(A) Random test: for definite-parity R, Gegenbauer-positivity of R  <=>  of (1+x)R  (several nu).
(B) Corollary: Dcrit(n) of c=0,a=0 (m^2=0) == Dcrit(n) of type-I stripped (c=1,a=0,m^2=0), n=2..60.
(C) Why the closed 4-tachyon amplitude cannot see D=26: residue = (open residue)^2, so the spin-0
    coefficient is c_0 = <R^2>/<1> > 0 identically; explicit check at D=27 at the closed level that
    contains the open level-2 scalar ghost."""
import sys, random; sys.path.insert(0, '.')
from rt1core import *
from fractions import Fraction as Fr
random.seed(1)
def times_x(c, nu):
    new = [Q(0)] * (len(c) + 1)
    for j, cj in enumerate(c):
        if cj == 0: continue
        d = 2 * (j + nu)
        new[j + 1] += cj * (j + 1) / d
        if j >= 1: new[j - 1] += cj * (j + 2 * nu - 1) / d
    return new
def monomials_to_geg(p, nu):
    out = [Q(0)] * len(p); cur = [Q(1)]
    for k, pk in enumerate(p):
        if k > 0: cur = times_x(cur, nu)
        for j, v in enumerate(cur): out[j] += pk * v
    return out
cnt = agree = 0; npos = 0
for trial in range(3000):
    deg = random.randint(1, 9); par = random.randint(0, 1)
    p = [Q(0)] * (deg + 1)
    for k in range(deg + 1):
        if k % 2 == par % 2: p[k] = Q(random.randint(-6, 9), random.randint(1, 5))
    if all(v == 0 for v in p): continue
    p[-1] = abs(p[-1]) + 1 if deg % 2 == par % 2 else p[-1]
    for nu in (Q(1, 4), Q(1, 2), Q(7, 2), Q(23, 2)):
        cR = monomials_to_geg(p, nu)
        q1 = [Q(0)] * (len(p) + 1)
        for k, v in enumerate(p): q1[k] += v; q1[k + 1] += v
        cS = monomials_to_geg(q1, nu)
        a = all(v >= 0 for v in cR); b = all(v >= 0 for v in cS)
        cnt += 1; agree += (a == b); npos += a
print(f"(A) parity lemma random test: {agree}/{cnt} agree (R positive in {npos} cases)")
mism = []
for n in range(2, 61):
    f0, k0 = open_facs(n, 0, 0, 0); f1, k1 = open_facs(n, 0, 0, 1)
    d0, _ = Dcrit_factors(f0, k0, tol=Fr(1, 10**9)); d1, _ = Dcrit_factors(f1, k1, tol=Fr(1, 10**9))
    if (d0 is None) != (d1 is None) or (d0 is not None and abs(d0 - d1) > Fr(1, 10**8)): mism.append((n, d0, d1))
print(f"(B) c=0,a=0 vs type-I Dcrit(n), n=2..60: mismatches = {mism}")
# (C)
fo, ko = open_facs(2, 1, -1, 0)
print("(C) open bosonic level-2 rho_j at D=27:", [float(v) for v in exact_rhos_factors(fo, ko, 27)])
fc, kc = closed_facs(2, 1, -1, -1)
print("    closed tachyon level-2 (=open^2) rho_j at D=27:", [float(v) for v in exact_rhos_factors(fc, kc, 27)])
for D in (27, 40, 57, 58):
    worst = None
    for n in range(1, 31):
        f, k = closed_facs(n, 1, -1, -1)
        ok, bad = ok_factors(f, k, D)
        if not ok: worst = (n, bad); break
    print(f"    closed bosonic tachyon, levels<=30 at D={D}: {'all positive' if worst is None else 'first failure '+str(worst)}")
