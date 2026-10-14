"""Exact transition points of D_max(a) for c=0 on-trajectory (from my own symbolic rho's)."""
import sympy as sp
a, D = sp.symbols('a D')
N2 = D*(a**2-2*a) + 8*a**2 + 14*a + 4                     # level 2 spin 0 numerator
N3 = (26*a+30) - D*(3-a)                                   # level 3 spin 0 (times (1-a)(1+a))
N4 = (D**2*a**4 - 4*D**2*a**3 - 4*D**2*a**2 + 16*D**2*a + 54*D*a**4 + 36*D*a**3 - 228*D*a**2 - 288*D*a - 64*D
      + 296*a**4 + 1336*a**3 + 2368*a**2 + 2000*a + 704)
N5 = (D**2*a**4 - 4*D**2*a**3 - 14*D**2*a**2 + 36*D**2*a + 45*D**2 + 90*D*a**4 + 120*D*a**3 - 800*D*a**2
      - 2000*D*a - 1250*D + 1304*a**4 + 8224*a**3 + 19464*a**2 + 20464*a + 8080)   # times -(a-1)
D3 = sp.solve(N3, D)[0]; D2 = sp.solve(N2, D)[0]
print('D3 =', D3, ' D2 =', sp.factor(D2))
r35 = sp.factor(sp.numer(sp.together(N5.subs(D, D3))))
print('level3 vs level5 :', r35)
r24 = sp.factor(sp.numer(sp.together(N4.subs(D, D2))))
print('level2 vs level4 :', r24)
r45 = sp.factor(sp.resultant(N4, N5, D))
print('level4 vs level5 resultant factors:', [ (sp.degree(f[0],a), f[1]) for f in sp.factor_list(r45)[1]])
def roots_in(p, lo, hi):
    return [r for r in sp.Poly(p, a).nroots(n=30) if r.is_real and lo < r < hi]
for name, p in [('3/5', r35), ('2/4', r24)]:
    for f, m in sp.factor_list(p)[1]:
        rr = roots_in(f, 0.3, 0.6)
        if rr: print(name, f, rr)
for f, m in sp.factor_list(r45)[1]:
    rr = roots_in(f, 0.3, 0.6)
    if rr: print('4/5 factor deg', sp.degree(f, a), rr)
# D at transitions; D_max values
def smallest_root_gt3(expr, aval):
    pr = sp.Poly(expr.subs(a, aval), D)
    rts = sorted([r for r in pr.nroots(n=30) if abs(sp.im(r)) < 1e-25 and sp.re(r) > 3], key=lambda z: sp.re(z))
    return rts[0] if rts else None
for av in [sp.Rational(-1,3), 0, sp.Rational(-3,5)]:
    print('D3 at', av, '=', D3.subs(a, av))
print('D3=3 at a =', sp.solve(sp.Eq(D3, 3), a))
