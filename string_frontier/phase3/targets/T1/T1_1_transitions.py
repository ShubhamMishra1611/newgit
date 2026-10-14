"""
T1 (i): exact boundary D_max(a) of the c=0 intercept family (delta=0) from levels 2, 3, 5.
Finds the level-5 spin-0 boundary curve and the transition points between the
level-3 (spin 0), level-5 (spin 0) and level-2 (spin 0) pieces.
"""
import sympy as sp
a, D, x = sp.symbols('a D x')
out = []
def log(*s):
    s = " ".join(str(z) for z in s); print(s); out.append(s)
def rho(n, j, c=0):
    kap = (n + 3*a) / sp.Integer(2)
    R = sp.Poly(sp.expand(sp.rf(a + kap*(x - 1) + 1, n - c) / sp.factorial(n)), x)
    nu = (D - 3) / 2
    rem = R; res = {}
    for jj in range(R.degree(), -1, -1):
        Cj = sp.Poly(sp.expand(sp.gegenbauer(jj, nu, x)), x)
        cj = sp.cancel(rem.coeff_monomial(x**jj) / Cj.coeff_monomial(x**jj))
        rem = sp.Poly(sp.cancel(sp.expand(rem.as_expr() - Cj.as_expr()*cj)), x)
        res[jj] = sp.factor(sp.cancel(cj / kap**jj))
    return res[j]
L2 = sp.numer(sp.together(rho(2, 0)))
L3 = sp.numer(sp.together(rho(3, 0)))
r5 = rho(5, 0)
L5 = sp.factor(sp.numer(sp.together(r5)))
log("level-5 spin-0 rho (c=0):", r5)
D2 = sp.solve(L2, D)[0]
D3 = sp.solve(sp.Poly(L3, D).as_expr() / ((1 - a)*(1 + a)), D)
D3 = [s for s in sp.solve(sp.cancel(L3 / ((a - 1)*(a + 1))), D)][0]
log("D_L2(a) =", sp.simplify(D2), "   D_L3(a) =", sp.simplify(D3))
# level 5 boundary: quadratic in D
P5 = sp.Poly(sp.cancel(L5 / sp.gcd(L5, (a-1)*(a+1)*(a+2)*(a+3))), D)
log("level-5 spin-0 numerator factorised:", L5)
log("  as polynomial in D:", P5.as_expr())
# transitions
f = sp.lambdify(a, sp.Poly(L5, D).as_expr(), 'mpmath')
import mpmath as mp
mp.mp.dps = 30
def L5at(av, Dv):
    return sp.N(L5.subs({a: av, D: Dv}), 30)
def D5(av):
    # smallest root D>3 of level-5 numerator where it turns negative
    poly = sp.Poly(L5.subs(a, av), D)
    rts = [r for r in sp.Poly(poly, D).nroots(n=30) if abs(sp.im(r)) < 1e-25 and sp.re(r) > 3]
    return min(sp.re(r) for r in rts) if rts else sp.oo
g1 = lambda av: D5(av) - D3.subs(a, av)
g2 = lambda av: D5(av) - D2.subs(a, av)
a1 = mp.findroot(lambda z: g1(sp.Float(str(z), 30)), 0.37)
a2 = mp.findroot(lambda z: g2(sp.Float(str(z), 30)), 0.53)
log(f"transition L3(j=0) -> L5(j=0) at a1 = {mp.nstr(a1, 20)},  D = {sp.N(D3.subs(a, sp.Float(str(a1),30)), 15)}")
log(f"transition L5(j=0) -> L2(j=0) at a2 = {mp.nstr(a2, 20)},  D = {sp.N(D2.subs(a, sp.Float(str(a2),30)), 15)}")
# exact algebraic characterisation: at a transition both numerators vanish -> resultant in D
res1 = sp.factor(sp.resultant(sp.cancel(L3/((a-1)*(a+1))), L5, D))
res2 = sp.factor(sp.resultant(L2, L5, D))
log("resultant_D(L3', L5) =", res1)
log("resultant_D(L2, L5)  =", res2)
for av in (sp.Rational(2, 5), sp.Rational(9, 20), sp.Rational(1, 2)):
    log(f"  check a={av}: D_L5 = {sp.N(D5(av), 12)}, D_L3 = {sp.N(D3.subs(a, av), 12)}, D_L2 = {sp.N(D2.subs(a, av), 12)}")
# left end of the region (D_max -> 3, and D_max = 4)
log("D_L3(a)=3 at a =", sp.solve(sp.Eq(D3, 3), a), ";  D_L3(a)=4 at a =", sp.solve(sp.Eq(D3, 4), a))
# monotonicity of the pieces
log("dD_L2/da =", sp.factor(sp.diff(D2, a)), ";  dD_L3/da =", sp.factor(sp.diff(D3, a)))
with open(__file__.replace('.py', '_output.txt'), 'w') as fh:
    fh.write("\n".join(out) + "\n")
