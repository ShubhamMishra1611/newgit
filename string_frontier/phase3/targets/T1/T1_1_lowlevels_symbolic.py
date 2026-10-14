"""
T1 (i): exact low-level partial-wave coefficients of the intercept family
   A_a(s,t) = Gamma(-alpha(s)) Gamma(-alpha(t)) / Gamma(c - alpha(s) - alpha(t)),  alpha(x) = a + x,
c = 0 (seed family; c=0,a=1 is the bosonic string) and c = 1 (superstring-type family; a=0 is type I stripped).
Residue at alpha(s) = n (s_n = n - a):  R_n(t) = (alpha(t)+1)_{n-c} / n!   (overall sign (-1)^c dropped).
External scalars of mass m^2 = -a + delta (delta = 0: externals on the leading trajectory).
kappa_n = (s_n - 4 m^2)/2 = (n + 3a - 4 delta)/2;   t = kappa (x-1).
Unitarity: rho_j = c_j / kappa^j >= 0  (c_j = coefficient of C_j^{(D-3)/2}(x)).
We print rho_j *exactly* as rational functions of (a, D) [and delta], factored.
"""
import sympy as sp
a, D, x, dl = sp.symbols('a D x delta')
out = []
def log(*s):
    s = " ".join(str(z) for z in s); print(s); out.append(s)

def rhos(n, c, delta=0):
    kap = (n + 3*a - 4*delta) / sp.Integer(2)
    t = kap*(x - 1)
    R = sp.expand(sp.rf(a + t + 1, n - c) / sp.factorial(n))
    nu = (D - 3) / 2
    P = sp.Poly(R, x)
    rem = P
    res = {}
    for j in range(P.degree(), -1, -1):
        Cj = sp.Poly(sp.expand(sp.gegenbauer(j, nu, x)), x)
        cj = sp.cancel(rem.coeff_monomial(x**j) / Cj.coeff_monomial(x**j))
        rem = sp.Poly(sp.cancel(sp.expand(rem.as_expr() - Cj.as_expr()*cj)), x)
        res[j] = sp.factor(sp.cancel(cj / kap**j))
    assert sp.simplify(rem.as_expr()) == 0
    return res

for c in (0, 1):
    log(f"\n===== family c={c}, externals on trajectory (delta=0): rho_j(a,D) =====")
    for n in range(1 + c, 5 + c):
        r = rhos(n, c)
        for j in sorted(r):
            log(f"  level n={n}, spin j={j}: rho = {r[j]}")
log("\n===== c=0 family, general external mass m^2 = -a + delta, levels 1-3 =====")
for n in (1, 2, 3):
    r = rhos(n, 0, dl)
    for j in sorted(r):
        log(f"  level n={n}, spin j={j}: rho = {r[j]}")

# boundary curves (c=0, delta=0)
log("\n===== c=0, delta=0: boundary curves from levels 1-3 =====")
r2, r3 = rhos(2, 0), rhos(3, 0)
log("  level 1: rho_0 = (1-a)/2 >= 0  ->  a <= 1 ;  rho_1 > 0 for all a.")
for name, expr in (("L2 j=0", r2[0]), ("L3 j=0", r3[0]), ("L3 j=1", r3[1])):
    num = sp.factor(sp.numer(sp.together(expr)))
    sol = sp.solve(sp.Eq(num, 0), D)
    log(f"  {name}: numerator {num};  zero at D = {[sp.simplify(s) for s in sol]}")
with open(__file__.replace('.py', '_output.txt'), 'w') as f:
    f.write("\n".join(out) + "\n")
