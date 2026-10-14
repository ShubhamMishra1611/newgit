import sympy as sp
x, D, a = sp.symbols('x D a')
def coeffs(n, aval):
    tv = -(n + 3*aval)/sp.Integer(2)*(1 - x)
    R = sp.Poly(sp.expand(sp.rf(aval + tv + 1, n)/sp.factorial(n)), x)
    lam = (D-3)/2
    rem = R; cs = {}
    for j in range(n, -1, -1):
        Cj = sp.Poly(sp.expand(sp.gegenbauer(j, lam, x)), x)
        cj = sp.cancel(rem.coeff_monomial(x**j)/Cj.coeff_monomial(x**j))
        cs[j] = sp.factor(cj)
        rem = sp.Poly(sp.cancel(sp.expand(rem.as_expr() - Cj.as_expr()*cj)), x)
    return cs
# symbolic in a and D at level 2,3
for n in (2, 3):
    cs = coeffs(n, a)
    for j, c in cs.items():
        print(n, j, c)
