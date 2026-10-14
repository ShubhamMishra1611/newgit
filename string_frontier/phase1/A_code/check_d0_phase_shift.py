"""One-loop D0-D0 eikonal phase shift in BFSS (SU(2) background X^1 = v t sigma3/2, X^2 = b sigma3/2).
Off-diagonal fluctuation 'mass^2' = r^2 + c v with r^2=b^2+v^2 t^2 (Euclidean: v -> i v gives cosh).
Weighted mode sum (bosons 8 x c=0, 2 x c=+-2; ghosts -2 x c=0; fermions -8 x c=+-1):
    f(x) = 6 + 2 cosh 2x - 8 cosh x ,  x = s v
Effective action  Gamma(b,v) = int_0^inf ds/s e^{-b^2 s} f(sv)/(2 sinh sv).
Compare term-by-term with  -int dtau V(r),  V = -sum_k c_k v^{2k+2}/r^{4k+3}  (one loop).
"""
import sympy as sp
x, s, v, b, u = sp.symbols('x s v b u', positive=True)
f = 6 + 2*sp.cosh(2*x) - 8*sp.cosh(x)
g = sp.series(f/(2*sp.sinh(x)), x, 0, 12).removeO()
print("f(x)/(2 sinh x) =", sp.expand(g))
results = {}
for n in [3, 5, 7, 9]:
    a_n = g.coeff(x, n)                     # coefficient of (sv)^n
    # int ds/s e^{-b^2 s} s^n v^n = Gamma(n) v^n / b^{2n}
    Gam_term = a_n*sp.gamma(n)*v**n/b**(2*n)
    # V_k = C v^{n+1}/r^{2n+1};  int dtau v^{n+1}/(b^2+v^2 tau^2)^{(2n+1)/2} = v^n/b^{2n} * I_n
    I_n = sp.integrate((1+u**2)**(-sp.Rational(2*n+1, 2)), (u, -sp.oo, sp.oo))
    C = sp.simplify(a_n*sp.gamma(n)/I_n)
    results[n] = C
    print(f"v^{n+1}/r^{2*n+1}: series coeff a_{n} = {a_n},  potential coefficient C = {C}")
assert results[3] == sp.Rational(15, 16), results[3]
assert results[5] == 0
print("CHECK: v^4/r^7 coefficient = 15/16 (BBPT hep-th/9706072 eq.13) and v^6/r^11 coefficient = 0 (BBPT remark) -- PASSED")
