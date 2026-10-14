"""FRONT E sanity checks (Phase 1). Run: python3 -I checks.py
(a) low-energy coefficients of Veneziano / Virasoro-Shapiro (zeta values)
(b) level-truncation identity, eq.(5) of arXiv:2406.02665
(c) LP reproducing the 3-sum-rule bound g3~ >= -10.61249 of arXiv:2011.02957 eq.(3.10)-(3.11), d=4
"""
import mpmath as mp
import numpy as np
from scipy.optimize import linprog
import sympy as sp

mp.mp.dps = 30
print("=== (a) zeta coefficients ===")
s, t = sp.symbols('s t')
# open: Gamma(1-s)Gamma(1-t)/Gamma(1-s-t) = exp(sum_k zeta(k)/k (s^k+t^k-(s+t)^k))
N = 6
expo = sum(sp.zeta(k)/k*(s**k + t**k - (s+t)**k) for k in range(2, N))
ser = sp.expand(sp.series(sp.exp(sp.Symbol('e')), sp.Symbol('e'), 0, 4).removeO().subs(sp.Symbol('e'), expo))
poly = sp.Poly(ser, s, t)
c11 = poly.coeff_monomial(s*t); c21 = poly.coeff_monomial(s**2*t); c31 = poly.coeff_monomial(s**3*t)
print("coef st      =", c11, "(expect -zeta2)")
print("coef s^2 t   =", c21, "(expect -zeta3)")
print("coef s^3 t   =", sp.nsimplify(c31), "=", sp.N(c31))
print("zeta3/zeta2 =", mp.zeta(3)/mp.zeta(2), " zeta4/zeta2 =", mp.zeta(4)/mp.zeta(2))
# numeric check of the exp identity
sv, tv = mp.mpf('0.13'), mp.mpf('-0.21')
lhs = mp.gamma(1-sv)*mp.gamma(1-tv)/mp.gamma(1-sv-tv)
rhs = mp.e**sum(mp.zeta(k)/k*(sv**k+tv**k-(sv+tv)**k) for k in range(2, 80))
print("exp-identity check:", lhs, rhs)
# closed (VS, SUSY-stripped): -G(-s)G(-t)G(-u)/(G(1+s)G(1+t)G(1+u)) = 1/(stu) + 2 zeta3 + ...
def VS(sv, tv):
    uv = -sv - tv
    return -mp.gamma(-sv)*mp.gamma(-tv)*mp.gamma(-uv)/(mp.gamma(1+sv)*mp.gamma(1+tv)*mp.gamma(1+uv))
for eps in [mp.mpf('1e-3'), mp.mpf('1e-4')]:
    sv, tv = eps*mp.mpf('0.7'), -eps*mp.mpf('0.3')
    uv = -sv-tv
    print("VS - 1/stu at |s|~%s:" % mp.nstr(eps, 2), mp.nstr(VS(sv, tv) - 1/(sv*tv*uv), 12), " 2*zeta3 =", mp.nstr(2*mp.zeta(3), 12))

print("=== (b) level truncation, 2406.02665 eq.(5) ===")
def AV(sv, tv):
    return mp.gamma(-sv)*mp.gamma(-tv)/mp.gamma(-sv-tv)
for k in [1, 2, 3, 5]:
    sv = mp.mpf('0.37')
    tv = -k + mp.mpf('1e-20')
    lhs = AV(sv, tv)
    rhs = sum(mp.rf(1-k, n)/mp.factorial(n)/(n - sv) for n in range(k))
    print("k=%d" % k, mp.nstr(lhs, 15), mp.nstr(rhs, 15))

print("=== (c) LP for g3~ lower bound, 2011.02957 (3.10)-(3.11), d=4 ===")
d = 4
# functionals (per heavy state, M=1, x = 1/m^2 in (0,1]):
#   g2: x^2 ; g3: (3 - 4/(d-2) J2) x^3 ; n4: J2 (2 J2 - (5d-4)) x^4 ;  J2 = J(J+d-3)
# find min alpha s.t. alpha*g2 + g3 + beta*n4 >= 0 for all states -> g3~ >= -alpha
xs = np.linspace(1e-4, 1.0, 4000)
Js = np.arange(0, 400, 2)
A_ub, b_ub = [], []
for J in Js:
    J2 = J*(J+d-3)
    for x in xs:
        g2 = x**2; g3 = (3 - 4/(d-2)*J2)*x**3; n4 = J2*(2*J2 - (5*d-4))*x**4
        # alpha*g2 + beta*n4 >= -g3   ->  -g2*alpha - n4*beta <= g3
        A_ub.append([-g2, -n4]); b_ub.append(g3)
res = linprog(c=[1, 0], A_ub=np.array(A_ub), b_ub=np.array(b_ub), bounds=[(None, None), (None, None)], method='highs')
print("status", res.status, " alpha* =", res.x[0], " beta* =", res.x[1], " (paper: 10.61249, 0.0671875)")
