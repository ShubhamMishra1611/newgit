"""R1: log-coefficient alpha_4 from 2009.12464 eq. (1.4),(1.9),(1.11),(1.13):
integrand = (1/2t)(1+q)/(1-q) (chi_bulk - chi_edge), q = e^{-t}; alpha_4 = coefficient of 1/t.
Characters as printed in (1.9)/(1.11) for d=3:
  s=1: bulk (6q^2-2q^3)/(1-q)^3, edge 2q/(1-q);   s=2: bulk (10q^3-6q^4)/(1-q)^3, edge (10q^2-2q^3)/(1-q)
Cross-check with the independent heat-kernel / trace-anomaly literature values:
  Maxwell on S^4: -31/45 = -4 a_Maxwell with a = 62/360 ; conformal scalar check: chi=(q+q^2)/(1-q)^3 (nu=i/2) -> -1/90.
"""
import sympy as sp
t = sp.symbols('t', positive=True)
q = sp.exp(-t)
def alpha(chi):
    integrand = (1+q)/(1-q)*chi/(2*t)
    return sp.nsimplify(sp.series(integrand, t, 0, 1).removeO().coeff(t, -1))
for name, cb, ce in [("photon", (6*q**2-2*q**3)/(1-q)**3, 2*q/(1-q)),
                     ("graviton", (10*q**3-6*q**4)/(1-q)**3, (10*q**2-2*q**3)/(1-q))]:
    ab, ae = alpha(cb), alpha(-ce)
    print(f"{name}: bulk {ab}, edge {ae}, total {ab+ae}")
# conformally coupled scalar in dS4: Delta = 1,2 -> chi = (q + q^2)/(1-q)^3 ; expected a-anomaly: -4a = -4/360 = -1/90
print("conformal scalar:", alpha((q+q**2)/(1-q)**3), " expected -1/90")
