# Reproduce ADLS 2009.12464 eq.(1.12)-(1.13): log-coefficient alpha_4 from t->0 expansion of
# integrand (1/2t)*(1+q)/(1-q)*chi_tot, q=e^{-t} (ell=1), chi_tot = chi_bulk - chi_edge, characters from (1.9),(1.11).
import sympy as sp
t=sp.symbols('t',positive=True)
q=sp.exp(-t)
def alpha(chi_b, chi_e):
    f=lambda chi: sp.series(sp.Rational(1,2)/t*(1+q)/(1-q)*chi, t, 0, 1).removeO().coeff(t,-1)
    return sp.nsimplify(f(chi_b)), sp.nsimplify(-f(chi_e))
# d=3 (4D): photon (3,1), graviton (3,2)
ph_b=(6*q**2-2*q**3)/(1-q)**3; ph_e=2*q/(1-q)
gr_b=(10*q**3-6*q**4)/(1-q)**3; gr_e=(10*q**2-2*q**3)/(1-q)
for name,(b,e) in {'photon':(ph_b,ph_e),'graviton':(gr_b,gr_e)}.items():
    ab,ae=alpha(b,e)
    print(name,'bulk',ab,'edge',ae,'total',ab+ae)
# conformally coupled / minimally coupled scalar d=3: chi = (q^{3/2+i nu}+q^{3/2-i nu})/(1-q)^3, m^2=2 -> nu=i/2
nu=sp.I/2
sc=(q**(sp.Rational(3,2)+sp.I*nu)+q**(sp.Rational(3,2)-sp.I*nu))/(1-q)**3
print('conformal scalar (m^2 l^2=2):', sp.nsimplify(sp.series(sp.Rational(1,2)/t*(1+q)/(1-q)*sc,t,0,1).removeO().coeff(t,-1)))
