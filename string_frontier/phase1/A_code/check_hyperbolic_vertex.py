"""Closed-string tachyon cubic coupling for the hyperbolic 3-vertex (Firat, arXiv:2102.03936, eqs. (1.5)-(1.8)).
kappa^2 V = -t^2 + (1/3) r^6 t^3  (alpha'=2), r = mapping radius of the equal-length (L=2 pi lambda) hyperbolic vertex.
Firat quotes C(lambda*) ~ 1.62187e8 at lambda* = arcsinh(1)/pi (eq. 1.8), and the minimal-area value ~1.602,
which is 6561/4096 (Yang-Zwiebach hep-th/0506077 eq. (1.1)).
Our reading (the extracted text was ambiguous on branch/sign; this is the reading that reproduces BOTH numbers):
  v(lam) = (1/2i) [ 2 logGamma(-i lam) - 2 logGamma(i lam) + log gamma((1-i lam)/2) - log gamma((1+3 i lam)/2) ],
  gamma(x) = Gamma(1-x)/Gamma(x), logs on the continuous (loggamma) branch;  ln r = (v - pi/2)/lam  (scale-factor integer l=-1).
"""
import mpmath as mp
mp.mp.dps = 30
I = mp.mpc(0, 1)
def lg(x): return mp.loggamma(1-x) - mp.loggamma(x)
def v(l):
    return mp.re((2*mp.loggamma(-I*l) - 2*mp.loggamma(I*l) + lg((1-I*l)/2) - lg((1+3*I*l)/2))/(2*I))
def C(l):
    return mp.exp(6*(v(l) - mp.pi/2)/l)/3
lam_star = mp.asinh(1)/mp.pi
print("lambda* = arcsinh(1)/pi =", mp.nstr(lam_star, 10))
print("C(lambda*) =", mp.nstr(C(lam_star), 8), "   [Firat eq.(1.8): 1.62187e8]")
for l in [0.5, 1, 2, 5, 10, 30, 100, 1000]:
    print(f"C({l}) =", mp.nstr(C(l), 10))
print("minimal-area 6561/4096 =", mp.nstr(mp.mpf(6561)/4096, 10), " (3*sqrt(3)/4)^6/3 =", mp.nstr((3*mp.sqrt(3)/4)**6/3, 10))
assert abs(C(lam_star)/mp.mpf('1.62187e8') - 1) < 1e-5
assert abs(C(1000) - mp.mpf(6561)/4096) < 2e-3
print("CHECK PASSED")
