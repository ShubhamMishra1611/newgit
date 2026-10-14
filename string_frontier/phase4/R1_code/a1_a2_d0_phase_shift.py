"""
R1 independent derivation of the one-loop D0-D0 effective action coefficients (A1) and the
probe/Aichelburg-Sexl three-loop target (A2).

One loop (derived here from the fluctuation spectrum, not copied):
  Euclidean background X1 = v tau, X2 = b.  Off-diagonal fluctuation operator H = -d^2 + b^2 + v^2 tau^2,
  spectrum b^2 + v(2n+1).  Background-gauge counting (complex off-diagonal modes):
    10 bosons (A0, X1..X9) - 2 ghosts: 6 with H, 1 with H+2v, 1 with H-2v
    16 fermion components -> det^{+4}(H+v) det^{+4}(H-v)
  ln Z = int ds/s e^{-s b^2} [6 + 2cosh 2sv - 8 cosh sv] / (2 sinh sv)
Matching to a Minkowski Lagrangian L = v^2/2 + sum_k c_k v^{2k}/r^{4k-1}:
  Wick rotation v_M = i v_E gives ln Z = int dtau sum_k (-1)^k c_k v_E^{2k}/r^{4k-1};
  int dtau (b^2+v^2tau^2)^{-q/2} = sqrt(pi) Gamma((q-1)/2)/Gamma(q/2) b^{1-q}/v ;
  int ds/s e^{-sb^2} s^m = Gamma(m) b^{-2m}.
Units: action (1/2)Tr Xdot^2 + (1/4)Tr[X,X]^2 (BBPT eq. (2) with M=R=1, g = 2R ... see report).
"""
import sympy as sp

x = sp.symbols('x')
f = 6 + 2*sp.cosh(2*x) - 8*sp.cosh(x)
ser = sp.series(f/(2*sp.sinh(x)), x, 0, 22).removeO()
print("f(x)/(2 sinh x) =", sp.expand(ser))
print()
print("One-loop Minkowski coefficients c_k of v^{2k}/r^{4k-1}  (N1=N2=1, units of BBPT eq.(2) with M=R=1):")
for m in range(3, 21, 2):
    a_m = ser.coeff(x, m)
    k = (m + 1)//2
    q = 2*m + 1
    # ln Z term: a_m v^m Gamma(m) b^{-2m}  ==  (-1)^k c_k v^{2k-1} b^{1-q} sqrt(pi) Gamma(m)/Gamma(q/2)
    c_k = sp.nsimplify(a_m*sp.gamma(sp.Rational(q, 2))/sp.sqrt(sp.pi)*(-1)**k)
    print(f"  v^{2*k}/r^{q}: c = {c_k}")

print()
print("BBPT eq.(12)-(13): L' = -p_- xdot^-,  xdot^- = (sqrt(1 - h v^2) - 1)/h,  h = 15 N1/(2 R^2 M^9 r^7), p_- = N2/R")
h, v, N1, N2, R, M, r = sp.symbols('h v N1 N2 R M r', positive=True)
xdot = (sp.sqrt(1 - h*v**2) - 1)/h
Lp = sp.series(-(N2/R)*xdot, h, 0, 5).removeO()
Lp = sp.expand(Lp.subs(h, sp.Rational(15, 2)*N1/(R**2*M**9*r**7)))
for term in sp.Add.make_args(Lp):
    print("  ", term)
