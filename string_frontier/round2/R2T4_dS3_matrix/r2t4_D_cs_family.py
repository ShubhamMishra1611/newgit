"""R2-T4 script D: is CEM (85) a member of the Chern-Simons family of all-loop S^3 answers?

CS family (ADLS 2009.12464 (H.24) with SO(4) normalization; CCFZ 2302.12281 (2.85)): renormalized levels
K = L +- i s, s = S0/2pi (L = 2 for ADLS (8.18) [bare l=0]; L = 0 for CCFZ (2.84) [delta-hat=0, Re c = 13]):
   |Z_L| = e^{S0} * (4/|K|) |sin(pi/K)|^2.
Compute the invariant coefficients d_2 (=c_1) and d_4 (=c_2) as functions of L and ask whether some real L
reproduces CEM (85)'s d_2 = 4pi^2(pi^2/3 - 2) and d_4 = 16 pi^4 - 8 pi^8/45 simultaneously.
"""
import sympy as sp
u, L = sp.symbols('u L', positive=True)
PI = sp.pi
s = 1/(2*PI*u)
K = L + sp.I*s; Kb = L - sp.I*s
logZ_minus_S0 = sp.log(4) - sp.Rational(1, 2)*sp.log(L**2 + s**2) + sp.log(sp.sin(PI/K)*sp.sin(PI/Kb))
reg = sp.expand_log(logZ_minus_S0 - 3*sp.log(u), force=True)   # remove -3 log S0 = +3 log u
ser = sp.series(sp.simplify(reg), u, 0, 6).removeO()
ser = sp.expand(ser)
d = {k: sp.simplify(ser.coeff(u, k)) for k in range(0, 6)}
d[0] = sp.simplify(ser.subs(u, 0))
for k in range(0, 6):
    print(f"d_{k}(L) =", sp.factor(d[k]))
d2_CEM = 4*PI**2*(PI**2/3 - 2)
d4_CEM = 16*PI**4 - 8*PI**8/45
Lsq = sp.solve(sp.Eq(d[2], d2_CEM), L)
print("\nL solving d_2(L) = d_2(CEM):", Lsq, " L^2 =", [sp.nsimplify(x**2) for x in Lsq])
for x in Lsq:
    print("  at that L: d_4(L) =", sp.N(d[4].subs(L, x), 12), " vs d_4(CEM) =", sp.N(d4_CEM, 12))
print("\nValues: d_2(L=0) =", sp.N(d[2].subs(L, 0), 10), " d_2(L=2) =", sp.N(d[2].subs(L, 2), 10), " d_2(CEM) =", sp.N(d2_CEM, 10))

# ---------------------------------------------------------------------------------------------------------
# Part 2: Blommaert-Tietto-Verlinde 2606.26241, eq. (5.32): rho_(0,1)(alpha+-=1) = 16 sinh^2(hbar/4) sinh^2(4pi^2/hbar)
# "= S_11^2 = e^{S_GH}", said to be "completely consistent" with CEM (5.7). BTV dictionary: c = 13 +- 6i(b^2 - b^-2),
# hbar = 4 pi b^2 (b real)  =>  BTV b^2 == CEM beta (CEM b^2 = i beta).  Compare log coefficients.
# Also BTV's description of CEM's count: integrate rho_disk(a+) rho_disk(a-) over [0,2pi]^2, rho_disk = sin(a/2) sinh(2 pi a/hbar) (5.31).
print("\n" + "="*90)
bt = sp.symbols('beta', positive=True)
absCEM = 16*bt*sp.sinh(PI*bt)**2*sp.sinh(PI/bt)**2/(1+bt**2)**2
BTV532 = 16*sp.sinh(PI*bt)**2*sp.sinh(PI/bt)**2            # hbar = 4 pi beta
print("BTV (5.32) / |CEM (5.7)| =", sp.simplify(BTV532/absCEM))
al = sp.symbols('alpha', positive=True)
I1 = sp.simplify(sp.integrate(sp.sin(al/2)*sp.sinh(2*PI*al/(4*PI*bt)), (al, 0, 2*PI)))
print("int_0^{2pi} rho_disk =", I1)
# log coefficients a in log X = S0 - a log S0 + ..., S0 = 2pi/beta (beta -> 0):
import mpmath as mp
mp.mp.dps = 40
def a_num(f):
    # numerical: a = d/d(log beta) [log f - 2pi/beta] at small beta
    g = lambda lb: mp.log(f(mp.e**lb)) - 2*mp.pi/mp.e**lb
    return mp.diff(g, mp.log(mp.mpf('1e-3')))
fCEM = lambda b_: 16*b_*mp.sinh(mp.pi*b_)**2*mp.sinh(mp.pi/b_)**2/(1+b_**2)**2
fBTV = lambda b_: 16*mp.sinh(mp.pi*b_)**2*mp.sinh(mp.pi/b_)**2
fINT = lambda b_: (2*b_**2*mp.sinh(mp.pi/b_)/(1+b_**2))**2
for name, f in [("|CEM (5.7)|", fCEM), ("BTV (5.32)", fBTV), ("[int rho_disk]^2 (BTV normalization)", fINT)]:
    print(f"  {name:40s}: log coefficient a (log X = S0 - a log S0 + ...) = {mp.nstr(a_num(f), 8)}")
print("=> the same sine/sinh skeleton gives a = 3, 2 or 4 depending on b-dependent normalization prefactors;\n"
      "   the ADLS value a = 3 is carried entirely by the h-basis Jacobian / e^{S0}(b) normalization, not by the skeleton.")
