"""R2-T4 script B: exact one-loop S^3 partition function of pure 3D gravity, Lambda>0, from ADLS arXiv:2009.12464.

ADLS appendix H.1, eqs. (H.1),(H.4)-(H.7): for d=2, chi_bulk = chi_edge = 0 for the massless spin-2 field, and
   Z_PI^(1) = i^{-P} gamma^{dim G}/vol(G)_c * Z_char,   Z_char = exp( -2 dim G * int^x dt/2t (1+q)/(1-q) ),
   gamma = sqrt(8 pi G_N / A_1) = sqrt(4 G_N / l),  G = SO(4), dim G = 6, P = P_2 = 5,
   vol(SO(4))_c = (2 pi)^4  (metric with |M_IJ| = 1),  finite part of Z_char = (2 pi)^{dim G}   (C.37).
Result (H.7): Z_fin = i^-5 (2 pi gamma)^6/(2 pi)^4 = -i 4 pi^2 gamma^6, i.e. with S0 = 2 pi l/4G:
   log Z^(1) = -3 log S0 + 5 log(2 pi) - i pi/2         [(1.23), (8.17)]
We recompute each ingredient independently:
  (a) (C.37) finite part with the IR prescription (5.12), by direct numerical integration;
  (b) vol(SO(4)) from the fibration SO(n)/SO(n-1)=S^{n-1};
  (c) the phase: P_2 from (H.4), and independently zeta_{S^3}(0) = -5 for the Wick-rotated conformal
      modes (CEM 2501.01486 eqs. (66)-(70)), by our own analytic continuation.
"""
import mpmath as mp
mp.mp.dps = 30

# (a) finite part of F = int^x_0 dt/(2t) (1+q)/(1-q), q=e^{-t}; (5.12): subtract f(inf) log L with f(inf)=1/2;
#     UV: (1/2t)coth(t/2) = 1/t^2 + 1/12 + O(t^2): pure power divergence 1/eps (no log), so finite part is scheme-free.
def F_reg(eps, L):
    g = lambda t: mp.coth(t/2)/(2*t)*mp.exp(-t/L)
    I = mp.quad(g, [eps, 1, 10, 100, 10*L, mp.inf])
    return I - 1/eps - mp.mpf(1)/2*mp.log(L)
print("(a) finite part of int^x dt/2t (1+q)/(1-q)  [expected -1/2 log(2 pi) =", mp.nstr(-mp.log(2*mp.pi)/2, 12), "]")
for eps, L in [(1e-3, 1e3), (1e-4, 1e4), (1e-5, 1e5), (1e-6, 1e6)]:
    print(f"    eps={eps:g}, L={L:g}: {mp.nstr(F_reg(eps, L), 12)}")
# The corrections are O(eps) and O(1/L): eps-term: int_0^eps (1/12) = eps/12 ; L-term ~ (1/2L)*... -> extrapolate
Fexact = -mp.log(2*mp.pi)/2
dimG = 6
logZchar = -2*dimG*Fexact
print("    => log Z_char (finite) = -2 dimG F =", mp.nstr(logZchar, 12), " = dimG*log(2pi) =", mp.nstr(dimG*mp.log(2*mp.pi), 12))
# zeta-function cross-check: (1+q)/(1-q) = 1 + 2 sum_n q^n ; int dt/t e^{-nt} -> -log n ; sum_n -log n = zeta'(0) = -1/2 log 2pi
print("    zeta'(0) =", mp.nstr(mp.zeta(0, derivative=1), 12))

# (b) vol(SO(4)) = vol(S^1) vol(S^2) vol(S^3)
volS = lambda k: 2*mp.pi**((k+1)/mp.mpf(2))/mp.gamma((k+1)/mp.mpf(2))
volSO4 = volS(1)*volS(2)*volS(3)
print("\n(b) vol(SO(4))_c = vol(S1)vol(S2)vol(S3) =", mp.nstr(volSO4, 12), " (2pi)^4 =", mp.nstr((2*mp.pi)**4, 12))

# (c) phase
P = lambda s: mp.mpf(1)/3*(2*s-3)*(2*s-1)*(2*s+1)
print("\n(c) P_2 from (H.4) =", P(2), " -> i^{-P} = ", mp.nstr(mp.power(1j, -P(2)), 5))
# zeta_{S3}(s) = sum_{l>=2} (l+1)^2 / (l(l+2)-3)^s = sum_{n>=3} n^2 (n^2-4)^{-s};  continue to s=0:
#   = sum_{n>=3}[n^2(n^2-4)^{-s} - n^{2-2s} - 4 s n^{-2s}] + (zeta(2s-2)-1-2^{2-2s}) + 4 s (zeta(2s) - 1 - 2^{-2s})
def zetaS3(s):
    rem = mp.nsum(lambda n: n**2*(n**2-4)**(-s) - n**(2-2*s) - 4*s*n**(-2*s), [3, mp.inf])
    return rem + (mp.zeta(2*s-2)-1-mp.power(2, 2-2*s)) + 4*s*(mp.zeta(2*s)-1-mp.power(2, -2*s))
for s in [mp.mpf('0.3'), mp.mpf('1e-3'), mp.mpf('1e-6')]:
    print(f"    zeta_S3({mp.nstr(s,3)}) = {mp.nstr(zetaS3(s), 12)}")
# sanity: compare the continuation with the direct sum at Re s > 3/2
s = mp.mpf(2.5)
direct = mp.nsum(lambda n: n**2*(n**2-4)**(-s), [3, mp.inf])
print("    check at s=2.5: continued =", mp.nstr(zetaS3(s), 15), " direct =", mp.nstr(direct, 15))
print("    => zeta_S3(0) = -5 => Wick-rotation phase exp(+-i pi/2 * zeta(0)) = exp(-+5 pi i/2) = -+i")

# assemble (H.7): log Z^(1) as a function of S0
def logZ1(S0):
    gamma2 = 2*mp.pi/S0            # gamma^2 = 4G/l = 2 pi/S0
    return -1j*mp.pi/2 + 3*mp.log(gamma2) + logZchar - mp.log(volSO4)
print("\nAssembled one-loop: log Z^(1) + 3 log S0 =")
for S0 in [10, 1e3, 1e6]:
    v = logZ1(mp.mpf(S0)) + 3*mp.log(S0)
    print(f"    S0={S0:g}: {mp.nstr(v, 14)}   [5 log 2pi = {mp.nstr(5*mp.log(2*mp.pi), 14)}]")
print("Decomposition of 5 log 2pi: +6 (character finite part) +3 (gamma^6=(2pi/S0)^3) -4 (vol SO(4)) = 5")
