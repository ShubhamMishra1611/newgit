"""
R2 checks of E4 (VS low-energy expansion) and E5 (closed-string lambda-family (g3~, g4~) table).
E4: series of M_VS = -G(-s)G(-t)G(-u)/[G(1+s)G(1+t)G(1+u)], u=-s-t, via sympy (exact).
E5: closed-string family of 2408.03362 (v1 eq. 8, satellite sum), rho=0:
    M = sum_k lam(3lam-1)^2 Gamma(3lam-1+k)/((3lam-1+2k) k!) * M_k,
    M_k = G(k+lam-lam s)G(k+lam-lam t)G(k+lam-lam u)/[G(k+2lam+lam s)G(k+2lam+lam t)G(k+2lam+lam u)]
  low energy  M = g0 + g2 (s^2+t^2+u^2) + g3 stu + g4 (s^2+t^2+u^2)^2 + ...  (v1 eq. 9; same form as CHVD 2011.02957 eq. 2.3)
  At t=0: M = g0 + 2 g2 s^2 + 4 g4 s^4;  coefficient of s^2 t (u=-s-t) = -g3.
  Units: M_gap = mu(1) = 1.  g3~ = g3/g2, g4~ = g4/g2.
  Cross-check g2, g4 against v1 eq.10 hypergeometrics.
"""
import sympy as sp
import mpmath as mp
mp.mp.dps = 25

# ---------------- E4 ----------------
s, t = sp.symbols('s t')
u = -s - t
eps = sp.symbols('eps')
# log of Gamma(1-x)/Gamma(1+x) = 2 gamma x + sum_{k odd>=3} 2 zeta(k) x^k / k  (exact identity)
def lograt(x, K=9):
    return 2*sp.EulerGamma*x + sum(2*sp.zeta(k)*x**k/k for k in range(3, K+1, 2))
expo = sp.expand(lograt(eps*s) + lograt(eps*t) + lograt(eps*u))
ser = sp.series(sp.exp(expo), eps, 0, 10).removeO()
# M_VS = 1/(stu) * Gamma(1-s)G(1-t)G(1-u)/[G(1+s)..]  (since -G(-x) = G(1-x)/x, three factors -> -1/(-stu)... check sign numerically below)
stu = s*t*u
poly = sp.expand(ser)
low = {}
for n in [3, 5, 6, 7]:
    low[n] = sp.factor(sp.expand(poly.coeff(eps, n)))
print("E4: [eps^3] =", low[3], "  -> /stu =", sp.simplify(low[3]/stu))
print("    [eps^5] / stu =", sp.factor(sp.simplify(low[5]/stu)), "  (expect zeta5*(s^2+t^2+u^2))", sp.simplify(low[5]/stu - sp.zeta(5)*(s**2+t**2+u**2)))
print("    [eps^6] / stu =", sp.factor(sp.simplify(low[6]/stu)), "  (expect 2 zeta3^2 stu)")
print("    [eps^7] / stu - zeta7/2 (s^2+t^2+u^2)^2 =", sp.simplify(low[7]/stu - sp.zeta(7)/2*(s**2+t**2+u**2)**2))
# numeric sign check
sv, tv = mp.mpf('0.07'), mp.mpf('-0.11'); uv = -sv-tv
MVS = -mp.gamma(-sv)*mp.gamma(-tv)*mp.gamma(-uv)/(mp.gamma(1+sv)*mp.gamma(1+tv)*mp.gamma(1+uv))
approx = 1/(sv*tv*uv) + 2*mp.zeta(3) + mp.zeta(5)*(sv**2+tv**2+uv**2) + 2*mp.zeta(3)**2*sv*tv*uv
print("    numeric: M_VS =", mp.nstr(MVS, 15), " truncated series =", mp.nstr(approx, 15))

# ---------------- E5 ----------------
DEG = 4
def pmul(A, B):
    C = {}
    for (i,j),x in A.items():
        for (k,l),y in B.items():
            if i+j+k+l <= DEG:
                C[(i+k,j+l)] = C.get((i+k,j+l),0)+x*y
    return C
def pexp(A):
    res = {(0,0):mp.mpf(1)}; term = {(0,0):mp.mpf(1)}
    for m in range(1, DEG+1):
        term = pmul(term, A); term = {k:v/m for k,v in term.items()}
        for k,v in term.items(): res[k] = res.get(k,0)+v
    return res
def lg_series(a, cs, ct):
    out = {}
    for m in range(1, DEG+1):
        der = mp.polygamma(m-1, a)/mp.factorial(m)
        for i in range(m+1):
            coef = mp.binomial(m,i)*cs**i*ct**(m-i)
            if coef != 0:
                out[(i,m-i)] = out.get((i,m-i),0)+der*coef
    return out
def Mk_series(lam, k):
    pre = lam*(3*lam-1)**2*mp.gamma(3*lam-1+k)/((3*lam-1+2*k)*mp.factorial(k)) \
          * mp.gamma(k+lam)**3/mp.gamma(k+2*lam)**3
    L = {}
    parts = [(lg_series(k+lam, -lam, 0), 1), (lg_series(k+lam, 0, -lam), 1), (lg_series(k+lam, lam, lam), 1),   # u=-s-t
             (lg_series(k+2*lam, lam, 0), -1), (lg_series(k+2*lam, 0, lam), -1), (lg_series(k+2*lam, -lam, -lam), -1)]
    for part, sg in parts:
        for kk,v in part.items(): L[kk] = L.get(kk,0)+sg*v
    E = pexp(L)
    return {kk: pre*v for kk,v in E.items()}
def family_g(lam):
    lam = mp.mpf(lam)
    c = {key: mp.nsum(lambda k: Mk_series(lam, int(k))[key], [0, mp.inf]) for key in [(0,0),(2,0),(2,1),(4,0)]}
    g0 = c[(0,0)]; g2 = c[(2,0)]/2; g3 = -c[(2,1)]; g4 = c[(4,0)]/4
    return g0, g2, g3, g4
def hyp_g2_g4(lam):
    lam = mp.mpf(lam); a=(1+lam)/2; b=(1+3*lam)/2
    g2 = mp.hyper([1,lam,lam,lam,a,a],[1+lam,1+lam,1+lam,b,b],1)
    g4 = mp.hyper([1,lam,lam,lam,lam,lam,a,a],[1+lam]*5+[b,b],1)/2
    g0 = 2*mp.hyper([1,lam,a,a],[1+lam,b,b],1)
    return g0, g2, g4

if __name__ == "__main__":
    print("\nE5: closed-string lambda family, rho=0, units mu(1)=1")
    for lam in ['0.5','1','2','5']:
        g0,g2,g3,g4 = family_g(lam)
        h0,h2,h4 = hyp_g2_g4(lam)
        print(f"  lam={lam}: g0={mp.nstr(g0,8)} (hyp {mp.nstr(h0,8)}), g2={mp.nstr(g2,8)} (hyp {mp.nstr(h2,8)}),"
              f" g4={mp.nstr(g4,8)} (hyp {mp.nstr(h4,8)});  g3~={mp.nstr(g3/g2,7)}  g4~={mp.nstr(g4/g2,7)}")
    print("  VS check: 2 zeta3^2/zeta5 =", mp.nstr(2*mp.zeta(3)**2/mp.zeta(5),7), " zeta7/(2 zeta5) =", mp.nstr(mp.zeta(7)/(2*mp.zeta(5)),7))
