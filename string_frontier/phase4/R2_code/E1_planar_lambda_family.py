"""
R2 independent check of claim E1.

Planar lambda-family of Cheung-Hillman-Remmen 2408.03362 (v2 numbering: eqs. 27-29, 34;
the phase-1 doc calls it "eq. 18" which is v1 numbering -- checked separately).

  A(s,t) = sum_{n>=0} Rhat(n,t)/(mu(n)-s),  mu(n)=(n+lam-1)/lam,
  Rhat(0,t)=rho/t,  Rhat(n,t) = ((1+lam)/2 + lam t)_{n-1} / ((1+3lam)/2)_{n-1}   (n>=1)
  closed form / satellite sum:
  A = -rho/(st) + sum_k lam(3lam-1) Gamma(k+lam-lam s)Gamma(k+lam-lam t)
                     / ((3lam-1+2k) k! Gamma(k+2lam-lam s-lam t))

Map to Elvang-Herderschee-Morales 2601.11705 (v1 eq. 1.3 / v2 eq. 1.3):
  F(s,u) = -g^2/(su) + sum a_{k,q} s^{k-q} u^q ,  planar adjacent channels (s,u)_EHM <-> (s,t)_CHR,
  g^2 <-> rho.  At lam=rho=1 A = -Gamma(-s)Gamma(-t)/Gamma(1-s-t) = EHM F with g=1 (checked below).

Nonlinear relations (2601.11705 v1 eq. 1.6 = v2 eq. 1.11; v2 requires extra 6-scalar parity condition):
  g^2 a20 = 2/5 a00^2,  g^2 a21 = 1/10 a00^2,  g^2 a31 = 2 g^2 a30 - a10 a00.
  rho-independent consequences:  a21/a20 = 1/4 ;   (2/5) a00 (2 a30 - a31) = a10 a20.

Two independent evaluations:
  (M1) dual-resonant sum over levels n (mpmath nsum with acceleration),
  (M2) bivariate Taylor expansion of the satellite sum (each Gamma ratio expanded with polygammas).
"""
import mpmath as mp
mp.mp.dps = 30

# ---------- M1: level sum ----------
def coeffs_levelsum(lam):
    lam = mp.mpf(lam)
    a0 = (1+lam)/2; b0 = (1+3*lam)/2
    mu = lambda n: (n+lam-1)/lam
    R0 = lambda n: mp.rf(a0, n-1)/mp.rf(b0, n-1)
    # t-derivatives of log Rhat: d/dt log (a0+lam t)_{n-1} = lam * sum_{j<n-1} 1/(a0+j+lam t)
    d1 = lambda n: lam*(mp.digamma(a0+n-1)-mp.digamma(a0))            # first log-derivative at t=0
    d2 = lambda n: lam**2*(mp.psi(1,a0+n-1)-mp.psi(1,a0))              # second log-derivative = -lam^2 sum 1/(a0+j)^2
    R1 = lambda n: R0(n)*d1(n)                                          # dR/dt
    R2h = lambda n: R0(n)*(d1(n)**2 + d2(n))/2                          # (1/2) d^2R/dt^2
    S = lambda f: mp.nsum(f, [1, mp.inf])
    # coefficient of s^p t^q  =  sum_n [t^q Rhat(n,t)] / mu^{p+1}
    c = {}
    c[(0,0)] = S(lambda n: R0(n)/mu(n))
    c[(1,0)] = S(lambda n: R0(n)/mu(n)**2)
    c[(0,1)] = S(lambda n: R1(n)/mu(n))
    c[(2,0)] = S(lambda n: R0(n)/mu(n)**3)
    c[(1,1)] = S(lambda n: R1(n)/mu(n)**2)
    c[(0,2)] = S(lambda n: R2h(n)/mu(n))
    c[(3,0)] = S(lambda n: R0(n)/mu(n)**4)
    c[(2,1)] = S(lambda n: R1(n)/mu(n)**3)
    c[(1,2)] = S(lambda n: R2h(n)/mu(n)**2)
    return c

# ---------- M2: satellite-sum Taylor expansion ----------
DEG = 4
def pmul(A, B):
    C = {}
    for (i,j),x in A.items():
        for (k,l),y in B.items():
            if i+j+k+l <= DEG:
                C[(i+k,j+l)] = C.get((i+k,j+l),0)+x*y
    return C
def pexp(A):
    # exp of polynomial with zero constant term
    res = {(0,0):mp.mpf(1)}; term = {(0,0):mp.mpf(1)}
    for m in range(1, DEG+1):
        term = pmul(term, A); term = {k:v/m for k,v in term.items()}
        for k,v in term.items(): res[k] = res.get(k,0)+v
    return res
def loggamma_series(a, cs, ct):
    # log Gamma(a + cs*s + ct*t) - log Gamma(a), to total degree DEG
    out = {}
    for m in range(1, DEG+1):
        der = mp.polygamma(m-1, a)/mp.factorial(m)
        # (cs s + ct t)^m
        for i in range(m+1):
            coef = mp.binomial(m,i)*cs**i*ct**(m-i)
            if coef != 0:
                out[(i,m-i)] = out.get((i,m-i),0)+der*coef
    return out
def satellite_term(lam, k):
    pre = lam*(3*lam-1)*mp.gamma(k+lam)**2/((3*lam-1+2*k)*mp.factorial(k)*mp.gamma(k+2*lam))
    L = {}
    for part, sgn in ((loggamma_series(k+lam, -lam, 0), 1),
                      (loggamma_series(k+lam, 0, -lam), 1),
                      (loggamma_series(k+2*lam, -lam, -lam), -1)):
        for kk,v in part.items(): L[kk] = L.get(kk,0)+sgn*v
    E = pexp(L)
    return {kk: pre*v for kk,v in E.items()}
def coeffs_satellite(lam):
    lam = mp.mpf(lam)
    keys = [(i,j) for i in range(DEG+1) for j in range(DEG+1) if i+j<=3]
    c = {}
    for key in keys:
        c[key] = mp.nsum(lambda k: satellite_term(lam, int(k))[key], [0, mp.inf])
    return c

def coeffs_hyp(lam):
    # closed forms: sum_{m>=0} (a0)_m/(b0)_m * lam^{p+1}/(m+lam)^{p+1} = (p+2)F(p+1)(1,a0,lam..;b0,1+lam..;1)
    lam = mp.mpf(lam); a0=(1+lam)/2; b0=(1+3*lam)/2
    out = {}
    for p in range(4):
        out[(p,0)] = mp.hyper([1,a0]+[lam]*(p+1), [b0]+[1+lam]*(p+1), 1)
    return out

def report(lam, c, tag):
    a00=c[(0,0)]; a10=c[(1,0)]; a20=c[(2,0)]; a21=c[(1,1)]; a30=c[(3,0)]; a31=c[(2,1)]
    print(f"  [{tag}] a00={mp.nstr(a00,12)} a10={mp.nstr(a10,12)} (a11={mp.nstr(c[(0,1)],12)})"
          f" a20={mp.nstr(a20,12)} a22={mp.nstr(c[(0,2)],12)} a21={mp.nstr(a21,12)}")
    print(f"        a30={mp.nstr(a30,12)} a31={mp.nstr(a31,12)} a32={mp.nstr(c[(1,2)],12)}")
    print(f"        a21/a20 = {mp.nstr(a21/a20,10)}   (SUSY+parity requires 0.25)")
    print(f"        rho needed by g^2 a20=2/5 a00^2 : {mp.nstr(mp.mpf(2)/5*a00**2/a20,10)};"
          f"  by g^2 a21=a00^2/10 : {mp.nstr(a00**2/(10*a21),10)}")
    lhs = mp.mpf(2)/5*a00*(2*a30-a31); rhs = a10*a20
    print(f"        third relation (rho-free form) (2/5)a00(2a30-a31) / (a10 a20) = {mp.nstr(lhs/rhs,10)} (needs 1)")
    return a21/a20

if __name__ == "__main__":
    # sanity: lam=1 satellite sum reproduces superstring Veneziano
    s0, t0 = mp.mpf('0.13'), mp.mpf('-0.27')
    sat = -1/(s0*t0) + mp.nsum(lambda k: 2*mp.gamma(k+1-s0)*mp.gamma(k+1-t0)/((2+2*k)*mp.factorial(k)*mp.gamma(k+2-s0-t0)), [0, mp.inf])
    ven = -mp.gamma(-s0)*mp.gamma(-t0)/mp.gamma(1-s0-t0)
    print("lam=1 satellite sum vs -G(-s)G(-t)/G(1-s-t):", mp.nstr(sat,15), mp.nstr(ven,15))
    print("zeta4/4 =", mp.nstr(mp.zeta(4)/4,12), " zeta4 =", mp.nstr(mp.zeta(4),12))
    for lam in ['0.25','0.5','0.9','1','1.1','2','4']:
        print(f"lambda = {lam}")
        r1 = report(lam, coeffs_levelsum(lam), "level sum ")
        r2 = report(lam, coeffs_satellite(lam), "satellite ")
        h = coeffs_hyp(lam)
        print("        hypergeometric closed forms a00,a10,a20,a30 =", [mp.nstr(h[(p,0)],12) for p in range(4)])
        print(f"        |diff ratio| = {mp.nstr(abs(r1-r2),3)}")
