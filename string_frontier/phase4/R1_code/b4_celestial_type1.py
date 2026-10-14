"""
R1 independent check of Stieberger-Taylor 1806.05688, eqs. (3.15), (3.30), (3.31), (3.32).

Definitions read from the fetched paper text:
  F_I(s,u) = Gamma(1-s)Gamma(1-u)/Gamma(1-s-u)                         (3.10)
  I(r,beta) = 1/2 int_0^inf w^{-beta-1} F_I(r w, -w) dw                 (3.15)
  I = 2 pi delta + (i pi/2)(-r)^beta / sinh(1/2 sum lam) * S(beta)       (3.30)
  S(beta) = sum_n Gamma(n - n/r) / (n! n^beta Gamma(-n/r))
Since beta = -(i/2) sum lambda, (i pi/2)/sinh(pi/2 sum lam) = pi/(2 sin pi beta) = (1/2)Gamma(beta)Gamma(1-beta),
which is the prefactor of (3.22).  [The printed sinh(1/2 sum lam) lacks the pi: harmless misprint.]

Physical prescription: s -> s + i0 <=> w -> w + i0 (r>1 real).  We integrate along a ray
w = t e^{i theta}, theta>0 (contour above all s-channel poles w = n/r) and check theta-independence.
For Re beta < 2 we use the subtracted form, which is the meromorphic continuation in beta:
  I_+(beta) = 1/2 [ int_0^{w0} w^{-b-1}(f-1) dw + int_{w0}^{inf} w^{-b-1} f dw - w0^{-b}/b ].
f-1 is computed from the zeta series of log f near w=0 to avoid cancellation.

R1 corrected closed form (derived in R1_verdicts.md, sec. B4):
  I_+(beta) = (pi / (2 sin pi beta)) * [ (r e^{-i pi})^beta S(beta) - U(beta) ]
  U(beta)  = sum_{m>=1} (-1)^{m-1} m^{-beta} Gamma(1+r m)/(m! Gamma(1+(r-1)m))   (u-channel poles, w=-m)
The u-sum diverges (terms ~ e^{mK}, K = r ln r - (r-1)ln(r-1)) but is alternating; we continue it with the
Lindelof integral  U = (1/2i) int_{c-i inf}^{c+i inf} phi(z)/sin(pi z) dz, c=1/2.
"""
import mpmath as mp
mp.mp.dps = 25

def logf_series(w, r, K=60):
    # log F_I(rw,-w) = sum_{k>=2} zeta(k)/k [ (rw)^k + (-w)^k - ((r-1)w)^k ]
    return mp.fsum(mp.zeta(k)/k*((r*w)**k + (-w)**k - ((r-1)*w)**k) for k in range(2, K))

def f(w, r):
    return mp.exp(mp.loggamma(1 - r*w) + mp.loggamma(1 + w) - mp.loggamma(1 - (r - 1)*w))

def fm1(w, r):
    if abs(w)*r < 0.2:
        return mp.expm1(logf_series(w, r))
    return f(w, r) - 1

def I_direct(r, b, theta=0.25, sign=+1):
    e = mp.exp(1j*sign*theta)
    w0 = 1*e
    g1 = lambda t: (t*e)**(-b-1)*fm1(t*e, r)*e
    g2 = lambda t: (t*e)**(-b-1)*f(t*e, r)*e
    a = mp.quad(g1, [0, 0.05, 0.2, 0.5, 1])
    c = mp.quad(g2, [1, 2, 4, 8, 16, 32, 64, mp.inf])
    return (a + c - w0**(-b)/b)/2

def S_sum(r, b, N=400):
    tot = mp.mpf(0)
    for n in range(1, N+1):
        x = mp.mpf(n)/r
        if abs(x - mp.nint(x)) < mp.mpf(10)**(-15):
            continue   # 1/Gamma(-n/r) = 0
        tot += mp.gamma(n - x)/(mp.factorial(n)*mp.gamma(-x)) * mp.mpf(n)**(-b)
    return tot

def I_paper(r, b, branch):
    """(3.30) with (-r)^beta = r^beta e^{branch*i*pi*beta}."""
    return mp.pi/(2*mp.sin(mp.pi*b)) * r**b*mp.exp(1j*mp.pi*branch*b) * S_sum(r, b)

def U_lindelof(r, b, c=mp.mpf('0.5')):
    phi = lambda z: z**(-b)*mp.exp(mp.loggamma(1 + r*z) - mp.loggamma(1 + z) - mp.loggamma(1 + (r - 1)*z))
    g = lambda y: phi(c + 1j*y)/mp.sin(mp.pi*(c + 1j*y))
    val = mp.quad(g, [-mp.inf, -20, -5, 0, 5, 20, mp.inf])   # dz = i dy
    return val*1j/(2j)

def I_corrected(r, b, sign=+1):
    # sign=+1: +i0 prescription -> (-r)^beta = r^beta e^{-i pi beta}
    return mp.pi/(2*mp.sin(mp.pi*b)) * (r**b*mp.exp(-1j*sign*mp.pi*b)*S_sum(r, b) - U_lindelof(r, b))

if __name__ == "__main__":
    print("== sanity: Res_{s=n} B(-s,1+s/r) numerically vs (3.32) ==")
    r = 2
    for n in [1, 3]:
        eps = mp.mpf('1e-12')
        B = lambda s: mp.gamma(-s)*mp.gamma(1 + s/r)/mp.gamma(1 - s + s/r)
        res_num = eps*B(n + eps)
        x = mp.mpf(n)/r
        print(f"  n={n}: numeric residue {mp.nstr(res_num,10)};  (3.32) gives {mp.nstr(mp.gamma(n-x)/(mp.factorial(n)*mp.gamma(-x)),10)}")
    print("== sanity: U via Lindelof vs closed forms at r=2 (binom(2m,m) generating functions) ==")
    print("  U(0) =", mp.nstr(U_lindelof(2, 0), 15), "  expected 1-1/sqrt5 =", mp.nstr(1 - 1/mp.sqrt(5), 15))
    print("  U(1) =", mp.nstr(U_lindelof(2, 1), 15), "  expected -2 ln(2/(1+sqrt5)) =", mp.nstr(-2*mp.log(2/(1 + mp.sqrt(5))), 15))
    print()
    print("== main table ==")
    for r in [2, mp.mpf('3.7')]:
        for b in [mp.mpf('-1.3'), mp.mpc('-0.5', '0.2'), mp.mpc('-0.25', '-0.4'), mp.mpc(0, '0.7'),
                  mp.mpf('0.5'), mp.mpf('1.0'), mp.mpf('1.5'), mp.mpf('-2.0')]:
            d1 = I_direct(r, b, 0.15); d2 = I_direct(r, b, 0.6); dm = I_direct(r, b, 0.3, sign=-1)
            print(f"r={mp.nstr(r,3)} beta={mp.nstr(b,4)}")
            print(f"   direct(+i0): theta=.15 {mp.nstr(d1,13)}   theta=.6 {mp.nstr(d2,13)}")
            print(f"   direct(-i0):           {mp.nstr(dm,13)}")
            if abs(mp.sin(mp.pi*b)) > 1e-10:
                fp = I_paper(r, b, -1)   # branch that has the right discontinuity for +i0
                print(f"   paper (3.30), (-r)^b=r^b e^(-i pi b): {mp.nstr(fp,13)}   -> mismatch {mp.nstr(d1-fp,6)}")
                cp = I_corrected(r, b, +1); cm = I_corrected(r, b, -1)
                print(f"   R1 corrected (s+u channels), +i0:   {mp.nstr(cp,13)}   -> mismatch {mp.nstr(d1-cp,3)}")
                print(f"   R1 corrected (s+u channels), -i0:   {mp.nstr(cm,13)}   -> mismatch {mp.nstr(dm-cm,3)}")
            else:
                print(f"   paper (3.30): pi/sin(pi beta) pole, residue factor (-r)^b S(b) = {mp.nstr((-r)**b*S_sum(r,b),10)}  (nonzero => spurious pole)")
                print(f"   R1 corrected: residue factor (-r)^b S - U = {mp.nstr((-r)**b*S_sum(r,b) - U_lindelof(r,b),10)}  (should be 0)")
    print()
    print("== residue at beta=0: true I has pole -1/(2 beta) => need (S(0) - U(0)) = -1 ; paper has S(0) ==")
    for r in [2, 3, mp.mpf('3.7'), 10]:
        S0 = S_sum(r, 0); U0 = U_lindelof(r, 0)
        print(f"  r={mp.nstr(r,3)}: S(0)={mp.nstr(S0,12)}  U(0)={mp.nstr(U0,12)}  S-U={mp.nstr(S0-U0,12)}")
