"""Low-energy coefficients of the Cheung-Hillman-Remmen VS deformation, arXiv:2408.03362v1 eqs.(2),(5),(9),(10),(13).
Heavy part only (rho=0 scalar interpretation). Units mu(1)=1, so M (first massive state) = 1.
Checks: (i) closed forms (10) vs direct dual-resonant sum; (ii) crossing: linear-in-s,t coefficients vanish;
(iii) lambda=1 -> VS: g0 = 2 zeta(3); (iv) lambda->0 limit -> scalar exchange g3/g2=3, g4/g2=1/2."""
import mpmath as mp
mp.mp.dps = 20

def closed(lam):
    a = (1+lam)/2; b = (1+3*lam)/2
    g0 = 2*mp.hyper([1, lam, a, a], [1+lam, b, b], 1)
    g2 = mp.hyper([1, lam, lam, lam, a, a], [1+lam, 1+lam, 1+lam, b, b], 1)
    f = lambda nu: mp.hyper([1, lam, lam, lam, (1+nu)/2, (1+nu)/2], [1+lam, 1+lam, 1+lam, b, b], 1)
    g3 = -4*lam*mp.diff(f, lam) + 3*mp.hyper([1, lam, lam, lam, lam, a, a], [1+lam]*4 + [b, b], 1)   # d/dlog nu = nu d/dnu
    g4 = mp.mpf(1)/2*mp.hyper([1] + [lam]*5 + [a, a], [1+lam]*5 + [b, b], 1)
    return g0, g2, g3, g4

def direct(lam, N=4000):
    """sum_{n>=1} R(n,t)[1/(mu-s)+1/(mu+s+t)], Taylor coefficients c[a][b] of s^a t^b up to total degree 4."""
    lam = mp.mpf(lam); a0 = (1+lam)/2; b0 = (1+3*lam)/2
    C = [[mp.mpf(0)]*5 for _ in range(5)]
    # running log-derivative sums for P(t) = (a0+lam t)_{n-1}/(b0)_{n-1}
    logP0 = mp.mpf(0); d = [mp.mpf(0)]*5   # d[j] = j-th t-derivative of log P at t=0
    for n in range(1, N+1):
        if n >= 2:
            k = n-2
            logP0 += mp.log((a0+k)/(b0+k))
            for j in range(1, 5):
                d[j] += (-1)**(j-1)*mp.factorial(j-1)*lam**j/(a0+k)**j
        # Taylor of R = exp(2 log P): coefficients r_0..r_4
        L = [2*logP0] + [2*d[j]/mp.factorial(j) for j in range(1, 5)]   # log R = sum L_j t^j
        # exp of series
        r = [mp.e**L[0], 0, 0, 0, 0]
        for m in range(1, 5):
            r[m] = sum(j*L[j]*r[m-j] for j in range(1, m+1))/m
        mu = (n+lam-1)/lam
        # 1/(mu-s) = sum s^a/mu^{a+1};  1/(mu+s+t) = sum_{p} (-1)^p (s+t)^p/mu^{p+1}
        K = [[mp.mpf(0)]*5 for _ in range(5)]
        for a in range(5):
            K[a][0] += 1/mu**(a+1)
        for p in range(5):
            for i in range(p+1):
                K[i][p-i] += (-1)**p*mp.binomial(p, i)/mu**(p+1)
        for a in range(5):
            for b in range(5-a):
                C[a][b] += sum(K[a][b-j]*r[j] for j in range(b+1))
    return C

for lam in [mp.mpf('0.5'), mp.mpf(1), mp.mpf(2), mp.mpf(5)]:
    g0, g2, g3, g4 = closed(lam)
    print("lambda=%s  closed: g0=%s g2=%s g3=%s g4=%s | g3/g2=%s g4/g2=%s" % (
        mp.nstr(lam, 3), *(mp.nstr(x, 10) for x in (g0, g2, g3, g4)), mp.nstr(g3/g2, 8), mp.nstr(g4/g2, 8)))
print("2 zeta(3) =", mp.nstr(2*mp.zeta(3), 12))
for lam in [1, 2]:
    C = direct(lam, N=3000)
    print("direct lambda=%d (N=3000): [s]=%s [t]=%s  g0~%s g2=[s^2]/2=%s g3=-[s^2 t]=%s g4=[s^4]/4=%s  check g2 via [t^2]/2=%s" % (
        lam, mp.nstr(C[1][0], 6), mp.nstr(C[0][1], 6), mp.nstr(C[0][0], 8), mp.nstr(C[2][0]/2, 10),
        mp.nstr(-C[2][1], 10), mp.nstr(C[4][0]/4, 10), mp.nstr(C[0][2]/2, 10)))
