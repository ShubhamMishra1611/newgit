"""
RB4 independent check of Stieberger-Taylor arXiv:1806.05688, Sec. 3.

Definitions (read from the fetched paper text, eqs. 3.2, 3.3, 3.10, 3.15):
  r = z12 z34/(z23 z41) > 1 (physical), 1/r = sin^2(theta/2) = -u/s
  F_I(s,u) = Gamma(1-s)Gamma(1-u)/Gamma(1-s-u),  s = alpha' s12 > 0,  u = -s/r
  I(r,beta) = 1/2 int_0^oo w^{-beta-1} F_I(r w, -w) dw,  beta = -(i/2) sum lambda_n
s-channel poles at w = n/r (on the integration path) -> i eps prescription.
I_plus  : s -> s + i0  (contour passes ABOVE the poles), equivalently r -> r + i0
I_minus : s -> s - i0

Method A: rotate the w contour onto the ray w = t e^{i theta}.  For I_plus any
          theta in (0, pi - theta_c) is allowed (F decays in the first quadrant
          and in most of the second; checked numerically by theta-independence).
          Analytic continuation to Re beta < 2 by subtracting F(0)=1 on [0,1]
          (F - 1 = zeta(2) r w^2 + O(w^3)).
Method B: the paper's own x-representation (3.18), valid for r<0 and continued
          to r>1 via Im r>0:   I_plus = -(r/2) Gamma(1-beta) int_0^1 dx/x (A - i0)^{beta-1},
          A = r ln x - ln(1-x); A has a simple zero at x* in (0,1); the path is
          deformed below x* (Im A < 0).  Valid for Re beta < 0.
"""
import mpmath as mp

mp.mp.dps = 30


def F(s, u):
    # Gamma(1-s)Gamma(1-u)/Gamma(1-s-u) via loggamma (principal branches are fine:
    # exp(loggamma) = Gamma on the whole plane)
    return mp.exp(mp.loggamma(1 - s) + mp.loggamma(1 - u) - mp.loggamma(1 - s - u))


def Fw(r, w):
    return F(r * w, -w)


def Fm1(r, w):
    """F(rw,-w) - 1 without cancellation: for |w| small use
    ln F = sum_{k>=2} zeta(k)/k (s^k + u^k - (s+u)^k), then expm1."""
    if abs(r * w) < mp.mpf('0.05'):
        L = mp.mpc(0)
        k = 2
        while True:
            term = mp.zeta(k) / k * ((r * w) ** k + (-w) ** k - ((r - 1) * w) ** k)
            L += term
            if abs(term) < mp.mpf(10) ** (-mp.mp.dps - 10) * abs(L) and k > 3:
                break
            k += 1
        return mp.expm1(L)
    return Fw(r, w) - 1


def I_ray(r, beta, theta):
    """1/2 int along w = t e^{i theta}, continued to Re beta < 2 (beta != 0)."""
    r = mp.mpf(r) if not isinstance(r, mp.mpc) else r
    beta = mp.mpmathify(beta)
    e = mp.exp(1j * theta)
    pref = mp.mpf(1) / 2 * e ** (-beta)  # (t e^{i th})^{-beta-1} e^{i th} dt
    f1 = lambda t: t ** (-beta - 1) * Fm1(r, t * e)
    f2 = lambda t: t ** (-beta - 1) * Fw(r, t * e)
    pts1 = [0, mp.mpf('0.25'), mp.mpf('0.5'), 1]
    pts2 = [1, 2, 4, 8, 16, 32, 64, 128, mp.inf]
    J = mp.quad(f1, pts1) + mp.quad(f2, pts2) - 1 / beta
    return pref * J


def I_plus(r, beta, theta=mp.pi / 2):
    return I_ray(r, beta, theta)


def I_minus(r, beta, theta=mp.pi / 2):
    return I_ray(r, beta, -theta)


def I_xrep_plus(r, beta, rho=None):
    """Method B (paper eq. 3.18 continued), Re beta < 0 only."""
    r = mp.mpf(r)
    beta = mp.mpmathify(beta)
    A = lambda x: r * mp.log(x) - mp.log(1 - x)
    xs = mp.findroot(A, mp.mpf('0.5'))
    if rho is None:
        rho = min(xs, 1 - xs) / 4

    def pw(a):
        # (A - i0)^{beta-1}, branch: arg in (-pi, 0] continuous from below
        if mp.im(a) == 0 and mp.re(a) < 0:
            lg = mp.log(-a) - 1j * mp.pi
        else:
            lg = mp.log(a)
        return mp.exp((beta - 1) * lg)

    g = lambda x: pw(A(x)) / x
    # real pieces
    # x < x*: substitute x = e^{-t} so the slow algebraic tail t^{beta-1} is integrated to t = oo
    gt = lambda t: pw(A(mp.exp(-t)))
    t0 = -mp.log(xs - rho)
    J = mp.quad(gt, [t0, t0 + 1, t0 + 4, t0 + 16, t0 + 64, mp.inf])
    J += mp.quad(g, [xs + rho, (1 + xs) / 2, 1 - (1 - xs) / 8, 1])
    # semicircle below x*: x = xs + rho e^{i phi}, phi from -pi to 0
    h = lambda ph: g(xs + rho * mp.exp(1j * ph)) * 1j * rho * mp.exp(1j * ph)
    J += mp.quad(h, [-mp.pi, -mp.pi / 2, 0])
    return -(r / 2) * mp.gamma(1 - beta) * J


if __name__ == '__main__':
    out = []
    def P(*a):
        s = ' '.join(str(x) for x in a)
        print(s); out.append(s)
    P('# Method A: theta-independence of rotated contour (I_plus), r=2, beta=-1.3')
    for th in [mp.mpf('0.05'), mp.mpf('0.3'), mp.pi / 4, mp.pi / 2, mp.mpf(2.0), mp.mpf(2.4)]:
        P('theta=%s' % mp.nstr(th, 5), mp.nstr(I_ray(2, -1.3, th), 20))
    P('# Method B (x-representation 3.18 with A - i0):')
    P('I_plus(x-rep) r=2 beta=-1.3', mp.nstr(I_xrep_plus(2, -1.3), 20))
    P('# I_minus (theta=-pi/2):', mp.nstr(I_minus(2, -1.3), 20))
    for (r, b) in [(2, -0.5), (3.7, -0.4 + 1.1j), (1.3, -2.2), (10, -0.25)]:
        a = I_plus(r, b); bb = I_xrep_plus(r, b)
        P('r=%s beta=%s  A:%s  B:%s  |A-B|=%s' % (r, b, mp.nstr(a, 18), mp.nstr(bb, 18), mp.nstr(abs(a - bb), 3)))
    open('rb4_direct_out.txt', 'w').write('\n'.join(out) + '\n')
