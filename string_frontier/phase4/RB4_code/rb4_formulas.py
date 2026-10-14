"""
RB4: the paper's pole-sum formulas (1806.05688 eqs. 3.30, 3.31, 3.32) vs direct integral,
and my own derivation of the missing u-channel term.

Paper (3.30), as printed (beta = -(i/2) sum lambda, so (1/2) sum lambda = i beta):
   I = 2 pi delta + (i pi/2) (-r)^beta  sinh(i beta)^{-1}  S_r(beta)
   S_r(beta) = sum_{n>=1} Gamma(n - n/r) / ( n! n^beta Gamma(-n/r) )
 'printed'   : sinh(i beta) = i sin(beta)
 'intended'  : sinh(i pi beta) = i sin(pi beta)   (only choice consistent with 3.31)
Paper (3.31): I = 2 pi delta + i pi/(1 - e^{-2 pi i beta}) sum_n Res_{s=n}{ (s/r)^{-beta} B(-s, 1+s/r) }
Paper (3.32): Res_{s=n} B(-s,1+s/r) = + Gamma(n-n/r)/(n! Gamma(-n/r))   [sign checked below]

My derivation (keyhole contour; branch cut of w^{-beta-1} on the positive w axis):
   2 I_+ - 2 e^{-2 pi i beta} I_-  = 2 pi i sum_{w=-m} Res  - (big-circle term)
   I_+ - I_- = - pi i sum_{w=n/r} Res          (s-channel poles, exact, convergent)
 => formally  I_+ = pi/(2 sin pi beta) [ (r e^{-i pi})^beta S_r(beta) + S_{1/r}(beta) ]
   with S_{1/r}(beta) = sum_{n>=1} (-1)^n binom(r n, n) n^{-beta}  (u-channel poles at w=-m)
 The u-channel sum diverges like C^n, C = r^r/(r-1)^(r-1) > 1, and the big circle does not
 vanish (F grows in a wedge around the negative w axis), so S_{1/r} must be defined by
 analytic continuation.  I define it by Abel/Mellin continuation of the generating function
   G(z) = sum_{n>=0} binom(r n, n) z^n = 1/(1 - r + r/B(z)),   B = 1 + z B^r   (GKP 5.61)
   S_{1/r}(beta) = 1/Gamma(beta) int_0^oo t^{beta-1} [G(-e^{-t}) - 1] dt   (continued in beta)
 and TEST whether the formula then reproduces the direct integral.
"""
import mpmath as mp
from rb4_direct import I_plus, I_minus, Fw

mp.mp.dps = 30


def S_r(r, beta, N=None):
    r = mp.mpf(r); beta = mp.mpmathify(beta)
    tot = mp.mpc(0)
    n = 1
    small = 0
    while True:
        a = -mp.mpf(n) / r
        # Gamma(n+a)/Gamma(a) = rf(a, n) (rising factorial) -> finite even if a is a pole
        term = mp.rf(a, n) / mp.factorial(n) * mp.power(n, -beta)
        tot += term
        if abs(term) < mp.mpf(10) ** (-mp.mp.dps - 5) * (1 + abs(tot)):
            small += 1
            if small > 30:
                break
        else:
            small = 0
        n += 1
        if N and n > N:
            break
    return tot


def paper_330(r, beta, variant='intended', branch=+1):
    """non-delta part of (3.30); (-r)^beta = r^beta e^{branch*i pi beta}"""
    beta = mp.mpmathify(beta)
    mr = mp.power(r, beta) * mp.exp(branch * 1j * mp.pi * beta)
    if variant == 'printed':
        den = mp.sinh(1j * beta)
    else:
        den = mp.sinh(1j * mp.pi * beta)
    return (1j * mp.pi / 2) * mr / den * S_r(r, beta)


def res_B_numeric(r, n):
    """Res_{s=n} B(-s, 1+s/r) by a small numerical contour integral."""
    r = mp.mpf(r)
    B = lambda s: mp.gamma(-s) * mp.gamma(1 + s / r) / mp.gamma(1 - s + s / r)
    rad = mp.mpf('1e-3')
    f = lambda ph: B(n + rad * mp.exp(1j * ph)) * rad * mp.exp(1j * ph) * 1j
    return mp.quad(f, [0, mp.pi, 2 * mp.pi]) / (2j * mp.pi)


def paper_331(r, beta, sign_332=+1):
    beta = mp.mpmathify(beta); r = mp.mpf(r)
    tot = mp.mpc(0)
    n = 1
    small = 0
    while True:
        term = (mp.mpf(n) / r) ** (-beta) * sign_332 * mp.rf(-mp.mpf(n) / r, n) / mp.factorial(n)
        tot += term
        small = small + 1 if abs(term) < mp.mpf(10) ** (-40) else 0   # some terms vanish exactly (integer n/r)
        if small > 30:
            break
        n += 1
    return 1j * mp.pi / (1 - mp.exp(-2j * mp.pi * beta)) * tot


# ---------- u-channel continuation ----------
def Bgen(r, z):
    # root of B = 1 + z B^r continuously connected to B(0)=1; for z in (-1,0): B in (0,1)
    if mp.im(z) == 0 and -1 <= z < 0:
        # unique root in (0,1); bracketing solver is robust for large r
        return mp.findroot(lambda B: B - 1 - z * B ** r, (mp.mpf(0), mp.mpf(1)), solver='anderson')
    return mp.findroot(lambda B: B - 1 - z * B ** r, mp.mpf(1) / (1 - z))


def G(r, z):
    B = Bgen(r, z)
    return 1 / (1 - r + r / B)


# --- truncated power series in t (lists of mpf), used for the Taylor data of h(t) ---
def ps_mul(a, b, N):
    return [mp.fsum(a[k] * b[n - k] for k in range(n + 1)) for n in range(N + 1)]


def ps_div(a, d, N):
    q = []
    for n in range(N + 1):
        q.append((a[n] - mp.fsum(d[k] * q[n - k] for k in range(1, n + 1))) / d[0])
    return q


def ps_log(a, N):
    L = [mp.log(a[0])]
    for n in range(1, N + 1):
        L.append((a[n] - mp.fsum(k * L[k] * a[n - k] for k in range(1, n)) / n) / a[0])
    return L


def ps_exp(L, N):
    E = [mp.exp(L[0])]
    for n in range(1, N + 1):
        E.append(mp.fsum(k * L[k] * E[n - k] for k in range(1, n + 1)) / n)
    return E


def h_taylor(r, N):
    """Taylor coefficients at t=0 of h(t) = G(-e^{-t}) - 1, G = B/(B(1-r)+r), B - 1 + e^{-t} B^r = 0."""
    r = mp.mpf(r)
    emt = [(-1) ** j / mp.factorial(j) for j in range(N + 1)]
    b0 = mp.findroot(lambda B: B - 1 + B ** r, (mp.mpf(0), mp.mpf(1)), solver='anderson')
    B = [b0] + [mp.mpf(0)] * N
    for it in range(12):  # Newton on series; order doubles each step
        Br = ps_exp([r * c for c in ps_log(B, N)], N)
        Br1 = ps_div(Br, B, N)
        Phi = [B[n] - (1 if n == 0 else 0) + ps_mul(emt, Br, N)[n] for n in range(N + 1)]
        dPhi = [(1 if n == 0 else 0) + r * x for n, x in enumerate(ps_mul(emt, Br1, N))]
        step = ps_div(Phi, dPhi, N)
        B = [B[n] - step[n] for n in range(N + 1)]
        if max(abs(x) for x in step) < mp.mpf(10) ** (-mp.mp.dps + 3):
            break
    den = [B[n] * (1 - r) + (r if n == 0 else 0) for n in range(N + 1)]
    Gs = ps_div(B, den, N)
    Gs[0] -= 1
    return Gs


def S_inv(r, beta, J=80, T0=mp.mpf(1)):
    """S_{1/r}(beta) := continuation of sum_{n>=1} (-1)^n binom(rn,n) n^{-beta}
    = 1/Gamma(beta) int_0^oo t^{beta-1} h(t) dt,  h(t) = G(-e^{-t}) - 1  (Re beta > 0),
    continued to all beta: on [0,T0] integrate the Taylor series of h at t=0 term by term
    (radius = |ln C + i pi| >= pi > T0)."""
    r = mp.mpf(r); beta = mp.mpmathify(beta)
    key = (r, J, mp.mp.dps)
    if key not in _hc_cache:
        _hc_cache[key] = h_taylor(r, J)
    hc = _hc_cache[key]
    if mp.im(beta) == 0 and mp.re(beta) <= 0 and mp.re(beta) == int(mp.re(beta)):
        k = int(-mp.re(beta))
        return (-1) ** k * mp.factorial(k) * hc[k]
    h = lambda t: G(r, -mp.exp(-t)) - 1
    Jv = mp.fsum(hc[j] * T0 ** (beta + j) / (beta + j) for j in range(J + 1))
    Jv += mp.quad(lambda t: t ** (beta - 1) * h(t), [T0, 2, 4, 8, 16, 32, 64, mp.inf])
    return Jv * mp.rgamma(beta)


_hc_cache = {}


def corrected(r, beta, sign=+1):
    """I_{+} (sign=+1) or I_{-} (sign=-1)"""
    beta = mp.mpmathify(beta)
    return mp.pi / (2 * mp.sin(mp.pi * beta)) * (
        mp.power(r, beta) * mp.exp(-sign * 1j * mp.pi * beta) * S_r(r, beta) + S_inv(r, beta))


if __name__ == '__main__':
    out = []
    def P(*a):
        s = ' '.join(str(x) for x in a)
        print(s); out.append(s)

    P('# 0. generating function check: G(z) vs direct series at z=-0.05, r=2.5')
    r = mp.mpf('2.5'); z = mp.mpf('-0.05')
    ser = mp.nsum(lambda n: mp.binomial(r * n, n) * z ** n, [0, mp.inf])
    P('series', mp.nstr(ser, 25), ' closed', mp.nstr(G(r, z), 25))
    P('# 0b. S_inv(beta=-k) vs (z d/dz)^k G at z=-1 (exact for integer beta<=0), r=2.5')
    for k in [0, 1, 2]:
        zd = mp.diff(lambda x: G(r, -mp.exp(x)), 0, k)  # (z d/dz)^k = d^k/dx^k with z=-e^x
        P(' k=%d  S_inv=%s   (zd/dz)^k G - delta_k0 =%s' % (k, mp.nstr(S_inv(r, -k), 20), mp.nstr(zd - (1 if k == 0 else 0), 20)))

    P('# 1. residue sign in (3.32): Res_{s=n} B(-s,1+s/r) numeric vs +Gamma(n-n/r)/(n!Gamma(-n/r))')
    for (r, n) in [(mp.mpf('2.5'), 1), (mp.mpf('2.5'), 2), (mp.mpf('3.3'), 3)]:
        num = res_B_numeric(r, n)
        pap = mp.rf(-n / r, n) / mp.factorial(n)
        P(' r=%s n=%d numeric=%s  paper(3.32)=%s' % (r, n, mp.nstr(num, 15), mp.nstr(pap, 15)))

    P('# 2. exact s-channel discontinuity: I_+ - I_- vs -pi i sum Res_w = -pi i r^beta S_r(beta)')
    for (r, b) in [(2, -1.3), (3.7, -0.4 + 1.1j), (5, 0.7j), (2.5, 1.5)]:
        d = I_plus(r, b) - I_minus(r, b)
        pr = -1j * mp.pi * mp.power(r, b) * S_r(r, b)
        P(' r=%s beta=%s  I+-I- =%s  pred=%s  |diff|=%s' % (r, b, mp.nstr(d, 15), mp.nstr(pr, 15), mp.nstr(abs(d - pr), 3)))

    P('# 3. paper (3.30)/(3.31) vs direct integral')
    pts = [(2, -1.3), (2, -0.5), (3.7, -0.4 + 1.1j), (5, 0.7j), (2.5, 0.3), (10, -0.25), (1.5, 1.5), (4, -2.7)]
    for (r, b) in pts:
        Ip = I_plus(r, b); Im = I_minus(r, b)
        P(' r=%s beta=%s' % (r, b))
        P('    direct I_+ =', mp.nstr(Ip, 15), '   direct I_- =', mp.nstr(Im, 15))
        for var in ['printed', 'intended']:
            for br in [+1, -1]:
                v = paper_330(r, b, var, br)
                P('    (3.30) %-8s (-r)^b=r^b e^{%+d i pi b}: %s  |-I+|=%s |-I-|=%s' % (
                    var, br, mp.nstr(v, 12), mp.nstr(abs(v - Ip), 3), mp.nstr(abs(v - Im), 3)))
        for sg in [+1, -1]:
            v = paper_331(r, b, sg)
            P('    (3.31) residue sign %+d: %s  |-I+|=%s |-I-|=%s' % (sg, mp.nstr(v, 12), mp.nstr(abs(v - Ip), 3), mp.nstr(abs(v - Im), 3)))
        cp = corrected(r, b, +1); cm = corrected(r, b, -1)
        P('    corrected (s+u):  I_+ %s |diff|=%s ;  I_- %s |diff|=%s' % (mp.nstr(cp, 15), mp.nstr(abs(cp - Ip), 3), mp.nstr(cm, 15), mp.nstr(abs(cm - Im), 3)))
    open('rb4_formulas_out.txt', 'w').write('\n'.join(out) + '\n')
