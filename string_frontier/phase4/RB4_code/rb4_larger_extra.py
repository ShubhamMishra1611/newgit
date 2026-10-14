"""
RB4 extra: (iii) s-/u-channel decomposition at large r (principal series beta=0.7i);
(iv) distributional large-r limit on the principal series: pair the full I(r, beta=iy)
     with a Gaussian test function phi(y)=exp(-y^2/2).  The physical distribution is the
     boundary value from Re beta<0 (where the w-integral converges), so
       <I_r, phi> = int dy I(r, b0 + i y) phi(y - i b0)   with b0 = -0.3 (no poles in -0.3<=Re b<0)
     Paper's claim I -> 2 pi delta(sum lambda) = pi delta(y) predicts <I_r,phi> -> pi.
"""
import mpmath as mp
from rb4_direct import Fm1, Fw
from rb4_formulas import S_r, S_inv
mp.mp.dps = 20
out = []
def P(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); out.append(s)

def I_plus_big(r, beta):
    r = mp.mpf(r); beta = mp.mpmathify(beta)
    e = mp.mpc(0, 1)
    f1 = lambda t: t ** (-beta - 1) * Fm1(r, t * e)
    f2 = lambda t: t ** (-beta - 1) * Fw(r, t * e)
    L = mp.log(r)
    pts1 = sorted(set([mp.mpf(0)] + [mp.mpf(10) ** (-k) for k in range(1, 12)] + [1 / L, 2 / L, 4 / L, mp.mpf(1)]))
    pts1 = [p for p in pts1 if p <= 1]
    J = mp.quad(f1, pts1) + mp.quad(f2, [1, 2, 4, 8, 16, 32, 64, mp.inf]) - 1 / beta
    return mp.mpf(1) / 2 * e ** (-beta) * J

P('# (iii) decomposition, beta=0.7i')
beta = mp.mpc(0, '0.7')
for r in [10, 100, 1000]:
    r = mp.mpf(r)
    st = mp.pi / (2 * mp.sin(mp.pi * beta)) * mp.power(r, beta) * mp.exp(-1j * mp.pi * beta) * S_r(r, beta)
    ut = mp.pi / (2 * mp.sin(mp.pi * beta)) * S_inv(r, beta)
    Ip = I_plus_big(r, beta)
    P(' r=%s s-term(paper,branch e^{-i pi b})=%s u-term=%s |s+u-I_+|=%s' % (r, mp.nstr(st, 10), mp.nstr(ut, 10), mp.nstr(abs(st + ut - Ip), 3)))

P('# (iv) <I_r, phi>, phi=exp(-y^2/2); YM/paper prediction for r->oo: pi =', mp.nstr(mp.pi, 8))
b0 = mp.mpf('-0.3')
for ex in [1, 2, 4, 8, 16]:
    r = mp.mpf(10) ** ex
    f = lambda y: I_plus_big(r, b0 + 1j * y) * mp.exp(-(y - 1j * b0) ** 2 / 2)
    val = mp.quad(f, [-12, -6, -3, -1, 0, 1, 3, 6, 12])
    # crude prediction of the weight left in the 'YM region' w < 1/ln r: (1/2) int_{-oo}^{-ln ln r} phi~(t) dt
    pred = mp.pi * mp.ncdf(-mp.log(mp.log(r)))
    P(' r=1e%-3d <I_r,phi> = %s   (pi*Phi(-ln ln r) = %s)' % (ex, mp.nstr(val, 10), mp.nstr(pred, 6)))
open('rb4_larger_extra_out.txt', 'w').write('\n'.join(out) + '\n')
