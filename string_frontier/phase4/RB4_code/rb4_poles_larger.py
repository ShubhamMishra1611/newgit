"""
RB4: (i) pole structure in beta of the true I vs the paper's formula,
     (ii) large-r (small-angle) behaviour of the true I vs paper's O(r^{beta-1}) expansion.
"""
import mpmath as mp
from rb4_direct import I_plus, I_minus, I_ray
from rb4_formulas import S_r, S_inv, corrected, paper_330

mp.mp.dps = 30
out = []
def P(*a):
    s = ' '.join(str(x) for x in a)
    print(s, flush=True); out.append(s)

P('# (i) S_r(k), S_{1/r}(k) at integer k; paper formula pi(-r)^b S_r/(2 sin pi b) has a pole at b=k unless S_r(k)=0')
for r in [mp.mpf('2.5'), mp.mpf(7)]:
    for k in [-3, -2, -1, 0, 1, 2, 3]:
        sr = S_r(r, k); si = S_inv(r, k)
        brk = mp.power(r, k) * (-1) ** k * sr + si   # bracket of corrected formula at beta=k
        P(' r=%s k=%+d  S_r=%s  S_1/r=%s   residue(paper)=%s  residue(corrected)=%s' % (
            r, k, mp.nstr(sr, 12), mp.nstr(si, 12),
            mp.nstr(mp.power(r, k) * (-1) ** k * sr / 2, 12), mp.nstr(brk / 2, 12)))
P(' true I: residue at beta=0 is -1/2 (from 1/2 int_0^1 w^{-beta-1} dw), no pole at beta=1 (F-1=O(w^2)),')
P(' residue at beta=2 = -zeta(2) r/2, no poles at beta<0 (exponential decay of F at large |w|).')
for r in [mp.mpf('2.5'), mp.mpf(7)]:
    eps = mp.mpf('1e-12')
    P(' r=%s  eps*I_+(eps-0? ) at beta=+-eps:' % r,
      mp.nstr(eps * I_plus(r, eps), 12), mp.nstr(-eps * I_plus(r, -eps), 12))
    P('   I_+ near beta=1 (1-1e-8, 1+1e-8):', mp.nstr(I_plus(r, 1 - mp.mpf('1e-8')), 12), mp.nstr(I_plus(r, 1 + mp.mpf('1e-8')), 12))
    P('   paper(intended,e^{-i pi b}) near beta=1:', mp.nstr(paper_330(r, 1 - mp.mpf('1e-8'), 'intended', -1), 12))
    P('   I_+ near beta=-2:', mp.nstr(I_plus(r, -2 + mp.mpf('1e-8')), 12), ' paper:', mp.nstr(paper_330(r, -2 + mp.mpf('1e-8'), 'intended', -1), 12))

P('# (ii) large r. True I_+ vs 1/2 Gamma(-beta) (ln r)^beta and vs paper leading term pi(-r)^{beta-1} zeta(beta)/(2 sin pi beta)')
def I_plus_big(r, beta):
    # rotated contour w = i t with extra breakpoints on the scale 1/ln r
    r = mp.mpf(r); beta = mp.mpmathify(beta)
    from rb4_direct import Fm1, Fw
    e = mp.mpc(0, 1)
    f1 = lambda t: t ** (-beta - 1) * Fm1(r, t * e)
    f2 = lambda t: t ** (-beta - 1) * Fw(r, t * e)
    L = mp.log(r)
    pts1 = sorted(set([mp.mpf(0)] + [mp.mpf(10) ** (-k) for k in range(1, 12)] + [1 / L, 2 / L, 4 / L, mp.mpf(1)]))
    pts1 = [p for p in pts1 if p <= 1]
    J = mp.quad(f1, pts1) + mp.quad(f2, [1, 2, 4, 8, 16, 32, 64, mp.inf]) - 1 / beta
    return mp.mpf(1) / 2 * e ** (-beta) * J

for beta in [mp.mpc(0, '0.7'), mp.mpf('-0.5'), mp.mpf('0.5'), mp.mpf('1.5')]:
    P(' beta =', beta)
    for ex in [1, 2, 4, 8, 16, 32]:
        r = mp.mpf(10) ** ex
        Ip = I_plus_big(r, beta)
        asym = mp.gamma(-beta) / 2 * mp.log(r) ** beta
        # refined: 1/2 int w^{-b-1} exp(-w (ln r + gamma - i pi + ln w)) dw  (Regge form F ~ Gamma(1+w) (r w e^{-i pi})^{-w})
        Lc = mp.log(r) + mp.euler - 1j * mp.pi
        refined = None
        if mp.re(beta) < 0:
            refined = mp.quad(lambda t: (t ** (-beta - 1) * mp.gamma(1 + 1j * t) * mp.exp(-1j * t * (mp.log(r * t) + 1j * mp.pi / 2 - 1j * mp.pi))) , [0, 1 / mp.log(r), 1, 10, mp.inf]) * mp.mpf(1) / 2 * (1j) ** (-beta)
        lead = mp.pi * (mp.power(r, beta - 1) * mp.exp(-1j * mp.pi * (beta - 1))) * mp.zeta(beta) / (2 * mp.sin(mp.pi * beta))
        P('   r=1e%-3d I_+=%s  I_+/[G(-b)(ln r)^b/2]=%s  |paper lead r^{b-1}|=%s %s' % (
            ex, mp.nstr(Ip, 10), mp.nstr(Ip / asym, 8), mp.nstr(abs(lead), 4),
            ('' if refined is None else ' Regge-form/I=' + mp.nstr(refined / Ip, 8))))

P('# (iii) decomposition at moderate r: s-channel term vs u-channel term of the corrected formula')
beta = mp.mpc(0, '0.7')
for r in [10, 100, 1000]:
    r = mp.mpf(r)
    st = mp.pi / (2 * mp.sin(mp.pi * beta)) * mp.power(r, beta) * mp.exp(-1j * mp.pi * beta) * S_r(r, beta)
    ut = mp.pi / (2 * mp.sin(mp.pi * beta)) * S_inv(r, beta)
    Ip = I_plus_big(r, beta)
    P(' r=%s s-term=%s u-term=%s sum-I_+=%s  |s|/|I|=%s' % (r, mp.nstr(st, 10), mp.nstr(ut, 10), mp.nstr(st + ut - Ip, 3), mp.nstr(abs(st / Ip), 5)))

open('rb4_poles_larger_out.txt', 'w').write('\n'.join(out) + '\n')
