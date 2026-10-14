"""
RB4 (supplementary to claim c): distributional large-r limit on the principal series.

Pair the full I(r, beta = i y) (y = Im beta = -(1/2) sum lambda) with phi(y) = exp(-y^2/2):
   <I_r, phi> = (1/2) int_0^oo dw/w F(rw + i0, -w) phihat(ln w),   phihat(t) = sqrt(2 pi) exp(-t^2/2).
(The distribution on Re beta = 0 is the boundary value from Re beta < 0, where the w-integral converges.)
Rotate the contour to w = i t: legitimate because phihat(ln w) is entire in ln w and F decays in the
first quadrant.  For F = 1 (Yang-Mills) the pairing is exactly pi, i.e. I_YM = pi delta(y) = 2 pi delta(sum lambda).
The claim 'I -> 2 pi delta(sum lambda) + O(r^{beta-1})' (1806.05688 text; 2405.01643 eq. 3.7) predicts
<I_r, phi> -> pi as r -> oo.
"""
import mpmath as mp
from rb4_direct import Fw, I_ray

mp.mp.dps = 20
out = []
def P(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); out.append(s)


def pair(r):
    r = mp.mpf(r)
    g = lambda t: Fw(r, 1j * t) * mp.sqrt(2 * mp.pi) * mp.exp(-(mp.log(t) + 1j * mp.pi / 2) ** 2 / 2) / t
    L = mp.log(r)
    pts = sorted(set([mp.mpf(0)] + [mp.mpf(10) ** k for k in range(-14, 4)] + [1 / L, 2 / L, 4 / L])) + [mp.inf]
    return mp.quad(g, pts) / 2


g1 = lambda t: mp.sqrt(2 * mp.pi) * mp.exp(-(mp.log(t) + 1j * mp.pi / 2) ** 2 / 2) / t
P('sanity, F=1 (YM):', mp.nstr(mp.quad(g1, [0, 1e-6, 1e-3, 1, 1e3, mp.inf]) / 2, 12), ' (pi =', mp.nstr(mp.pi, 12), ')')
for ex in ['0.5', 1, 2, 4, 8, 16, 32, 64]:
    mp.mp.dps = 25 + 2 * int(mp.ceil(mp.mpf(ex)))   # precision must exceed log10(r) digits
    r = mp.mpf(10) ** mp.mpf(ex)
    val = pair(r)
    P(' r=1e%-4s <I_r,phi> = %s   |<I_r,phi>|/pi = %s' % (ex, mp.nstr(val, 10), mp.nstr(abs(val) / mp.pi, 6)))

# independent cross-check at r=10: integrate I_+(r, b0+iy) phi(y - i b0) over y on Re beta = b0 = -0.3
mp.mp.dps = 12
r = mp.mpf(10); b0 = mp.mpf('-0.3')
f = lambda y: I_ray(r, b0 + 1j * y, mp.pi / 2) * mp.exp(-(y - 1j * b0) ** 2 / 2)
chk = mp.quad(f, [-9, -4, -1.5, 0, 1.5, 4, 9])
P(' cross-check r=10 via shifted y-contour int dy I_+(r,-0.3+iy) phi(y+0.3i):', mp.nstr(chk, 8))
open('rb4_distributional_out.txt', 'w').write('\n'.join(out) + '\n')
