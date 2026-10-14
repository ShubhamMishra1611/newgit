"""RB4: growth wedge of F(rw,-w) around the negative w axis (where the u-channel poles sit).
Prediction from Stirling: on w = R e^{i(pi-phi)}, |F| ~ exp(R (K cos phi - pi sin phi)), K = r ln r - (r-1) ln(r-1),
so F grows for |phi| < arctan(K/pi) and decays outside."""
import mpmath as mp
from rb4_direct import Fw
mp.mp.dps = 20
for r in [2, 10]:
    r = mp.mpf(r)
    K = r * mp.log(r) - (r - 1) * mp.log(r - 1)
    print('r=%s K=%s  arctan(K/pi)=%s deg   arctan(pi/K)=%s deg' % (r, mp.nstr(K, 5), mp.nstr(mp.degrees(mp.atan(K / mp.pi)), 4), mp.nstr(mp.degrees(mp.atan(mp.pi / K)), 4)))
    for phideg in [10, 20, 30, 45, 60, 75]:
        phi = mp.radians(phideg)
        vals = []
        for R in [20, 40, 80]:
            w = R * mp.exp(1j * (mp.pi - phi))
            vals.append(mp.log(abs(Fw(r, w))) / R)
        print('   phi=%2d deg  (1/R) ln|F| at R=20,40,80: %s   Stirling slope K cos - pi sin = %s' % (
            phideg, [mp.nstr(v, 4) for v in vals], mp.nstr(K * mp.cos(phi) - mp.pi * mp.sin(phi), 4)))
