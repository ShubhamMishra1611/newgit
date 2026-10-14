"""
R1: large-r (forward) behaviour of the celestial type-I integral I(r,beta) of 1806.05688 (3.15).
Claim in 1806.05688: I - 2 pi delta is a small-angle series in r^{-1}: leading term
    (1/2) Gamma(beta)Gamma(1-beta) (-r)^{beta-1} zeta(beta)          [(3.22), k=1, S(beta-1,1)=zeta(beta)]
R1 finding: the missing u-channel term -(pi/(2 sin pi beta)) U(beta) behaves like
    -(pi/(2 sin pi beta)) (ln r)^beta / Gamma(1+beta)  = -Gamma(1-beta) (ln r)^beta / (2 beta)
(leading Hankel-contour asymptotics of the Lindelof integral), which dominates r^{beta-1} for Re beta < 1.
Here we test numerically with the direct (+i0) integral.
"""
import mpmath as mp
from b4_celestial_type1 import I_direct, U_lindelof
mp.mp.dps = 20

for b in [mp.mpc('-0.5', '0.3'), mp.mpc(0, '0.7')]:
    print(f"beta = {mp.nstr(b,4)}")
    for r in [10, 100, 1000, 10000]:
        d = I_direct(r, b, 0.3)
        lead_paper = 0.5*mp.gamma(b)*mp.gamma(1 - b)*(r*mp.exp(-1j*mp.pi))**(b - 1)*mp.zeta(b)
        U = U_lindelof(r, b)
        uterm = -mp.pi/(2*mp.sin(mp.pi*b))*U
        uasym = -mp.gamma(1 - b)*mp.log(r)**b/(2*b)
        print(f"  r={r:6d}: direct={mp.nstr(d,8)}  |direct|={mp.nstr(abs(d),5)}"
              f"  paper k=1 term={mp.nstr(lead_paper,5)} (|.|={mp.nstr(abs(lead_paper),3)})"
              f"  u-term={mp.nstr(uterm,6)}  (ln r)^b asym={mp.nstr(uasym,6)}"
              f"  direct-uterm={mp.nstr(d-uterm,5)}")

print()
print("Regge-type asymptotics: f(w) ~ Gamma(1+w)(r w e^{-i pi})^{-w}  =>  I ~ (1/2) int w^{-b-1} Gamma(1+w) (r w e^{-i pi})^{-w} dw")
for b in [mp.mpc('-0.5', '0.3'), mp.mpc(0, '0.7')]:
    for r in [100, 10000, 10**6]:
        d = I_direct(r, b, 0.3)
        e = mp.exp(1j*0.3)
        g = lambda t: (t*e)**(-b-1)*(mp.gamma(1 + t*e)*(r*t*e*mp.exp(-1j*mp.pi))**(-t*e) - (1 if abs(t) < 1 else 0))*e
        reg = (mp.quad(g, [0, 0.01, 0.1, 1, 10, mp.inf]) - e**(-b)/b)/2
        lead = 0.5*mp.gamma(-b)*(mp.log(r) - mp.log(mp.log(r)) - 1j*mp.pi)**b
        print(f"  beta={mp.nstr(b,3)} r=1e{int(mp.log10(r))}: direct={mp.nstr(d,6)}  Regge-integral={mp.nstr(reg,6)}  (1/2)Gamma(-b)(ln r - lnln r - i pi)^b={mp.nstr(lead,6)}")
