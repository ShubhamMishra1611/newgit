"""
R2 independent checks of D1 (KKLT toy, hep-th/0301240 Sec. III), D2 (KPV hep-th/0112197 eqs. 40,45; BDGL 1809.06861),
D3 (g_s M^2 of the five Table-2 vacua of 2406.13751).
V_F derived here from e^K (K^{rho rhobar}|D W|^2 - 3|W|^2), K=-3 ln(-i(rho-rhobar)), W=W0+A e^{i a rho}, rho=i sigma:
  V_F = a^2A^2 e^{-2a s}/(6 s) + a A e^{-a s} W/(2 s^2),  W = W0 + A e^{-a s}  (equals the KKLT Sec. IIIB formula).
"""
import mpmath as mp
mp.mp.dps = 30
W0, A, a = mp.mpf('-1e-4'), mp.mpf(1), mp.mpf('0.1')

def VF(s):
    e = mp.exp(-a*s); W = W0 + A*e
    return a**2*A**2*e**2/(6*s) + a*A*e*W/(2*s**2)

# -------- D1 --------
scr = mp.findroot(lambda s: W0 + A*mp.exp(-a*s)*(1 + 2*a*s/3), 110)
VAdS = -a**2*A**2*mp.exp(-2*a*scr)/(6*scr)
print("D1: sigma_cr =", mp.nstr(scr, 8), " V_AdS(eq.14) =", mp.nstr(VAdS, 6), " V_F(sigma_cr) =", mp.nstr(VF(scr), 6))

import numpy as np
from scipy.optimize import brentq
W0f, Af, af = -1e-4, 1.0, 0.1
def dVF_f(x):
    e = np.exp(-af*x); W = W0f + Af*e
    dT1 = af**2*Af**2*(-2*af*e**2/(6*x) - e**2/(6*x**2))
    dT2 = (af*Af/2)*((-af*e*W - af*Af*e**2)/x**2 - 2*e*W/x**3)
    return dT1 + dT2
def VF_f(x):
    e = np.exp(-af*x); W = W0f + Af*e
    return af**2*Af**2*e**2/(6*x) + af*Af*e*W/(2*x**2)
def crit_points(D, n):
    D = float(D)
    dV = lambda x: dVF_f(x) - n*D/x**(n+1)
    xs = np.linspace(60, 600, 20001); ys = dV(xs)
    pts = []
    for i in np.where(np.sign(ys[:-1]) != np.sign(ys[1:]))[0]:
        r = brentq(dV, xs[i], xs[i+1], xtol=1e-12)
        v = VF_f(r) + D/r**n
        is_min = ys[i] < 0 < ys[i+1]
        pts.append((mp.mpf(r), mp.mpf(v), 1 if is_min else -1))
    return pts
def ds_window(n, lo, hi):
    def vmin(D):
        m = [v for (_, v, t) in crit_points(D, n) if t > 0]
        return min(m) if m else None
    # lower edge: bisection on sign of V at local min
    a_, b_ = lo
    for _ in range(60):
        m = (a_+b_)/2
        if vmin(m) is not None and vmin(m) > 0: b_ = m
        else: a_ = m
    Dl = b_
    a_, b_ = hi
    for _ in range(60):
        m = (a_+b_)/2
        if vmin(m) is not None: a_ = m
        else: b_ = m
    return Dl, a_

D = mp.mpf('3e-9')
pts = crit_points(D, 3)
for (s, v, v2) in pts:
    # canonical phi = sqrt(3/2) ln sigma :  d^2V/dphi^2 = (2/3)(s^2 V'' + s V')
    Vfun = lambda x: VF(x) + D/x**3
    v2 = mp.diff(Vfun, s, 2)
    Vpp_can = mp.mpf(2)/3*(s**2*v2 + s*mp.diff(Vfun, s))
    print(f"   D=3e-9 crit pt: sigma={mp.nstr(s,6)} V={mp.nstr(Vfun(s),5)} min? {v2>0}  V''_can/V={mp.nstr(Vpp_can/Vfun(s),5)}")
Dl3, Du3 = ds_window(3, (1e-9, 3.5e-9), (5e-9, 1.2e-8))
print("   dS window D/sigma^3: [", mp.nstr(Dl3, 4), ",", mp.nstr(Du3, 4), "]")
Dl2, Du2 = ds_window(2, (1e-11, 4e-11), (5e-11, 1.5e-10))
print("   dS window D/sigma^2: [", mp.nstr(Dl2, 4), ",", mp.nstr(Du2, 4), "]")

# -------- D2: KPV --------
b02 = mp.mpf('0.93266'); b04 = b02**2
def Veff(psi, pM):
    V2 = mp.sqrt(b04*mp.sin(psi)**4 + (mp.pi*pM - psi + mp.sin(2*psi)/2)**2)/mp.pi
    return V2 - (2*psi - mp.sin(2*psi))/(2*mp.pi)
def has_local_min(pM):
    f = lambda p: mp.diff(lambda x: Veff(x, pM), p)
    xs = [mp.mpf(i)/400 for i in range(1, 1257)]
    prev = f(xs[0])
    for x0, x1 in zip(xs[:-1], xs[1:]):
        cur = f(x1)
        if prev < 0 and cur > 0:
            return mp.findroot(f, (x0+x1)/2)
        prev = cur
    return None
lo, hi = mp.mpf('0.05'), mp.mpf('0.12')
for _ in range(40):
    m = (lo+hi)/2
    if has_local_min(m) is not None: lo = m
    else: hi = m
print("\nD2: KPV critical p/M =", mp.nstr(lo, 5), " ; psi_min(p/M=0.03) =", mp.nstr(has_local_min(mp.mpf('0.03')), 5),
      " ; 2 pi (p/M)/b0^4 =", mp.nstr(2*mp.pi*mp.mpf('0.03')/b04, 5))
I0 = mp.quad(lambda x: (x*mp.coth(x) - 1)/mp.sinh(x)**2*(mp.sinh(2*x) - 2*x)**(mp.mpf(1)/3), [0, 1, 5, mp.inf])
cpp = mp.cbrt(2)/I0
cp = mp.mpf('1.18')
Mmin = mp.mpf(8)/3*mp.sqrt(mp.pi*cp*cpp)
print("    I(0) =", mp.nstr(I0, 6), " c'' = 2^(1/3)/I(0) =", mp.nstr(cpp, 6), " sqrt(gs)M_min = (8/3)sqrt(pi c' c'') =",
      mp.nstr(Mmin, 5), " -> gs M^2 >=", mp.nstr(Mmin**2, 5), " ; Junghans 64 pi c'c''/9 =", mp.nstr(64*mp.pi*cp*cpp/9, 5))

# -------- D3 --------
tab = [(16, 1.051), (16, 0.913), (18, 0.796), (20, 0.808), (16, 0.746)]   # 2406.13751 Table 2: (M, g_s M)
print("\nD3: 2406.13751 Table 2  g_s M^2 =", [round(M*gM, 2) for M, gM in tab], "; p/M =", [round(1/M, 4) for M, _ in tab])
