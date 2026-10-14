"""Quick sanity checks for Front D (de Sitter vacua). Phase 1 only.
1. KKLT toy model (hep-th/0301240, eqs. (11)-(14) and unnumbered potential in Sec. III B).
2. KPV NS5 potential (hep-th/0112197, eqs. (9), (40), (45)): critical p/M.
3. Conifold destabilization bound (1809.06861 eqs. (2.10),(2.18),(2.33),(3.3); 2201.03572 eq. (2.9)).
"""
import numpy as np
from scipy.optimize import brentq, minimize_scalar
from scipy.integrate import quad

print("=== 1. KKLT toy model ===")
W0, A, a = -1e-4, 1.0, 0.1
# DW=0 condition (KKLT eq. 13): W0 = -A e^{-a s}(1 + 2 a s/3)
f = lambda s: W0 + A*np.exp(-a*s)*(1 + 2*a*s/3)
s_cr = brentq(f, 10, 500)
VAdS = -a**2*A**2*np.exp(-2*a*s_cr)/(6*s_cr)  # eq. (14)
print(f"sigma_cr = {s_cr:.4f}  (KKLT text: ~113)")
print(f"V_AdS    = {VAdS:.4e}  (eq. 14)")
def VF(s):  # unnumbered eq. in KKLT Sec. III B, without uplift
    e = np.exp(-a*s)
    return a*A*e/(2*s**2)*(s*a*A*e/3 + W0 + A*e)
print(f"check V_F(sigma_cr) = {VF(s_cr):.4e}")
def scan(D, n):
    V = lambda s: VF(s) + D/s**n
    ss = np.linspace(50, 2000, 200001)
    vv = V(ss)
    dv = np.diff(vv)
    idx_min = np.where((dv[:-1] < 0) & (dv[1:] > 0))[0]
    idx_max = np.where((dv[:-1] > 0) & (dv[1:] < 0))[0]
    out = []
    for i in idx_min:
        r = minimize_scalar(V, bracket=(ss[i-5], ss[i+1], ss[i+5]))
        out.append(("min", r.x, r.fun))
    for i in idx_max:
        r = minimize_scalar(lambda s: -V(s), bracket=(ss[i-5], ss[i+1], ss[i+5]))
        out.append(("max", r.x, -r.fun))
    return V, out
# canonical field phi = sqrt(3/2) ln sigma (KKLT Sec. IV): dV/dphi = sigma V'(sigma) sqrt(2/3)
def swamp(V, s, h=1e-3):
    d1 = (V(s+h)-V(s-h))/(2*h); d2 = (V(s+h)-2*V(s)+V(s-h))/h**2
    # phi = sqrt(3/2) ln s  =>  d/dphi = sqrt(2/3) s d/ds
    Vp = np.sqrt(2/3)*s*d1
    Vpp = (2/3)*(s**2*d2 + s*d1)
    return Vp/V(s), Vpp/V(s)
for (D, n) in [(3e-9, 3), (0.0, 3)]:
    V, ext = scan(D, n)
    print(f"D={D}, uplift D/sigma^{n}:")
    for kind, s, v in ext:
        g, h = swamp(V, s)
        print(f"   {kind}: sigma={s:.3f}  V={v:.4e}  |V'|/V={abs(g):.2e}  V''/V={h:.3e} (canonical, Planck units)")
# range of D giving a dS minimum (n=3), and the analogous range for n=2
for n in (3, 2):
    Ds = np.logspace(-11, -7, 400)
    ok = []
    for D in Ds:
        V, ext = scan(D, n)
        mins = [e for e in ext if e[0] == "min"]
        if mins and mins[0][2] > 0: ok.append(D)
    if ok: print(f"n={n}: dS minimum for D in [{min(ok):.3e}, {max(ok):.3e}] (grid)")

print("=== 2. KPV NS5 potential ===")
b0sq = 0.93266  # KPV eq. (9)
b04 = b0sq**2
def Veff(psi, r):  # KPV eqs. (40),(45) in units of A0; r = p/M
    V2 = np.sqrt(b04*np.sin(psi)**4 + (np.pi*r - psi + 0.5*np.sin(2*psi))**2)/np.pi
    return V2 - (2*psi - np.sin(2*psi))/(2*np.pi)
def has_min(r):
    ps = np.linspace(1e-4, 3.0, 300001)  # avoid round-off noise near psi=pi
    dv = np.diff(Veff(ps, r))
    return np.any((dv[:-1] < 0) & (dv[1:] > 0))
lo, hi = 0.01, 0.2
assert has_min(0.03) and not has_min(0.15)
for _ in range(50):
    mid = 0.5*(lo+hi)
    if has_min(mid): lo = mid
    else: hi = mid
print(f"critical p/M = {lo:.5f}  (KPV text: ~8%)")
r = 0.03
ps = np.linspace(1e-4, 3.0, 300001); v = Veff(ps, r); dv = np.diff(v)
i = np.where((dv[:-1] < 0) & (dv[1:] > 0))[0][0]
print(f"p/M=0.03: psi_min={ps[i+1]:.4f}, small-psi estimate 2 pi p/(b0^4 M)={2*np.pi*r/b04:.4f} (eq. 46)")
print(f"   V(psi_min)/A0 = {v[i+1]:.6f}, eq.(47) estimate p/M(1-8pi^2 p^2/(3 b0^12 M^2)) = {r*(1-8*np.pi**2*r**2/(3*b0sq**6)):.6f}")

print("=== 3. Conifold destabilization bound ===")
I = lambda T: quad(lambda x: (x/np.tanh(x)-1)/np.sinh(x)**2*(np.sinh(2*x)-2*x)**(1/3), T, 60.0, limit=400)[0]  # integrand ~ x e^{-4x/3}; tail beyond 60 negligible
I0 = I(1e-6)
cpp = 2**(1/3)/I0
cp = 1.18  # 1809.06861 eq. (2.18), taken from the literature, not recomputed here
Mmin = 8/3*np.sqrt(np.pi*cp*cpp)
print(f"I(0) = {I0:.5f};  c'' = 2^(1/3)/I(0) = {cpp:.4f}  (1809.06861 eq. 2.33: ~1.75)")
print(f"sqrt(g_s) M_min = (8/3) sqrt(pi c' c'') = {Mmin:.3f}  (eq. 3.3: ~6.8)")
print(f"g_s M^2 >= 64 pi c' c''/9 = {64*np.pi*cp*cpp/9:.2f}  (2201.03572 eq. 2.9: ~46.1)")

print("=== 4. Control numbers for the five candidate dS vacua (2505.00149 Table 1 / 2406.13751) ===")
# (ID, M, g_s, g_s M as printed); p = 1 anti-D3 in all examples
table = [(1, 16, 0.0657, 1.051), (2, 16, 0.0571, 0.913), (3, 18, 0.0442, 0.796),
         (4, 20, 0.0404, 0.808), (5, 16, 0.0466, 0.746)]
for ID, M, gs, gsM in table:
    gsM2 = gs*M**2
    rate_exp = 27*b0sq**6*gs*M**6/(512*np.pi*1**3)   # KPV eq. (65), p = 1
    print(f"Ex{ID}: g_sM={gs*M:.3f} (table {gsM}), g_sM^2={gsM2:.1f}, p/M={1/M:.4f}, "
          f"g_sM^2/46.1={gsM2/46.1:.2f}, KPV decay exponent={rate_exp:.3g}")
print("=== 5. Discreteness of warped uplift vs width of KKLT-toy dS window ===")
for gsM in (0.8, 1.0, 2.0, 5.0):
    print(f"g_sM={gsM}: K->K+1 rescales e^(4A_min) ~ exp(-8 pi K/(3 g_s M)) by {np.exp(-8*np.pi/(3*gsM)):.3g}; "
          f"KKLT-toy window D_max/D_min (n=3) ~ 2.58")
