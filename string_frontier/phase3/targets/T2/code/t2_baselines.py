"""T2 baselines (re-derived here; cf. phase1/D_dS_checks/quick_checks.py).
1. KKLT hep-th/0301240 eqs.(13),(14): sigma_cr, V_AdS for W0=-1e-4, A=1, a=0.1.
2. KPV hep-th/0112197 eqs.(40),(45), b0^2=0.93266 (eq.9): critical p/M.
3. BDGL 1809.06861 eqs.(2.10),(2.18),(2.33),(3.3); Junghans 2201.03572 eq.(2.9).
"""
import numpy as np
from scipy.optimize import brentq
from scipy.integrate import quad

print("=== 1. KKLT ===")
W0, A, a = -1e-4, 1.0, 0.1
s_cr = brentq(lambda s: W0 + A*np.exp(-a*s)*(1 + 2*a*s/3), 10, 500)
VAdS = -a**2*A**2*np.exp(-2*a*s_cr)/(6*s_cr)
print(f"sigma_cr = {s_cr:.4f}   (target 113.589)")
print(f"V_AdS    = {VAdS:.4e} (target -1.997e-15)")

print("=== 2. KPV ===")
b04 = 0.93266**2
def V(psi, r):
    return np.sqrt(b04*np.sin(psi)**4 + (np.pi*r - psi + 0.5*np.sin(2*psi))**2)/np.pi - (2*psi - np.sin(2*psi))/(2*np.pi)
def has_min(r):
    ps = np.linspace(1e-4, 3.0, 300001); dv = np.diff(V(ps, r))
    return np.any((dv[:-1] < 0) & (dv[1:] > 0))
lo, hi = 0.01, 0.2
for _ in range(50):
    mid = 0.5*(lo+hi); lo, hi = (mid, hi) if has_min(mid) else (lo, mid)
print(f"(p/M)_crit = {lo:.5f}   (target ~0.0805; KPV text 'about 8%')")

print("=== 3. BDGL / Junghans ===")
f = lambda x: (x/np.tanh(x) - 1)/np.sinh(x)**2*(np.sinh(2*x) - 2*x)**(1/3)
I0 = quad(f, 1e-8, 60, limit=400)[0]
cpp = 2**(1/3)/I0; cp = 1.18
print(f"I(0) = {I0:.5f}; c'' = 2^(1/3)/I(0) = {cpp:.4f} (BDGL eq.2.33: ~1.75)")
print(f"Sigma = 2^(2/3) I(0) = {2**(2/3)*I0:.5f} (2406.13751 eq.B.8: 1.13983); 2/Sigma = {2/(2**(2/3)*I0):.4f} = c''")
print(f"sqrt(g_s) M_min = (8/3) sqrt(pi c' c'') = {8/3*np.sqrt(np.pi*cp*cpp):.3f} (eq.3.3: 6.8)")
print(f"g_s M^2 >= 64 pi c' c''/9 = {64*np.pi*cp*cpp/9:.2f} with c''={cpp:.4f};  = {64*np.pi*cp*1.75/9:.2f} with c''=1.75 (2201.03572 eq.2.9: 46.1)")
