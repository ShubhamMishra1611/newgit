"""RT2 independent implementation of KPV + alpha'^2 potentials (claim 6).
Normalised potential v = V/(4 pi mu5 M/g_s).
 KPV  : HSV 2208.02826 eq (2.6)  -> v = pi p/M + S - (psi - sin2psi/2)
 HSV  : eq (2.11),(2.12)          -> v_KPV - f(psi) S/(gsM)^2,  f = c1 + 2c2 cot^2 + c2 cot^4, c1=8.825, c2=1.891
 SV   : 2212.07437 eq (3.5) (LaTeX checked on ar5iv), Table 1 constants c1..c8
 S(psi) = sqrt(b0^4 sin^4 psi + (pi p/M - psi + sin(2psi)/2)^2), b0^2 = 2 a0^{1/2}/6^{1/3}
 anti-D3 tension correction (HSV eq 2.15; TASI 2512.17095 eq 6.122): 1 - c/(gsM)^2, c = 5.9241
"""
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq, minimize_scalar
a0 = quad(lambda x: (x/np.tanh(x)-1)/np.sinh(x)**2*(np.sinh(2*x)-2*x)**(1/3), 0, 50)[0]
b0sq = 2*np.sqrt(a0)/6**(1/3); b04 = b0sq**2
c1, c2, c3, c4, c5, c6, c7, c8 = 8.825, 1.891, 2.448, 8.696, 32.61, 2.174, 14.64, 18.65
print('a0=%.5f b0^2=%.5f' % (a0, b0sq))

def S(psi, q):  # q = p/M
    return np.sqrt(b04*np.sin(psi)**4 + (np.pi*q - psi + 0.5*np.sin(2*psi))**2)
def v_kpv(psi, q):
    return np.pi*q + S(psi, q) - (psi - 0.5*np.sin(2*psi))
def fcurv(psi):
    ct = 1/np.tan(psi); return c1 + 2*c2*ct**2 + c2*ct**4
def v_hsv(psi, q, x):  # x = g_s M
    return v_kpv(psi, q) - fcurv(psi)*S(psi, q)/x**2
def sv_bracket(psi, q, x):
    ct = 1/np.tan(psi); sn = np.sin(psi); A = np.pi*q - (psi - 0.5*np.sin(2*psi))
    return 1 + (c3 - c1 + (c4-2*c2)*ct**2 - c2*ct**4 + c5*ct**4/sn**4*A**2 - c6*ct**3/sn**2*A)/x**2
def v_sv(psi, q, x):
    ct = 1/np.tan(psi); sn = np.sin(psi)
    return S(psi, q)*sv_bracket(psi, q, x) + (np.pi*q - (psi - 0.5*np.sin(2*psi)))*(1 + c7/x**2 + c8*ct/(x**2*sn))

PSI = np.linspace(1e-3, np.pi-1e-3, 40001)
def extrema(f, lo=1e-3, hi=np.pi-1e-3):
    v = f(PSI); d = np.sign(np.diff(v)); idx = np.where(np.diff(d) != 0)[0]+1
    mins, maxs = [], []
    for i in idx:
        a, b = PSI[max(i-2, 0)], PSI[min(i+2, len(PSI)-1)]
        if d[i-1] < 0 < d[i]:
            r = minimize_scalar(f, bounds=(a, b), method='bounded', options={'xatol': 1e-12}); mins.append((r.x, r.fun))
        else:
            r = minimize_scalar(lambda t: -f(t), bounds=(a, b), method='bounded', options={'xatol': 1e-12}); maxs.append((r.x, -r.fun))
    return mins, maxs

# --- baseline KPV critical p/M
def has_min_kpv(q):
    m, _ = extrema(lambda t: v_kpv(t, q)); return any(t < np.pi/2 for t, _ in m)
qc = brentq(lambda q: 0.5 - has_min_kpv(q), 0.06, 0.1, xtol=1e-6)
print('KPV critical p/M = %.5f (KPV/HSV: ~0.08)' % qc)

# --- HSV shark fin: smallest g_sM with a metastable minimum below the equator, scanning p/M
def hsv_min(q, x):
    m, M_ = extrema(lambda t: v_hsv(t, q, x)); m = [mm for mm in m if mm[0] < np.pi/2]
    return m
best = None
for q in np.arange(0.070, 0.0861, 0.0005):
    lo, hi = 5.0, 60.0
    if not hsv_min(q, hi): continue
    for _ in range(40):
        mid = 0.5*(lo+hi)
        if hsv_min(q, mid): hi = mid
        else: lo = mid
    if best is None or hi < best[1]: best = (q, hi)
print('HSV shark-fin tip (ours): p/M = %.4f, g_sM = %.3f   [HSV text: 0.0826, 20.363]' % best)
for x in [25., 30., 50., 100.]:
    qs = np.arange(0.02, 0.09, 0.0002); ok = [q for q in qs if hsv_min(q, x)]
    hsvfit = 0.0313227 + 2.73449/x - 72.7164/x**2 + 1182.92/x**3 - 8223.59/x**4
    print('  g_sM=%5.1f: HSV min region p/M in [%.4f, %.4f] (ours); HSV eq.2.13 lower edge %.4f' % (x, min(ok), max(ok), hsvfit))
# HSV Fig.4: max of v positive
def maxpos_thresh(q):
    f = lambda x: max([mv for _, mv in extrema(lambda t: v_hsv(t, q, x))[1]] + [-1e9])
    return brentq(f, 1.0, 20.0, xtol=1e-6)
print('HSV Fig.4 thresholds (max V>0): p/M=1e-4 -> g_sM>%.3f ; p/M=3 -> g_sM>%.3f  [HSV: 5.238, 2.086]' % (maxpos_thresh(1e-4), maxpos_thresh(3.0)))

# --- SV checks: zero-energy line and smallest g_sM with zero-energy metastable min
def sv_meta(q, x):
    m, _ = extrema(lambda t: v_sv(t, q, x)); m = [mm for mm in m if mm[0] < np.pi/2]
    return m
def zero_line(x):
    def g(q):
        m = sv_meta(q, x)
        return m[0][1] if m else np.nan
    qs = np.linspace(0.001, 0.08, 400); vals = np.array([g(q) for q in qs])
    ok = np.where(np.isfinite(vals[:-1]) & np.isfinite(vals[1:]) & (np.sign(vals[:-1]) != np.sign(vals[1:])))[0]
    if len(ok) == 0: return None
    i = ok[0]; return brentq(g, qs[i], qs[i+1], xtol=1e-7)
for x in [3.6, 4., 5., 10., 20.]:
    zl = zero_line(x)
    print('  SV zero-energy line g_sM=%5.1f: p/M=%s ; SV eq.4.1: %.5f' % (x, ('%.5f' % zl) if zl else 'none', 0.1029/x**1.0909))
for x in [3.3, 3.4, 3.5, 3.6, 3.7]:
    qs = np.arange(0.015, 0.035, 0.0002)
    pos = [q for q in qs if sv_meta(q, x) and sv_meta(q, x)[0][1] > 0]
    print('  SV g_sM=%.1f: metastable min with V_up>0 for p/M in %s' % (x, ('[%.4f,%.4f]' % (min(pos), max(pos))) if pos else 'none'))

# --- the five vacua
vac = [('Ex1', 1.0511, 16), ('Ex2', 0.9134, 16), ('Ex3', 0.7961, 18), ('Ex4', 0.8075, 20), ('Ex5', 0.7455, 16)]
print('\nFive vacua (p=1):')
print(' id   g_sM  p/M     psi_KPV  f(psi)/(gsM)^2  c1/(gsM)^2  cot-part share  SVbracket  c_D3/(gsM)^2  R_S3^2/a\'  R_NS5^2/a\'(b0^4 gsM sin^2)  HSV minima<pi/2  HSV maxima  SV minima (psi, v/v_D3)')
for name, x, Mm in vac:
    q = 1/Mm
    r = minimize_scalar(lambda t: v_kpv(t, q), bounds=(0.05, 1.2), method='bounded'); pk = r.x
    fc = fcurv(pk)/x**2; c1p = c1/x**2
    br = sv_bracket(pk, q, x)
    mH, MH = extrema(lambda t: v_hsv(t, q, x)); mS, MS = extrema(lambda t: v_sv(t, q, x))
    vD3 = 2*np.pi*q
    print(' %s %6.3f %6.4f %8.3f %12.1f %12.1f %12.2f %11.1f %11.2f %9.2f %12.3f      %s  %s  %s' % (
        name, x, q, pk, fc, c1p, 1-c1/fcurv(pk), br, 5.9241/x**2, b0sq*x, b04*x*np.sin(pk)**2,
        [round(t, 3) for t, _ in mH if t < np.pi/2], [round(t, 3) for t, _ in MH], [(round(t, 3), round(v/vD3, 1)) for t, v in mS]))
for q in [1/16, 1/20]:
    pk = minimize_scalar(lambda t: v_kpv(t, q), bounds=(0.05, 1.2), method='bounded').x
    F = fcurv(pk)
    print('p/M=%.4f: f(psi_KPV)=%.2f -> f/(gsM)^2<=0.1 needs g_sM>=%.1f ; <=1 needs %.1f ; anti-D3 c/(gsM)^2<=0.1 needs %.1f' % (q, F, np.sqrt(F/0.1), np.sqrt(F), np.sqrt(5.9241/0.1)))
print('anti-D3 tension correction 1 - 5.9241/(gsM)^2 changes sign at g_sM = %.3f' % np.sqrt(5.9241))
