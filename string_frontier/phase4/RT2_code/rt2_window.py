"""RT2 independent check of claim 7 (dS window vs flux steps, expected yield).
KKLT single modulus (hep-th/0301240): K=-3ln(T+Tbar), W=W0+A e^{-aT}; F-term potential along real sigma
  V_F = a A e^{-a s}/(2 s^2) * (s a A e^{-a s}/3 + W0 + A e^{-a s})   [KKLT eq.(12)-(13) form]
Uplift D/s^n. Xi := D/s_AdS^n / |V_AdS| (uplift at the AdS point over |V_AdS|).
Window: Xi_min (uplifted local min reaches V=0) .. Xi_max (local min disappears).
"""
import numpy as np, json
from scipy.optimize import brentq, minimize_scalar
def VF(s, W0, A, a): e = A*np.exp(-a*s); return a*e/(2*s**2)*(s*a*e/3 + W0 + e)
def sAdS(W0, A, a): return brentq(lambda s: W0 + A*np.exp(-a*s)*(1+2*a*s/3), 1e-3, 1e4)
def window(x, n, a=0.1, A=1.0):
    s0 = x/a; W0 = -A*np.exp(-x)*(1+2*x/3); VA = VF(s0, W0, A, a)
    def localmin(D):
        f = lambda s: VF(s, W0, A, a) + D/s**n
        ss = np.linspace(0.3*s0, 6*s0, 20001); v = f(ss); d = np.sign(np.diff(v))
        idx = np.where((d[:-1] < 0) & (d[1:] > 0))[0]
        if len(idx) == 0: return None
        r = minimize_scalar(f, bounds=(ss[idx[0]-1], ss[idx[0]+2]), method='bounded', options={'xatol': 1e-10*s0}); return r.fun
    Xi = lambda D: D/s0**n/abs(VA)
    Dsc = abs(VA)*s0**n
    Dmax = brentq(lambda t: 1.0 if localmin(t*Dsc) is not None else -1.0, 1.0, 20.0, xtol=1e-7)*Dsc
    Dmin = brentq(lambda t: localmin(t*Dsc), 0.5, Dmax/Dsc*0.999999, xtol=1e-9)*Dsc
    # mass of sigma at SUSY AdS in units of 1/l_AdS: m^2 = V''/(2 K_TTbar), K_TTbar=3/(4 s^2); l_AdS^2 = 3/|V_AdS|
    h = 1e-4*s0; Vpp = (VF(s0+h, W0, A, a) - 2*VA + VF(s0-h, W0, A, a))/h**2
    m2 = Vpp/(2*3/(4*s0**2)); ml = np.sqrt(m2*3/abs(VA))
    return Xi(Dmin), Xi(Dmax), ml, s0, W0, VA
s0 = sAdS(-1e-4, 1.0, 0.1); x0 = 0.1*s0
print('KKLT baseline: sigma_AdS = %.4f (KKLT: 113.589), x=%.4f' % (s0, x0))
for n in [2, 3]:
    lo, hi, ml, s, W0, VA = window(x0, n)
    print(' n=%d: Xi window [%.3f, %.3f]  D window [%.4g, %.4g]  m*l_AdS=%.2f  V_AdS=%.4g' % (n, lo, hi, lo*abs(VA)*s**n, hi*abs(VA)*s**n, ml, VA))
xs = [1, 2, 2.5, 4, 8, 11.36, 20, 30]
for x in xs:
    lo, hi, ml, *_ = window(x, 2); print('  x=%5.2f: Xi in [%.3f, %.3f] w=%.3f m*l_AdS=%.2f' % (x, lo, hi, np.log(hi/lo), ml))
J = json.load(open('/home/user/newgit/string_frontier/phase3/targets/T2/data_extracted_data_dS.json'))
mlad = np.array([min(r) for r in J['km mass spectrum in units of 1/lAdS']])
gsM = np.array(J['gsM'])
print('\nper-vacuum calibration (match lightest Kahler mass * l_AdS):')
ws = []
for i in range(5):
    xc = brentq(lambda x: window(x, 2)[2] - mlad[i], 0.5, 40)
    lo, hi, ml, *_ = window(xc, 2); w = np.log(hi/lo); ws.append(w)
    step = (4*np.pi/3)/gsM[i]
    print(' Ex%d: m l_AdS=%.2f -> x_cal=%.3f Xi window [%.3f, %.3f] factor %.2f w=%.3f | step dlnXi per dK\'=1: %.3f (factor %.0f) | w/step=%.3f' % (i+1, mlad[i], xc, lo, hi, hi/lo, w, step, np.exp(step), w/step))
ws = np.array(ws)
print('\nexpected yields under log-uniform Xi on [0.1,10] (assumption):')
for lab, w in [('mean calibrated w', ws.mean()), ('min w', ws.min()), ('max w', ws.max()), ('x=11.36 w', np.log(window(11.36, 2)[1]/window(11.36, 2)[0])), ('paper 1<~Xi<~3', np.log(3))]:
    P = w/np.log(100); print('  %-20s w=%.3f P=%.3f expected of 396: %.0f ; 5/expected=%.3f' % (lab, w, P, 396*P, 5/(396*P)))

# --- precursor alignment Xi of the five vacua (2406.13751 (B.16) and exact from repo V_AdS)
cabs = lambda v: abs(complex(v['re'], v['im']))
gsS = np.array(J['SUSY gs']); zS = np.array([cabs(v) for v in J['SUSY zcf']]); W0S = np.array([cabs(v) for v in J['SUSY W0']])
VtS = np.array(J['SUSY volTilde']); VEA = np.array(J['corrected CY volume AdS']); VAdS = np.array(J['V AdS'])
M = np.array([J['M vector'][i][0] for i in range(5)], float)
zeta_c = 114.0368; eta = 3*zeta_c/128
XiB16 = zS**(4/3)/W0S**2*VEA**(2/3)*VtS**(1/3)/(gsS*M)**2*zeta_c
Xiex = eta*zS**(4/3)/(gsS*M**2*VtS**(2/3))/VEA**(4/3)/abs(VAdS)
print('\nXi precursor (B.16):', np.round(XiB16, 3)); print('Xi precursor exact (repo V_AdS):', np.round(Xiex, 3))
