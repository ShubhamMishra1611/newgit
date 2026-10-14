"""T2 task (4): width of the dS window in the uplift coefficient vs. flux discreteness; expected yield vs 5/396.
KKLT single-modulus model (hep-th/0301240 eqs.11-14 and the unnumbered potential of Sec.III B):
  K = -3 ln(T+Tbar), W = W0 + A e^{-aT};  V_F(s) = (aA e^{-as}/(2 s^2)) (s a A e^{-as}/3 + W0 + A e^{-as})
  uplift V_up = D / s^n, n = 2 (anti-D3 in a warped throat, V ~ V_E^{-4/3} ~ s^{-2}; as used in 2406.13751 eq.3.31)
  and n = 3 (original KKLT eq.16).  Xi := V_up(s_AdS)/|V_AdS| (alignment parameter of 2406.13751 eq.3.34).
In terms of x = a s the window [Xi_min, Xi_max] depends only on x_AdS and n.
"""
import json, os
import numpy as np
from scipy.optimize import brentq, minimize_scalar

def window(xA, n):
    a, A = 1.0, 1.0
    W0 = -A*np.exp(-xA)*(1 + 2*xA/3)              # SUSY condition (KKLT eq.13)
    VF = lambda s: a*A*np.exp(-a*s)/(2*s**2)*(s*a*A*np.exp(-a*s)/3 + W0 + A*np.exp(-a*s))
    VAdS = -a**2*A**2*np.exp(-2*xA)/(6*xA)      # (KKLT eq.14) at s = xA/a
    assert abs(VF(xA) - VAdS) < 1e-10*abs(VAdS)
    def minimum(Xi):
        D = Xi*abs(VAdS)*xA**n
        V = lambda s: VF(s) + D/s**n
        ss = np.linspace(0.3*xA, 6*xA, 60001); vv = V(ss); dv = np.diff(vv)
        idx = np.where((dv[:-1] < 0) & (dv[1:] > 0))[0]
        if len(idx) == 0: return None
        i = idx[0]
        r = minimize_scalar(V, bounds=(ss[i], ss[i+2]), method='bounded', options={'xatol': 1e-12})
        return r.x, r.fun
    # Xi_min: V(min) = 0 ; Xi_max: minimum disappears
    f0 = lambda Xi: minimum(Xi)[1]
    lo, hi = 0.5, 30.0
    for _ in range(60):
        m = 0.5*(lo+hi); lo, hi = (m, hi) if minimum(m) is not None else (lo, m)
    Xi_max = lo
    Xi_min = brentq(f0, 0.5, Xi_max*(1-1e-9)) if f0(Xi_max*(1-1e-9)) > 0 else np.nan
    # canonical mass of the modulus at the AdS min in AdS units (phi = sqrt(3/2) ln s)
    h = 1e-4*xA
    d2 = (VF(xA+h) - 2*VF(xA) + VF(xA-h))/h**2
    m2 = (2/3)*xA**2*d2           # d^2V/dphi^2 at a critical point = (2/3) s^2 V''(s)
    mL = np.sqrt(m2*3/abs(VAdS))  # m * l_AdS with l_AdS^2 = 3/|V|
    return Xi_min, Xi_max, mL

print("=== Baseline: KKLT W0=-1e-4, A=1, a=0.1 (x_AdS = 11.3589) ===")
xK = 0.1*113.5891
for n in (3, 2):
    a = window(xK, n)
    print(f"n={n}: Xi window [{a[0]:.3f}, {a[1]:.3f}]  -> D window [{a[0]*1.9966e-15*113.5891**n:.3e}, {a[1]*1.9966e-15*113.5891**n:.3e}]"
          f"  (phase1 grid: n=3 [2.994e-9, 7.713e-9], n=2 [2.637e-11, 9.604e-11])")

print("\n=== Window vs x_AdS = a*sigma (n=2) ===")
tab = []
for xA in (1, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 10, 11.36, 15, 20, 25, 30, 40):
    lo, hi, mL = window(xA, 2); tab.append((xA, lo, hi, mL))
    print(f"x_AdS={xA:6.2f}: Xi in [{lo:.3f}, {hi:.3f}], ln(Xi_max/Xi_min) = {np.log(hi/lo):.3f}, m*l_AdS = {mL:.2f}")
tab = np.array(tab)

# calibrate x_AdS to each example's lightest Kahler mass in AdS units (repo: km mass spectrum in units of 1/lAdS)
here = os.path.dirname(os.path.abspath(__file__))
R = json.load(open(os.path.join(here, '..', 't2_heavy_conifold_results.json')))
D = json.load(open(os.path.join(here, '..', 'data_extracted_data_dS.json')))
mLmin = [min(abs(np.array(D['km mass spectrum in units of 1/lAdS'][i], float))) for i in range(5)]
print("\n=== Per-example window (toy calibrated to lightest Kahler mass in AdS units) and K'-step ===")
print("ID  m_min*l_AdS  x_cal  Xi_min  Xi_max  w=ln(ratio) | Xi_prec(exact) Xi_prec(B.16) in-window? | dlnXi per dK'=1  w/step | c-ratio per dK'")
summ = []
for i, r in enumerate(R):
    g = lambda x: window(x, 2)[2] - mLmin[i]
    xc = brentq(g, 0.6, 60)
    lo, hi, _ = window(xc, 2); w = np.log(hi/lo)
    step = 4*np.pi/3/(r['gD']*r['M'])        # Xi ∝ |z|^{4/3}, |z| ∝ exp(-pi K'/(g_s M)) (n_cf = 2, eq.3.11)
    inw = lo <= r['Xi_exact'] <= hi
    summ.append((w, step))
    print(f"Ex{r['ID']}  {mLmin[i]:8.2f}  {xc:6.2f}  {lo:.3f}  {hi:.3f}  {w:.3f}      | {r['Xi_exact']:.2f}  {r['Xi_B16']:.2f}  {'yes' if inw else 'NO '}"
          f"   | {step:.2f}  {w/step:.3f}  | {np.exp(-step):.2e}")

print("\n=== Expected yield vs reported numbers (2406.13751 Table 1, Fig.3 caption, Sec.5) ===")
N_antiD3, N_aligned, N_dS_flux, N_dS_vac = 33371, 396, 5, 30
print(f"aligned fraction of anti-D3 PFVs: {N_aligned}/{N_antiD3} = {N_aligned/N_antiD3:.4f}")
for label, w in (("toy n=2, x=11.36", np.log(window(11.36, 2)[1]/window(11.36, 2)[0])),
                 ("toy n=2, calibrated range", np.mean([s[0] for s in summ])),
                 ("paper's stated 1<~Xi<~3", np.log(3.0)),
                 ("empirical precursor range of the 5 (B.16: 1.52-3.93)", np.log(3.927/1.521))):
    P = w/np.log(100.0)
    print(f"{label:55s}: w = {w:.3f} -> P(in window | log-uniform on [0.1,10]) = {P:.3f}; expected ~{P*N_aligned:.0f} of 396; observed {N_dS_flux} -> efficiency {N_dS_flux/(P*N_aligned):.3f}")
print(f"observed per-flux-choice yield 5/396 = {5/396:.4f}; 30 vacua / 396 = {30/396:.4f}")
print(f"density of aligned points in ln Xi (if uniform): {N_aligned/np.log(100):.0f} per e-fold -> mean spacing {np.log(100)/N_aligned:.4f}")
