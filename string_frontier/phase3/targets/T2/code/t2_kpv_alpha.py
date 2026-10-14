"""T2 task (3): anti-D3/NS5 (KPV) metastability with known alpha'^2 corrections at g_s M ~ 1.
Potentials, all in units of 4 pi mu5 M / g_s (overall 1/g_s^2 factor does not affect shape):
  F(psi) = pi p/M - psi + sin(2psi)/2,  S(psi) = sqrt(b0^4 sin^4 psi + F^2),  b0^2 = 0.93266
  KPV  (hep-th/0112197 eqs.40,45; HSV 2208.02826 eq.2.6):  v_KPV = S + F
  HSV  (2208.02826 eqs.2.11-2.12, curvature only):        v_HSV = v_KPV - S (c1 + 2c2 cot^2 + c2 cot^4)/(g_sM)^2
  SV   (2212.07437 eq.3.5 with A.16-A.23, all known alpha'^2 terms):
        v_SV = S [1 + (c3-c1 + (c4-2c2) cot^2 - c2 cot^4 + c5 cot^4 F^2/sin^4 - c6 cot^3 F/sin^2)/(g_sM)^2]
             + F [1 + c7/(g_sM)^2 + c8 cot/((g_sM)^2 sin)]
  (Our reading of the text-extracted eq.3.5: the combination written 'pi p/M - psi - sin(2psi)/2' is
   pi p/M - (psi - sin(2psi)/2) = F, and the c7,c8 bracket multiplies the full F, as in eqs. A.20, A.22.
   This reading is validated below against SV's quoted zero-energy line, eq. 4.1.)
"""
import numpy as np
from scipy.optimize import brentq, minimize_scalar

b0sq = 0.93266; b04 = b0sq**2
c1, c2, c3, c4, c5, c6, c7, c8 = 8.825, 1.891, 2.448, 8.696, 32.61, 2.174, 14.64, 18.65

def parts(psi, r):
    F = np.pi*r - psi + 0.5*np.sin(2*psi)
    S = np.sqrt(b04*np.sin(psi)**4 + F**2)
    return F, S
def v_kpv(psi, r, g=None):
    F, S = parts(psi, r); return S + F
def v_hsv(psi, r, g):
    F, S = parts(psi, r); ct = 1/np.tan(psi)
    return S + F - S*(c1 + 2*c2*ct**2 + c2*ct**4)/g**2
def v_sv(psi, r, g):
    F, S = parts(psi, r); ct = 1/np.tan(psi); sn = np.sin(psi)
    br = (c3 - c1 + (c4 - 2*c2)*ct**2 - c2*ct**4 + c5*ct**4*F**2/sn**4 - c6*ct**3*F/sn**2)/g**2
    return S*(1 + br) + F*(1 + c7/g**2 + c8*ct/(g**2*sn))

def local_minima(fun, r, g, lo=0.02, hi=np.pi/2+0.3, n=40001):
    ps = np.linspace(lo, hi, n); v = fun(ps, r, g); dv = np.diff(v)
    idx = np.where((dv[:-1] < 0) & (dv[1:] > 0))[0]
    out = []
    for i in idx:
        res = minimize_scalar(lambda x: fun(x, r, g), bounds=(ps[i], ps[i+2]), method='bounded', options={'xatol': 1e-10})
        out.append((res.x, res.fun))
    return out

print("=== Baseline: leading-order KPV ===")
def kpv_has_min(r): return len(local_minima(v_kpv, r, None, 1e-3, 3.0, 300001)) > 0
lo, hi = 0.01, 0.2
for _ in range(45):
    m = 0.5*(lo+hi); lo, hi = (m, hi) if kpv_has_min(m) else (lo, m)
print(f"(p/M)_crit(KPV) = {lo:.5f}")

print("\n=== Validation 1: HSV curvature-corrected potential (2208.02826 Fig.3 / eq.2.13) ===")
def hsv_min_exists(r, g):
    return any(x > 0.05 for x, _ in local_minima(v_hsv, r, g, 0.02, 2.5, 20001))
# tip of 'shark fin': minimal g_s M (over p/M) for which an interior minimum exists.
gs_grid = np.concatenate([np.arange(10, 40, 0.25), np.arange(40, 201, 5.0)])
best = (None, 1e9)
for r in np.linspace(0.074, 0.090, 65):
    okg = [g for g in gs_grid if hsv_min_exists(r, g)]
    if okg and min(okg) < best[1]: best = (r, min(okg))
r0 = best[0]
glo, ghi = best[1]-0.25, best[1]
for _ in range(30):
    gm = 0.5*(glo+ghi); glo, ghi = (glo, gm) if hsv_min_exists(r0, gm) else (gm, ghi)
print(f"shark-fin tip (ours, grid dp/M=0.00025): p/M = {r0:.4f}, g_sM = {ghi:.3f}   [HSV text: p/M=0.0826, g_sM=20.363]")
for g in (25.0, 30.0, 50.0, 100.0):
    fit = 0.0313227 + 2.73449/g - 72.7164/g**2 + 1182.92/g**3 - 8223.59/g**4
    rlo, rhi = 0.0, 0.0826
    # lower boundary: smallest r with a minimum at this g
    rs = np.linspace(0.02, 0.085, 651); ok = [r for r in rs if hsv_min_exists(r, g)]
    print(f"g_sM={g:5.1f}: lower edge of min-region (ours) p/M = {min(ok):.4f}; HSV fit eq.2.13 = {fit:.4f}")

# HSV Fig.4: boundary g_sM above which the maximum of V_HSV is positive (asymptotes 5.238 / 2.086)
def maxval(r, g):
    ps = np.linspace(0.01, np.pi-0.01, 200001); v = v_hsv(ps, r, g); dv = np.diff(v)
    mx = ps[1:-1][(dv[:-1] > 0) & (dv[1:] < 0)]
    return max(v_hsv(m, r, g) for m in mx) if len(mx) else np.nan
for r in (1e-4, 3.0):
    gs_ = np.linspace(1.5, 8, 131); vals = [maxval(r, g) for g in gs_]
    root = next((brentq(lambda g: maxval(r, g), gs_[k], gs_[k+1]) for k in range(len(gs_)-1)
                 if np.isfinite(vals[k]) and np.isfinite(vals[k+1]) and vals[k] < 0 <= vals[k+1]), None)
    print(f"HSV Fig.4 check: p/M={r}: max V>0 for g_sM > {root:.3f}  [HSV text: 5.238 (p/M<<0.1), 2.086 (p/M>>0.1)]")

print("\n=== Validation 2: SV full alpha'^2 potential, zero-energy line (2212.07437 eq.4.1) ===")
def sv_min(r, g):
    m = local_minima(v_sv, r, g, 0.005, 1.6, 40001)
    return m[0] if m else None
for g in (4.0, 5.0, 10.0, 20.0):
    # find r where V_up = v_sv(psi_min) changes sign
    f = lambda r: (sv_min(r, g) or (np.nan, np.nan))[1]
    rs = np.linspace(0.002, 0.08, 157); vals = [f(r) for r in rs]
    root = None
    for k in range(len(rs)-1):
        if np.isfinite(vals[k]) and np.isfinite(vals[k+1]) and vals[k]*vals[k+1] < 0:
            root = brentq(f, rs[k], rs[k+1]); break
    print(f"g_sM={g:5.1f}: zero-energy p/M (ours) = {root if root is None else round(root,5)};  eq.4.1: {0.1029/g**1.0909:.5f}")

# SV threshold: smallest g_sM with a zero-energy metastable vacuum (SV text: g_sM ~ 3.6 at p/M ~ 0.025)
for g in (3.4, 3.6, 3.8):
    sg = [(r, m[1]) for r in np.linspace(0.005, 0.05, 181) for m in [sv_min(r, g)] if m is not None]
    pos = [r for r, v in sg if v > 0]
    print(f"SV check g_sM={g}: V_up>0 minima exist for p/M in " + (f"[{min(pos):.4f},{max(pos):.4f}]" if pos else "none"))

print("\n=== Candidate parameters: p = 1, p/M = 1/16..1/20, g_sM in [0.7, 1.1] ===")
print("columns: g_sM, p/M | KPV psi_min, v_KPV/(2pi p/M) | at psi_KPV: delta_curv, SV-DBI bracket, SV-CS bracket,"
      " (v_SV-v_KPV)/v_KPV | R_S3^2/alpha'=b0^2 g_sM, Rt^2=g_sM sin^2 psi | HSV interior min? | SV min: psi, Rt^2, V_up/(2pi p/M)")
rows = []
for g in (0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00, 1.05, 1.10):
    for M in (16, 18, 20):
        r = 1/M
        pk = local_minima(v_kpv, r, None, 1e-3, 3.0, 300001)[0]
        psi = pk[0]; F, S = parts(psi, r); ct = 1/np.tan(psi); sn = np.sin(psi)
        dcurv = (c1 + 2*c2*ct**2 + c2*ct**4)/g**2
        dsv = (c3 - c1 + (c4 - 2*c2)*ct**2 - c2*ct**4 + c5*ct**4*F**2/sn**4 - c6*ct**3*F/sn**2)/g**2
        dcs = c7/g**2 + c8*ct/(g**2*sn)
        rel = (v_sv(psi, r, g) - v_kpv(psi, r))/v_kpv(psi, r)
        hs = [m for m in local_minima(v_hsv, r, g, 0.02, 2.5, 20001)]
        sm = sv_min(r, g)
        line = (f"{g:4.2f} 1/{M:2d} | {psi:.3f} {pk[1]/(2*np.pi*r):.3f} | {dcurv:7.1f} {dsv:8.1f} {dcs:7.1f} {rel:8.1f} |"
                f" {b0sq*g:.2f} {g*np.sin(psi)**2:.3f} | {'yes' if hs else 'no '} | ")
        line += (f"{sm[0]:.3f} {g*np.sin(sm[0])**2:.4f} {sm[1]/(2*np.pi*r):+.3f}" if sm else "none")
        print(line)
        rows.append((g, M, psi, dcurv, dsv, dcs, rel, sm))

print("\n=== The five vacua (g_sM, M from 2406.13751 Table 2) ===")
for ID, g, M in [(1, 1.051, 16), (2, 0.913, 16), (3, 0.796, 18), (4, 0.808, 20), (5, 0.746, 16)]:
    r = 1/M; psi = local_minima(v_kpv, r, None, 1e-3, 3.0, 300001)[0][0]
    F, S = parts(psi, r); ct = 1/np.tan(psi)
    dcurv = (c1 + 2*c2*ct**2 + c2*ct**4)/g**2
    sm = sv_min(r, g); hs = local_minima(v_hsv, r, g, 0.02, 2.5, 20001)
    print(f"Ex{ID}: p/M={r:.4f} g_sM={g}: psi_KPV={psi:.3f}, R_NS5^2/alpha' = b0^4 g_sM sin^2psi = {b04*g*np.sin(psi)**2:.3f},"
          f" curvature corr. factor = {dcurv:.1f}, constant term c1/(g_sM)^2 alone = {c1/g**2:.1f}; HSV interior min: {'yes' if hs else 'no'};"
          f" SV min: " + (f"psi={sm[0]:.3f}, V_up/V_KPV-antiD3 = {sm[1]/(2*np.pi*r):+.3f}" if sm else "none"))
# minimal g_sM at which the alpha'^2 corrections at the KPV minimum are <= 10% (curvature factor), p/M=1/16
for M in (16, 20):
    r = 1/M; psi = local_minima(v_kpv, r, None, 1e-3, 3.0, 300001)[0][0]; ct = 1/np.tan(psi)
    num = c1 + 2*c2*ct**2 + c2*ct**4
    print(f"p/M=1/{M}: curvature factor at psi_KPV = {num:.1f}/(g_sM)^2 -> <=0.1 needs g_sM >= {np.sqrt(num/0.1):.1f}; <=1 needs g_sM >= {np.sqrt(num):.1f}")

print("\n=== Full-range scan psi in (0.001, pi-0.001) of the truncated potentials at the five vacua ===")
for ID, g, M in [(1, 1.051, 16), (2, 0.913, 16), (3, 0.796, 18), (4, 0.808, 20), (5, 0.746, 16)]:
    r = 1/M; ps = np.linspace(1e-3, np.pi-1e-3, 400001)
    out = []
    for name, fun in (("HSV", v_hsv), ("SV", v_sv)):
        v = fun(ps, r, g); dv = np.diff(v)
        mins = ps[1:-1][(dv[:-1] < 0) & (dv[1:] > 0)]; maxs = ps[1:-1][(dv[:-1] > 0) & (dv[1:] < 0)]
        out.append(f"{name}: minima psi={np.round(mins,3).tolist()} (V/V_antiD3={[round(float(fun(m,r,g))/(2*np.pi*r),1) for m in mins]}),"
                   f" maxima psi={np.round(maxs,3).tolist()}")
    print(f"Ex{ID} (g_sM={g}, p/M=1/{M}): " + " | ".join(out))
print("Interpretation: psi>pi/2 is past the equator (anti-D3 charge annihilated; SV: F2 terms out of control there).")
