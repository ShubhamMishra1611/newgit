import sys, time
import numpy as np
sys.path.insert(0, "/home/user/newgit/string_frontier/round2/review/RR2_code")
from rr2_casimir_indep import *

np.set_printoptions(precision=6)
print("=" * 70); print("A. Casimir density checks")
for (L, m) in [(10., 0.01), (100., 0.005), (300., 0.002), (50., 0.05), (1000., 1e-4)]:
    b = rho_bessel(L, m); t = rho_thermal(L, m)
    print(f"L={L:7.1f} m={m:.1e}  Bessel={b:.10e} thermal={t:.10e} rel.diff={b/t-1:.2e}")
print("massless limit: rho_bessel(L,1e-9)/rho0 =", rho_bessel(100., 1e-9 / 100., kmax=200000) / rho0(100.))
# kmax truncation
for kmax in (50, 200, 2000, 20000):
    print(f"kmax={kmax:6d}: rho(L=100,m=1e-4)/rho0 = {rho_bessel(100., 1e-4, kmax)/rho0(100.):.12f}")

print("\nIMV eq.3.11 small-mR series vs exact (and with the m^3 term restored)")
for mR in (0.003, 0.01, 0.03, 0.1, 0.3):
    R = 1.0; m = mR / R; L = 2 * PI * R
    exact = rho_bessel(L, m)
    imv = PI**2 / L**4 * (1 / 90 - mR**2 / 6 + mR**4 / 48)
    hiT = PI**2 / (90 * L**4) - m**2 / (24 * L**2) + m**3 / (12 * PI * L)  # high-T expansion
    # next term of high-T expansion: + m^4/(64 pi^2) [ln(m^2 L^2) - ln a_b], ln a_b = ln(16 pi^2) + 3/2 - 2 gamma_E
    lnab = np.log(16 * PI**2) + 1.5 - 2 * np.euler_gamma
    hiT4 = hiT + m**4 / (64 * PI**2) * (np.log(m**2 * L**2) - lnab)
    print(f"mR={mR:5.3f}: IMV3.11 rel.err={imv/exact-1:+.3e}  (m^3 term)/exact={m**3/(12*PI*L)/exact:+.3e}"
          f"  with m^3: {hiT/exact-1:+.3e}  with m^3+m^4log: {hiT4/exact-1:+.3e}")

print("=" * 70); print("B. bosons + Lambda only: R_max")
lam = 2.6e-11
Ls = np.geomspace(10, 2000, 200001)
V = (lam - 4 * rho0(Ls)) / Ls**2
Rmax_num = Ls[np.argmax(V)] / (2 * PI)
print(f"numerical R_max = {Rmax_num:.4f} eV^-1 ; (120 pi^2 Lambda)^(-1/4) = {(120*PI**2*lam)**-0.25:.4f} eV^-1"
      f" = {(120*PI**2*lam)**-0.25*1e9:.4e} GeV^-1 ; 1/(2pi R_max) = {1/(2*PI*(120*PI**2*lam)**-0.25)*1e3:.4f} meV")

print("=" * 70); print("C. IMV reproduction, Lambda = 2.6e-11 eV^4, PDG-2016 splittings (assumed)")
d21p = 7.53e-5; d32NO = 2.44e-3; d32IO = 2.51e-3
d31NO = d32NO + d21p          # PDG 2016 quotes Dm32 for NO
t0 = time.time()
for order, d3l in (("NO", d31NO), ("IO", -d32IO)):
    a = onset(lam, d21p, d3l, order, 4, "AdS"); s = onset(lam, d21p, d3l, order, 4, "any")
    print(f"Dirac {order}: dS onset (first min) = {s*1e3:.4f} meV ; AdS onset = {a*1e3:.4f} meV")
    for nf in (2,):
        r, Lb, nl = rhomin_AdS(0.0, d21p, d3l, order, nf)
        print(f"   Majorana {order} m_l=0: rho_min = {r:.4e} eV^4 = {r/lam:.3f} x Lambda(2.6e-11) ; L*={Lb:.2f}")
# variant: if one (mis)reads 2.44e-3 as Dm31 for NO
a = onset(lam, d21p, d32NO, "NO", 4, "AdS"); s = onset(lam, d21p, d32NO, "NO", 4, "any")
print(f"  [variant NO with Dm31=2.44e-3] dS onset {s*1e3:.4f} meV, AdS onset {a*1e3:.4f} meV")
print(f"  (time {time.time()-t0:.1f}s)")

print("=" * 70); print("D. Planck 2018 rho_Lambda")
G = 6.67430e-11; hbar = 1.054571817e-34; c = 299792458.0; eV = 1.602176634e-19
Mpc = 3.0856775814913673e22
MP = np.sqrt(hbar * c / (8 * PI * G)) * c**2 / eV         # reduced Planck mass in eV
h, OL = 0.6736, 0.6847
H0 = 100e3 * h / Mpc * hbar / eV                            # eV
rhoc = 3 * H0**2 * MP**2
rhoL = OL * rhoc
print(f"M_P(red) = {MP:.6e} eV ; H0 = {H0:.5e} eV ; rho_c = {rhoc:.5e} eV^4 = {rhoc/h**2:.5e} h^2")
print(f"rho_Lambda = {rhoL:.5e} eV^4 = ({rhoL**0.25*1e3:.4f} meV)^4 = {rhoL/MP**4:.4e} M_P^4")
for (hh, oo, lab) in ((0.6766, 0.6889, "TT,TE,EE+lowE+lens+BAO"), (0.6727, 0.6834, "TT,TE,EE+lowE")):
    print(f"   alt. Planck column {lab}: rho_L = {oo*3*(100e3*hh/Mpc*hbar/eV)**2*MP**2:.4e}")

print("=" * 70); print("E. Update: NuFIT 6.0 (IC24 w/ SK-atm) + Planck 2018")
d21 = 7.49e-5; d31 = 2.513e-3; d32 = -2.484e-3
res = {}
for order, d3l in (("NO", d31), ("IO", d32)):
    a = onset(rhoL, d21, d3l, order, 4, "AdS"); s = onset(rhoL, d21, d3l, order, 4, "any")
    ms_a = masses(a, d21, d3l, order); ms_0 = masses(0.0, d21, d3l, order)
    res[order] = a
    print(f"Dirac {order}: dS onset {s*1e3:.4f} meV ; AdS onset {a*1e3:.4f} meV ;"
          f" Sigma(m_l=0)={ms_0.sum():.5f} eV ; Sigma(onset)={ms_a.sum():.5f} eV")
    r, Lb, nl = rhomin_AdS(0.0, d21, d3l, order, 2)
    print(f"   Majorana {order} m_l=0: rho_min/rho_obs = {r/rhoL:.3f}  (grid local maxima of B: {nl})")
# 3 sigma corners
print("3-sigma corners (Dm21 in [6.92,8.05]e-5; NO Dm31 in [2.451,2.578]e-3; IO Dm32 in [-2.547,-2.421]e-3)")
for order, rng in (("NO", (2.451e-3, 2.578e-3)), ("IO", (-2.547e-3, -2.421e-3))):
    out = []
    for x21 in (6.92e-5, 8.05e-5):
        for x3 in rng:
            a = onset(rhoL, x21, x3, order, 4, "AdS")
            out.append((a, masses(a, x21, x3, order).sum(), masses(0, x21, x3, order).sum(), x21, x3))
    am = [o[0] for o in out]; sm = [o[1] for o in out]; s0 = [o[2] for o in out]
    print(f"  {order}: AdS onset range {min(am)*1e3:.4f}-{max(am)*1e3:.4f} meV ;"
          f" Sigma_max range {min(sm):.5f}-{max(sm):.5f} ; Sigma_min range {min(s0):.5f}-{max(s0):.5f}")
    for o in out:
        print(f"      Dm21={o[3]:.3e} Dm3l={o[4]:+.3e}: m*={o[0]*1e3:.4f} meV Sigma*={o[1]:.5f}")

print("=" * 70); print("F. rho_min(m1)/rho_obs, Dirac NO (AdS avoided iff rho_L > rho_min)")
for m1 in (0.0, 1e-3, 2e-3, 4e-3, 6e-3, 7.644e-3, 7.7e-3, 1e-2, 2.3e-2, 5e-2):
    r, Lb, nl = rhomin_AdS(m1, d21, d31, "NO", 4)
    rr = max(r, 0.0)
    print(f"m1={m1*1e3:7.3f} meV: rho_min = {r:+.4e} eV^4 = {r/rhoL:+.4e} rho_obs ;"
          f" rho_min^(1/4)/m1 = {(rr**0.25/m1) if m1>0 and rr>0 else float('nan'):.3f}  L*={Lb:.2f}  nloc={nl}")
print("Lambda scaling of AdS onset (Dirac NO): m1*/rho^(1/4)")
for fac in (0.5, 1, 2, 5):
    a = onset(fac * rhoL, d21, d31, "NO", 4, "AdS", hi=0.2)
    print(f"  rho = {fac} rho_obs: m1* = {a*1e3:.4f} meV ; m1*/rho^(1/4) = {a/(fac*rhoL)**0.25:.3f}")

print("=" * 70); print("G. Brute-force V(L) grid cross-check (independent of the B/C maximisation)")
for order, d3l, ml_list in (("NO", d31, (6.40e-3, 6.55e-3, 7.60e-3, 7.70e-3)), ("IO", d32, (2.04e-3, 2.08e-3, 2.52e-3, 2.55e-3))):
    for ml in ml_list:
        ex = V_extrema(rhoL, masses(ml, d21, d3l, order), 4)
        print(f"  Dirac {order} m_l={ml*1e3:.3f} meV: minima (L, F_min) = {[(round(x,3), f'{y:+.3e}') for x,y in ex]}")
for order, d3l in (("NO", d31), ("IO", d32)):
    ex = V_extrema(rhoL, masses(0.0, d21, d3l, order), 2)
    print(f"  Majorana {order} m_l=0: minima = {[(round(x,3), f'{y:+.3e}') for x,y in ex]}")

print("=" * 70); print("H. Monotonicity / uniqueness scans")
for order, d3l in (("NO", d31), ("IO", d32)):
    mls = np.linspace(0, 0.03, 61)
    bmax = [rhomin_AdS(x, d21, d3l, order, 4)[0] for x in mls]
    cmax = [rhomin_any(x, d21, d3l, order, 4)[0] for x in mls]
    nlB = [rhomin_AdS(x, d21, d3l, order, 4)[2] for x in mls[::6]]
    nlC = [rhomin_any(x, d21, d3l, order, 4)[2] for x in mls[::6]]
    print(f"  {order}: max_L B monotone increasing in m_l: {np.all(np.diff(bmax) > 0)} ;"
          f" max_L C monotone: {np.all(np.diff(cmax) > 0)} ; #local maxima of B on grid {nlB} ; of C {nlC}")
    print(f"      C_max >= B_max everywhere: {np.all(np.array(cmax) >= np.array(bmax))}")
