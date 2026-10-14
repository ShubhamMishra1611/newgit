"""
R2-T2 part B: update the 3D bounds with NuFIT 6.0 (arXiv:2410.05380, Table 1, 'IC24 with SK
atmospheric data' variant) and Planck 2018 (arXiv:1807.06209, Table 2, TT,TE,EE+lowE+lensing:
H0 = 67.36, Omega_Lambda = 0.6847), and translate into Sum m_nu for comparison with DESI DR2
(arXiv:2503.14744 eq. 15: Sum m_nu < 0.0642 eV (95%, LCDM); w0wa: < 0.163 eV; lightest mass
m_l < 0.023 eV).
"""
import numpy as np
from r2t2_casimir import status, threshold, sum_masses, spectrum, classify

out = []
def log(s=''):
    print(s); out.append(s)

# Planck 2018 vacuum energy density in eV^4
Mp = 2.435e27                    # reduced Planck mass [eV]
H100 = 2.1332e-33                # 100 km/s/Mpc in eV
h, OL = 0.6736, 0.6847
rho_c = 3 * (h * H100)**2 * Mp**2
rhoL = OL * rho_c
log(f"Planck 2018: rho_crit = {rho_c:.4e} eV^4, rho_Lambda = {rhoL:.4e} eV^4 = ({rhoL**0.25*1e3:.4f} meV)^4 ; in reduced Planck units rho_L/Mp^4 = {rhoL/Mp**4:.3e}")

# NuFIT 6.0 (IC24 + SK-atm): best fit and 3 sigma ranges
dm21 = (7.49e-5, 6.92e-5, 8.05e-5)
dm31_NO = (2.513e-3, 2.451e-3, 2.578e-3)
dm32_IO = (-2.484e-3, -2.547e-3, -2.421e-3)

log("\n=== 3D thresholds, NuFIT 6.0 best fit, Planck 2018 rho_Lambda ===")
res = {}
for nature in ('Majorana', 'Dirac'):
    for ordering, dm3l in (('NO', dm31_NO[0]), ('IO', dm32_IO[0])):
        s0, _ = status(0.0, ordering, dm21[0], dm3l, rhoL, nature)
        tA = threshold(ordering, dm21[0], dm3l, rhoL, nature, 'AdS')
        tV = threshold(ordering, dm21[0], dm3l, rhoL, nature, 'any')
        res[(nature, ordering)] = (tV, tA)
        line = f"  {nature:8s} {ordering}: m_light=0 -> {s0:4s}"
        if tA is not None and tA > 0:
            line += f"; dS-min onset {1e3*tV:.3f} meV; AdS onset {1e3*tA:.3f} meV"
            line += f"; => conjecture requires m_light < {1e3*tA:.3f} meV, i.e. Sum m_nu < {sum_masses(tA, ordering, dm21[0], dm3l):.5f} eV"
            line += f" (Sum at m_light=0: {sum_masses(0.0, ordering, dm21[0], dm3l):.5f} eV)"
        log(line)

log("\n=== robustness: AdS onset over NuFIT 6.0 3-sigma corners of (dm21, |dm3l|) ===")
for ordering, dmset in (('NO', dm31_NO), ('IO', dm32_IO)):
    vals = []
    for a in dm21:
        for b in dmset:
            t = threshold(ordering, a, b, rhoL, 'Dirac', 'AdS')
            vals.append((t, sum_masses(t, ordering, a, b)))
    ts = [1e3*v[0] for v in vals]; ss = [v[1] for v in vals]
    log(f"  Dirac {ordering}: AdS onset m_light in [{min(ts):.3f}, {max(ts):.3f}] meV ; Sum m_nu bound in [{min(ss):.5f}, {max(ss):.5f}] eV")

log("\n=== robustness: dependence on rho_Lambda (Dirac NO, best-fit splittings) ===")
for f in (0.5, 0.9, 1.0, 1.1, 2.0, 5.0):
    t = threshold('NO', dm21[0], dm31_NO[0], f * rhoL, 'Dirac', 'AdS')
    log(f"  rho_L x {f:4.1f}: AdS onset m1 = {1e3*t:.3f} meV ; (m1/rho_L^1/4 = {t/(f*rhoL)**0.25:.4f})")

log("\n=== Majorana: how large would rho_Lambda have to be to avoid AdS3 (lightest mass = 0)? ===")
for ordering, dm3l in (('NO', dm31_NO[0]), ('IO', dm32_IO[0])):
    lo, hi = rhoL, rhoL
    while status(0.0, ordering, dm21[0], dm3l, hi, 'Majorana')[0] == 'AdS':
        hi *= 2
    lo = hi / 2
    for _ in range(40):
        mid = np.sqrt(lo * hi)
        if status(0.0, ordering, dm21[0], dm3l, mid, 'Majorana')[0] == 'AdS':
            lo = mid
        else:
            hi = mid
    log(f"  Majorana {ordering}, m_light = 0: AdS3 disappears for rho_Lambda > {hi:.4e} eV^4 = {hi/rhoL:.2f} x rho_obs")
log("  (IMV fig. 5 statement: 'around 7 times higher than the experimental Lambda_4' in NO, 'several orders of magnitude' in IO)")

log("\n=== Lower bound on rho_Lambda as a function of the lightest Dirac mass (NO): rho_min(m1) ===")
for m1 in (1e-3, 2e-3, 4e-3, 7.7e-3, 1e-2, 2.3e-2):
    lo, hi = 1e-20, 1e-6
    for _ in range(60):
        mid = np.sqrt(lo * hi)
        if status(m1, 'NO', dm21[0], dm31_NO[0], mid, 'Dirac')[0] == 'AdS':
            lo = mid
        else:
            hi = mid
    log(f"  m1 = {1e3*m1:6.2f} meV: rho_min = {hi:.4e} eV^4 = {hi/rhoL:.3e} rho_obs ; rho_min^(1/4)/m1 = {hi**0.25/m1:.4f}")

log("\n=== Comparison with DESI DR2 (2503.14744) ===")
tNO = res[('Dirac', 'NO')][1]; tIO = res[('Dirac', 'IO')][1]
log(f"  conjecture (3D, SM + Dirac): NO: m1 < {1e3*tNO:.2f} meV, Sum in [{sum_masses(0,'NO',dm21[0],dm31_NO[0]):.4f}, {sum_masses(tNO,'NO',dm21[0],dm31_NO[0]):.4f}] eV")
log(f"                               IO: m3 < {1e3*tIO:.2f} meV, Sum in [{sum_masses(0,'IO',dm21[0],dm32_IO[0]):.4f}, {sum_masses(tIO,'IO',dm21[0],dm32_IO[0]):.4f}] eV")
log("  DESI DR2 + CMB, LCDM: Sum < 0.0642 eV (95%); m_l < 0.023 eV;   w0waCDM: Sum < 0.163 eV (95%)")

with open('r2t2_update_output.txt', 'w') as f:
    f.write("\n".join(out) + "\n")
