"""
R2-T2: Standard Model on a circle, radion potential, and the non-SUSY AdS conjecture.

Setup (ESTABLISHED, Arkani-Hamed-Dubovsky-Nicolis-Villadoro hep-th/0703067;
Ibanez-Martin-Lozano-Valenzuela (IMV) arXiv:1706.05392 eqs. 3.1-3.12):
  4D -> 3D on a circle of radius R (circumference L = 2 pi R).  In the 3D Einstein
  frame the radion potential is
      V(R) = (2 pi r^3 / R^2) * F(R),   F(R) = rho_Lambda + sum_i (-1)^{s_i} n_i rho(R; m_i)
  (IMV eq. 3.12; r = fixed reference radius, irrelevant for stationary points and signs),
  with the one-loop Casimir energy density per degree of freedom (IMV eq. 3.10)
      rho(R; m) = sum_{k>=1} 2 m^4/(2 pi)^2 * K_2(2 pi R m k)/(2 pi R m k)^2
                -> pi^2/(90 L^4)  as m -> 0.
  Bosons contribute -rho, fermions with periodic boundary conditions +rho.
  Massless bosons: photon (2) + graviton (2) = 4 dof.  Neutrinos: Majorana 2 dof each,
  Dirac 4 dof each.  Everything heavier (electron: m_e R ~ 1e7 at the relevant R) is
  exponentially negligible.

Stationary points: V'(R) = 0  <=>  R F'(R) = 2 F(R).  A local minimum is AdS3 iff F < 0
there, dS3 iff F > 0.

The non-SUSY AdS conjecture (Ooguri-Vafa 1610.01533; CONJECTURE) forbids stable
non-SUSY AdS vacua; applied here (with the stability caveats discussed in IMV sec. 2)
it forbids any AdS3 minimum.

Units: eV throughout.
"""
import numpy as np
import mpmath as mp
from scipy.special import kn

def rho_dof(R, m, kmax=200):
    """Casimir energy density per dof (positive number), circle radius R, mass m [eV, eV^-1]."""
    L = 2 * np.pi * R
    if m == 0.0:
        return np.pi**2 / (90.0 * L**4)
    x = m * L
    k = np.arange(1, kmax + 1)
    # sum_k m^2 K2(k x)/(2 pi^2 k^2 L^2); kn underflows harmlessly for large k x
    terms = kn(2, k * x) / k**2
    return m**2 / (2 * np.pi**2 * L**2) * np.sum(terms)

def F_of_R(R, rho_L, masses, n_per, n_bos=4):
    """F(R) = rho_Lambda - n_bos*rho(0) + n_per * sum_i rho(m_i)."""
    out = rho_L - n_bos * rho_dof(R, 0.0)
    for m in masses:
        out += n_per * rho_dof(R, m)
    return out

def V_of_R(R, *args, **kw):
    return F_of_R(R, *args, **kw) / R**2

def classify(rho_L, masses, n_per, Rgrid=None):
    """Return list of (R, type, F) for local minima of V(R) on a log grid, refined."""
    if Rgrid is None:
        Rgrid = np.logspace(0, 6, 6001)   # 1 .. 1e6 eV^-1  (0.2 um .. 20 cm)
    V = np.array([V_of_R(R, rho_L, masses, n_per) for R in Rgrid])
    mins = []
    for i in range(1, len(Rgrid) - 1):
        if V[i] < V[i - 1] and V[i] < V[i + 1]:
            # refine with golden-section on [R_{i-1}, R_{i+1}]
            a, b = Rgrid[i - 1], Rgrid[i + 1]
            g = (np.sqrt(5) - 1) / 2
            c, d = b - g * (b - a), a + g * (b - a)
            for _ in range(80):
                if V_of_R(c, rho_L, masses, n_per) < V_of_R(d, rho_L, masses, n_per):
                    b = d
                else:
                    a = c
                c, d = b - g * (b - a), a + g * (b - a)
            R0 = 0.5 * (a + b)
            F0 = F_of_R(R0, rho_L, masses, n_per)
            mins.append((R0, 'AdS' if F0 < 0 else 'dS', F0))
    return mins

def spectrum(mlight, ordering, dm21, dm3l):
    """Neutrino masses from the lightest mass; dm3l = Dm^2_31 (NO, >0) or Dm^2_32 (IO, <0)."""
    if ordering == 'NO':
        m1 = mlight
        m2 = np.sqrt(m1**2 + dm21)
        m3 = np.sqrt(m1**2 + dm3l)
    else:
        m3 = mlight
        m2 = np.sqrt(m3**2 - dm3l)          # dm3l = m3^2 - m2^2 < 0
        m1 = np.sqrt(m2**2 - dm21)
    return [m1, m2, m3]

def status(mlight, ordering, dm21, dm3l, rho_L, nature):
    n_per = 4 if nature == 'Dirac' else 2
    ms = classify(rho_L, spectrum(mlight, ordering, dm21, dm3l), n_per)
    types = [t for _, t, _ in ms]
    if 'AdS' in types:
        return 'AdS', ms
    if 'dS' in types:
        return 'dS', ms
    return 'none', ms

def threshold(ordering, dm21, dm3l, rho_L, nature, kind, lo=0.0, hi=0.05, tol=1e-6):
    """Smallest lightest-mass at which a minimum of type `kind` (AdS or any) appears (bisection,
    assuming monotonic onset as the mass grows, which we check on a scan)."""
    def has(m):
        s, _ = status(m, ordering, dm21, dm3l, rho_L, nature)
        return (s == 'AdS') if kind == 'AdS' else (s in ('AdS', 'dS'))
    if has(lo):
        return 0.0
    if not has(hi):
        return None
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        if has(mid):
            hi = mid
        else:
            lo = mid
        if hi - lo < tol:
            break
    return 0.5 * (lo + hi)

def sum_masses(mlight, ordering, dm21, dm3l):
    return sum(spectrum(mlight, ordering, dm21, dm3l))

if __name__ == '__main__':
    out = []
    def log(s=''):
        print(s); out.append(s)

    log("=== check: Casimir density per dof, massless limit and IMV eq. 3.11 small-mR expansion ===")
    R = 50.0
    L = 2 * np.pi * R
    log(f"rho(R,0)*L^4 = {rho_dof(R,0.0)*L**4:.10f}   pi^2/90 = {np.pi**2/90:.10f}")
    for mR in (0.01, 0.05, 0.1):
        m = mR / R
        exact = rho_dof(R, m)
        series = np.pi**2/(90*L**4) - np.pi**2/(6*L**4)*mR**2 + np.pi**2/(48*L**4)*mR**4
        log(f"mR={mR}: exact {exact*L**4:.10f}  IMV(3.11) series {series*L**4:.10f}  rel.diff {abs(exact-series)/exact:.2e}")
    # massless bosons + Lambda: maximum at R_max = (1/(120 pi^2 Lambda))^(1/4) (IMV eq. 3.6)
    rhoL_imv = 2.6e-11  # eV^4  (IMV: 2.6e-47 GeV^4)
    Rg = np.logspace(1, 13, 20001)
    Vb = np.array([(rhoL_imv - 4*rho_dof(r, 0.0))/r**2 for r in Rg])
    log(f"bosons+Lambda only: numerical R_max = {Rg[np.argmax(Vb)]:.4e} eV^-1 ; IMV eq. 3.6 (1/(120 pi^2 Lambda))^(1/4) = {(1/(120*np.pi**2*rhoL_imv))**0.25:.4e} eV^-1 (IMV eq. 3.7 quotes 7.55e10 GeV^-1 = 7.55e1 eV^-1)")

    log("\n=== (A) reproduce IMV 3D results: Lambda = 2.6e-11 eV^4, PDG-2016-type splittings (assumed: dm21 = 7.53e-5, |dm32| = 2.44e-3 NO / 2.51e-3 IO) ===")
    dm21_p = 7.53e-5
    dm31_NO_p = 2.44e-3 + dm21_p
    dm32_IO_p = -2.51e-3
    for nature in ('Majorana', 'Dirac'):
        for ordering, dm3l in (('NO', dm31_NO_p), ('IO', dm32_IO_p)):
            s0, ms0 = status(0.0, ordering, dm21_p, dm3l, rhoL_imv, nature)
            tA = threshold(ordering, dm21_p, dm3l, rhoL_imv, nature, 'AdS')
            tV = threshold(ordering, dm21_p, dm3l, rhoL_imv, nature, 'any')
            log(f"  {nature:8s} {ordering}: m_light=0 -> {s0:4s};  any-minimum onset m* = {tV if tV is None else f'{1e3*tV:.3f} meV'};  AdS onset m* = {tA if tA is None else f'{1e3*tA:.3f} meV'}")
    log("  IMV (sec. 3.4, Table 1) report: Majorana -> AdS always; Dirac NO: none < 6.7, dS 6.7-7.73, AdS > 7.73 meV; Dirac IO: none < 2.1, dS 2.1-2.56, AdS > 2.56 meV")
    log("  ADNV (hep-th/0703067) report: Dirac NO: none < 7.1, dS 7.1-8.3, AdS > 8.3 meV; IO: none < 2.5, dS 2.5-3.1, AdS > 3.1 meV")

    with open(__file__.replace('.py', '_partA_output.txt'), 'w') as f:
        f.write("\n".join(out) + "\n")
