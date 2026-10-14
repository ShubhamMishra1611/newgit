"""
RR2 independent re-derivation of the SM-on-a-circle radion potential bounds.
Written from first principles BEFORE reading r2t2_casimir.py / r2t2_update.py.

Units: eV.  L = 2 pi R = circumference.

Physics (derived here, see RR2_verdicts.md section A):
 * 4D -> 3D Einstein frame: ds4^2 = (r/R)^2 g3 + R^2 dphi^2  =>  V3(R) = 2 pi r^3 rho4(R) / R^2.
 * One-loop Casimir energy density per dof of a field periodic on the circle =
   minus the thermal free-energy density at T = 1/L (Euclidean rotation; T^t_t of the
   thermal state with the compact direction relabelled):
      rho_dof(L,m) = -f_B(T=1/L, m) = -(T/2pi^2) Int p^2 ln(1 - exp(-E/T)) dp  (>0)
                   = (m^4/2pi^2) sum_k K2(kmL)/(kmL)^2                          (Bessel form)
   bosons enter with -rho, periodic fermions with +rho, antiperiodic fermions with -(7/8-type) rho.
 * F(L) = rho_L - B(L),  B(L) = 4 rho(L,0) - sum_i n_i rho(L,m_i)   [n_i = 4 Dirac, 2 Majorana]
   V ~ F/L^2.  V'(L) = 2 (C(L) - rho_L)/L^3 with C = B - L B'/2.
 * AdS min exists  <=>  min_L F < 0  <=>  rho_L < max_L B(L)        (F>0 at both ends)
 * any min exists  <=>  rho_L < max_L C(L)    (C -> -inf at small L, C -> 0+ at large L,
                                              so C crosses rho_L upward (min) then downward (max))
   -> both onsets reduce to 1D MAXIMISATIONS; no grid search for minima of V needed.
"""
import numpy as np
from scipy.special import kn, kv
from scipy.integrate import quad
from scipy.optimize import brentq, minimize_scalar

PI = np.pi

# ---------------------------------------------------------------- Casimir pieces
def rho0(L):
    return PI**2 / (90.0 * L**4)

def c0(L):
    return 3.0 * rho0(L)

def rho_bessel(L, m, kmax=2000):
    if m == 0.0:
        return rho0(L)
    k = np.arange(1, kmax + 1)
    y = k * m * L
    return m**4 / (2 * PI**2) * np.sum(kv(2, y) / y**2)

def c_bessel(L, m, kmax=2000):
    """c = rho - L rho'/2 per dof (uses d/dy[K2(y)/y^2] = -K3(y)/y^2)."""
    if m == 0.0:
        return c0(L)
    k = np.arange(1, kmax + 1)
    y = k * m * L
    return m**4 / (2 * PI**2) * np.sum(kv(2, y) / y**2 + kv(3, y) / (2 * y))

def rho_thermal(L, m):
    """Independent route: -f_B(T=1/L) by direct momentum integral."""
    T = 1.0 / L
    f = lambda p: p**2 * np.log1p(-np.exp(-np.sqrt(p * p + m * m) / T))
    val, err = quad(f, 0, 60 * T + 0 * m, limit=400, epsabs=0, epsrel=1e-12)
    return -T / (2 * PI**2) * val

# ---------------------------------------------------------------- spectrum
def masses(ml, d21, d3l, order):
    """NuFIT convention: d3l = Dm^2_31 (NO) or Dm^2_32 (IO, negative)."""
    if order == "NO":
        m1 = ml; m2 = np.sqrt(ml**2 + d21); m3 = np.sqrt(ml**2 + d3l)
    else:
        m3 = ml; m2 = np.sqrt(ml**2 - d3l); m1 = np.sqrt(m2**2 - d21)
    return np.array([m1, m2, m3])

def B(L, ms, nf):
    return 4 * rho0(L) - nf * sum(rho_bessel(L, m) for m in ms)

def C(L, ms, nf):
    return 4 * c0(L) - nf * sum(c_bessel(L, m) for m in ms)

def maxfun(fun, ms, nf, Llo=1e-1, Lhi=1e5, N=1500):
    """Global max of fun(L) on a log grid, then golden refinement around the best cell.
    Returns (max value, argmax L, number of local maxima found on the grid)."""
    Ls = np.geomspace(Llo, Lhi, N)
    vals = np.array([fun(L, ms, nf) for L in Ls])
    i = int(np.argmax(vals))
    nloc = int(np.sum((vals[1:-1] > vals[:-2]) & (vals[1:-1] > vals[2:])))
    a, b = Ls[max(i - 1, 0)], Ls[min(i + 1, N - 1)]
    res = minimize_scalar(lambda lx: -fun(np.exp(lx), ms, nf), bounds=(np.log(a), np.log(b)),
                          method="bounded", options={"xatol": 1e-12})
    return -res.fun, np.exp(res.x), nloc

def rhomin_AdS(ml, d21, d3l, order, nf):
    return maxfun(B, masses(ml, d21, d3l, order), nf)

def rhomin_any(ml, d21, d3l, order, nf):
    return maxfun(C, masses(ml, d21, d3l, order), nf)

def onset(rhoL, d21, d3l, order, nf, which="AdS", lo=1e-5, hi=0.05):
    f = (lambda ml: rhomin_AdS(ml, d21, d3l, order, nf)[0] - rhoL) if which == "AdS" else \
        (lambda ml: rhomin_any(ml, d21, d3l, order, nf)[0] - rhoL)
    return brentq(f, lo, hi, xtol=1e-9)

# ---------------------------------------------------------------- brute-force cross-check of V
def V_extrema(rhoL, ms, nf, Llo=1.0, Lhi=1e5, N=200000):
    """Brute force: dense grid of V(L) = F/L^2 (vectorised Bessel sum), list local minima."""
    Ls = np.geomspace(Llo, Lhi, N)
    k = np.arange(1, 401)[:, None]
    Bv = 4 * rho0(Ls)
    for m in ms:
        y = k * m * Ls[None, :]
        Bv = Bv - nf * m**4 / (2 * PI**2) * np.sum(kv(2, y) / y**2, axis=0)
    V = (rhoL - Bv) / Ls**2
    idx = np.where((V[1:-1] < V[:-2]) & (V[1:-1] < V[2:]))[0] + 1
    return [(Ls[i], V[i] * Ls[i]**2) for i in idx]   # (L_min, F at min: sign = sign of V)
