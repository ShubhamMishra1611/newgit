"""RT1 independent partial-wave core (written from scratch, not derived from pwcore.py).

Kinematics (equal-mass scalars m^2, mostly-minus, s+t+u=4m^2):
  CM: s = 4(p^2+m^2), t = -2p^2(1-x), x = cos(theta)  =>  t = kappa (x-1), kappa = 2p^2 = (s-4m^2)/2.
  Below threshold (s<4m^2) kappa<0; x := 1 + t/kappa is the analytic continuation of cos(theta).

Unitarity at a pole s = M^2 (factorisation, real cubic couplings g_j):
  Res_{s=M^2} = sum_j g_j^2 N_j (s-4m^2)^j C_j^{nu}(x),  N_j>0,  nu=(D-3)/2
  so with R(x) = sum_j c_j C_j^nu(x) the condition is rho_j := c_j / kappa^j >= 0.
  (derivation + explicit Feynman-rule check in rt1_subthreshold.py)

Gegenbauer coefficients are obtained by multiplying linear factors (b + k x) one at a time
using the three-term recurrence  x C_j = [(j+1) C_{j+1} + (j+2nu-1) C_{j-1}] / (2(j+nu)).
Works with Fraction (exact) or mpmath mpf (high precision).
"""
from fractions import Fraction as Fr
import mpmath as mp


def gegen_from_linear_factors(factors, nu, zero=None):
    """factors: list of (b,k) meaning polynomial prod (b + k x). Returns list c_j (j=0..deg)."""
    if zero is None:
        zero = Fr(0) if isinstance(nu, Fr) else mp.mpf(0)
    c = [zero + 1]
    for (b, k) in factors:
        new = [zero] * (len(c) + 1)
        for j, cj in enumerate(c):
            if cj == 0:
                continue
            new[j] += b * cj
            if k != 0:
                d = 2 * (j + nu)
                new[j + 1] += k * cj * (j + 1) / d
                if j >= 1:
                    new[j - 1] += k * cj * (j + 2 * nu - 1) / d
        c = new
    return c


def open_factors(n, a, m2, c=0, num=Fr):
    """Residue of Gamma(-al_s)Gamma(-al_t)/Gamma(c-al_s-al_t) at al_s=n, al(x)=a+x (alpha'=1):
       R_n = (al_t+1)_{n-c}/n!  (overall (-1)^c dropped).  Returns (factors, kappa)."""
    a = num(a); m2 = num(m2)
    sn = n - a
    kap = (sn - 4 * m2) / 2
    facs = []
    for i in range(1, n - c + 1):
        # al_t + i = a + i + kap(x-1)
        facs.append((a + i - kap, kap))
    return facs, kap


def rho_signs(c, kap):
    """signs of rho_j = c_j/kap^j (kap != 0)."""
    out = []
    for j, cj in enumerate(c):
        s = (cj > 0) - (cj < 0)
        if kap < 0 and j % 2 == 1:
            s = -s
        out.append(s)
    return out


def rho_vals(c, kap):
    return [cj / kap ** j for j, cj in enumerate(c)]


def level_ok(factors, kap, D, num=Fr, tol=0):
    """True if all rho_j >= -tol at dimension D."""
    nu = (num(D) - 3) / 2
    if kap == 0:
        # residue independent of x: only j=0
        prod = num(1)
        for b, k in factors:
            prod *= b
        return prod >= -tol, [prod]
    c = gegen_from_linear_factors(factors, nu)
    r = rho_vals(c, kap)
    return all(x >= -tol for x in r), r


def mp_level_min(factors, kap, D):
    """min_j rho_j normalised by max|rho| (mpmath) -> sign indicator, and argmin j."""
    nu = (mp.mpf(D) - 3) / 2
    fac = [(mp.mpf(b.numerator) / b.denominator if isinstance(b, Fr) else mp.mpf(b),
            mp.mpf(k.numerator) / k.denominator if isinstance(k, Fr) else mp.mpf(k)) for b, k in factors]
    kp = mp.mpf(kap.numerator) / kap.denominator if isinstance(kap, Fr) else mp.mpf(kap)
    c = gegen_from_linear_factors(fac, nu, zero=mp.mpf(0))
    r = [cj / kp ** j for j, cj in enumerate(c)]
    scale = max(abs(x) for x in r)
    rn = [x / scale for x in r]
    jmin = min(range(len(rn)), key=lambda j: rn[j])
    return rn[jmin], jmin, rn


def EPS():
    # relative threshold below which a normalised rho is treated as an (identically) zero coefficient
    return mp.mpf(10) ** (-(mp.mp.dps * 2) // 3)


def Dcrit_level(factors, kap, Dlo=3, Dhi=200, iters=70, dps=None):
    """Largest D such that all rho_j >= 0 (uses monotonicity in D). mpmath bisection.
       Returns (Dcrit, spin that fails just above)."""
    if kap == 0:
        return mp.inf, None
    lo = mp.mpf(Dlo) + mp.mpf(10) ** -12
    hi = mp.mpf(Dhi)
    v, j, _ = mp_level_min(factors, kap, hi)
    if v >= -EPS():
        return mp.inf, None
    v, j, _ = mp_level_min(factors, kap, lo)
    if v < -EPS():
        return mp.mpf(3), j
    for _ in range(iters):
        mid = (lo + hi) / 2
        v, j, _ = mp_level_min(factors, kap, mid)
        if v >= -EPS():
            lo = mid
        else:
            hi = mid
    _, j, rn = mp_level_min(factors, kap, hi)
    # report the most negative spin just above Dcrit
    return lo, j
