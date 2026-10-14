"""
pwcore.py -- exact Gegenbauer partial-wave machinery for T1.

Conventions (OUR ATTEMPT, standard):
  A residue polynomial R(t) at a pole s = s_n, external masses m (all equal), is expanded as
      R(t(x)) = sum_j c_j C_j^{nu}(x),   nu = (D-3)/2,   x = cos(theta),
      t = kappa (x - 1),  kappa = (s_n - 4 m^2)/2.
  Unitarity (factorization on positive-norm states) requires
      rho_j := c_j / kappa^j >= 0      (all j)
  [for s_n > 4m^2 this is c_j >= 0; for a bound state below threshold, kappa<0, the
   physical coupling-squared is c_j/kappa^j: (|p||p'|)^j C_j(cos) is the polynomial structure].
  Writing R(t) = sum_k rt_k y^k with y = t + kappa (= kappa x), and
      x^k = k!/2^k sum_m (nu+k-2m)/(m! (nu)_{k-m+1}) C^{nu}_{k-2m}(x)   (DLMF 18.18.? standard),
  one finds rho_j = (nu+j)/(nu)_{j+1} * sum_m rt_{j+2m} kappa^{2m} (j+2m)!/(2^{j+2m} m! (nu+j+1)_m).
  The prefactor is >0 for nu>0, so   sign(rho_j) = sign(P_j(nu)) with the POLYNOMIAL
      P_j(nu) = sum_{m=0}^{M} rt_{j+2m} kappa^{2m} (j+2m)!/(2^{j+2m} m!) prod_{i=m}^{M-1}(nu+j+1+i),
  M = floor((deg R - j)/2).  Everything is exact (python-flint fmpq / fmpq_poly);
  real roots are isolated rigorously with arb (fmpz_poly.complex_roots).
"""
from math import factorial
from flint import fmpq, fmpq_poly, fmpz_poly, arb

T = fmpq_poly([0, 1])  # the variable


def poch_poly(c, n):
    """prod_{k=0}^{n-1} (t + c + k) as fmpq_poly in t (c rational)."""
    p = fmpq_poly([1])
    c = fmpq(c)
    for k in range(n):
        p = p * fmpq_poly([c + k, 1])
    return p


def shift_poly(R, kappa):
    """Return R(y - kappa) as polynomial in y."""
    coeffs = R.coeffs()
    out = fmpq_poly([0])
    lin = fmpq_poly([-fmpq(kappa), 1])
    for c in reversed(coeffs):
        out = out * lin + fmpq_poly([c])
    return out


def sign_polys(R, kappa):
    """dict j -> P_j(nu) (fmpq_poly); sign(P_j(nu)) = sign(rho_j) for nu > 0."""
    kappa = fmpq(kappa)
    Ry = shift_poly(R, kappa)
    deg = Ry.degree()
    rt = [Ry[k] for k in range(deg + 1)]
    k2 = kappa * kappa
    P = {}
    for j in range(deg + 1):
        M = (deg - j) // 2
        Q = fmpq_poly([1])
        acc = fmpq_poly([0])
        for m in range(M, -1, -1):
            k = j + 2 * m
            if rt[k] != 0:
                coef = rt[k] * k2 ** m * fmpq(factorial(k), 2 ** k * factorial(m))
                acc = acc + Q * coef
            Q = Q * fmpq_poly([j + m, 1])  # (nu + j + m) = factor i=m-1 for the next m
        P[j] = acc
    return P


def rho_values(R, kappa, D):
    """Exact rho_j (including positive prefactor, WITHOUT the kappa^j) at rational D.
    c_j = kappa^j * rho_j are the actual Gegenbauer coefficients of R(t(x))."""
    nu = (fmpq(D) - 3) / 2
    kappa = fmpq(kappa)
    Ry = shift_poly(R, kappa)
    deg = Ry.degree()
    rt = [Ry[k] for k in range(deg + 1)]
    out = {}
    for j in range(deg + 1):
        s = fmpq(0)
        poch = fmpq(1)  # (nu+j+1)_m
        for m in range(0, (deg - j) // 2 + 1):
            k = j + 2 * m
            s += rt[k] * kappa ** (2 * m) * fmpq(factorial(k), 2 ** k * factorial(m)) / poch
            poch *= (nu + j + 1 + m)
        pre = (nu + j)
        pj = fmpq(1)
        for i in range(j + 1):
            pj *= (nu + i)
        out[j] = s * pre / pj
    return out


def gegenbauer_coeffs(R, kappa, D):
    """Actual coefficients c_j of R(t(x)) in C_j^{(D-3)/2}(x)."""
    rv = rho_values(R, kappa, D)
    kappa = fmpq(kappa)
    return {j: v * kappa ** j for j, v in rv.items()}


def first_neg_root(P):
    """For polynomial P(nu): returns (status, root) where
       status='neg0' if P<0 just above nu=0 (fails for all D>3, up to roots at 0),
       status='ok'   if P>=0 for all nu>0 (root=None),
       status='root' if P>0 just above 0 and first sign change to negative at root (arb)."""
    if P.is_zero():
        return 'ok', None
    # remove factors nu^k (roots at nu=0 do not matter for nu>0)
    coeffs = P.coeffs()
    i0 = 0
    while coeffs[i0] == 0:
        i0 += 1
    Pr = fmpq_poly(coeffs[i0:])
    # sign at nu=0+ is sign of Pr(0)
    s0 = 1 if Pr[0] > 0 else -1
    roots = []
    if Pr.degree() >= 1:
        for r, mult in Pr.numer().complex_roots():
            if r.imag.is_zero() and mult % 2 == 1:
                re = r.real
                if re > 0:
                    roots.append(re)
                elif not (re < 0):
                    raise ValueError("root ball straddles 0")
    roots.sort(key=lambda z: float(z.mid()))
    if s0 < 0:
        return 'neg0', None
    if not roots:
        return 'ok', None
    return 'root', roots[0]


def level_crit(R, kappa):
    """Critical nu for one level: positivity of all rho_j holds for 0 < nu <= nu_c.
    Returns (nu_c as arb or None=infinity or 0 if fails near nu=0, j_responsible)."""
    P = sign_polys(R, kappa)
    best, jbest = None, None
    for j, Pj in P.items():
        st, r = first_neg_root(Pj)
        if st == 'neg0':
            return arb(0), j
        if st == 'root':
            if best is None or float(r.mid()) < float(best.mid()):
                best, jbest = r, j
    return best, jbest


def level_ok_at(R, kappa, D, strict=False):
    """Exact check: all rho_j >= 0 at rational D. Returns (ok, list of (j, rho_j) negative, zeros)."""
    nu = (fmpq(D) - 3) / 2
    P = sign_polys(R, kappa)
    neg, zer = [], []
    for j, Pj in P.items():
        v = Pj(nu)
        if v < 0:
            neg.append(j)
        elif v == 0 and not Pj.is_zero():
            zer.append(j)
    return (len(neg) == 0), neg, zer


def D_of_nu(nu):
    return 3 + 2 * nu
