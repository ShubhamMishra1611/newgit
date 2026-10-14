"""Claim 6/7 (closed): exact per-level Dcrit for
  VS stripped        : residue [(t+1)_{n-1}]^2, massless, kappa=n/2
  four-dilaton (2102.08951 sec 3.6, ar5iv eq.(58)): residue P*[(t+1)_{n-1}]^2, P=(st)^2+(tu)^2+(su)^2 at s=n
  closed bosonic tachyon (KLT: open^2): residue [(t+2)_n]^2, m^2=-1, kappa=(n+3)/2
  S-families (VS-ratio, al_s+al_t+al_u=S): residue [(al_t+1)_{n-1-S}]^2, m^2=(S-3a)/4
  CHR closed lambda family (2408.03362 eq.(13)), massless, mu(n)=(n+lam-1)/lam
Usage: python3 rt1_closed.py NMAX
"""
import sys; sys.path.insert(0, '.')
from rt1core import *
from fractions import Fraction as Fr
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 40

def times_x(c, nu):
    new = [Q(0)] * (len(c) + 1)
    for j, cj in enumerate(c):
        if cj == 0: continue
        d = 2 * (j + nu)
        new[j + 1] += cj * (j + 1) / d
        if j >= 1:
            new[j - 1] += cj * (j + 2 * nu - 1) / d
    return new

def multiply_poly(c, pcoeffs, nu):
    """multiply Gegenbauer series c by polynomial sum_k p_k x^k."""
    out = [Q(0)] * (len(c) + len(pcoeffs) - 1)
    cur = c
    for k, pk in enumerate(pcoeffs):
        if k > 0:
            cur = times_x(cur, nu)
        if pk != 0:
            for j, v in enumerate(cur):
                out[j] += pk * v
    return out

def dilaton_ok(n, D):
    D = q(D); nu = (D - 3) / 2
    kap = Q(n, 2)
    facs = [(Q(i) - kap, kap) for i in range(1, n)]
    facs = facs + facs
    c = gegen_from_linear_factors(facs, nu, zero=Q(0))
    # P(x) with s=n, t=kap(x-1), u=-n-t
    import sympy as sp
    x = sp.symbols('x'); t = sp.Rational(n, 2) * (x - 1); s = n; u = -n - t
    P = sp.Poly(sp.expand((s*t)**2 + (t*u)**2 + (s*u)**2), x)
    pc = [Q(int(sp.Rational(P.coeff_monomial(x**k)).p), int(sp.Rational(P.coeff_monomial(x**k)).q)) for k in range(P.degree() + 1)]
    c = multiply_poly(c, pc, nu)
    r = [cj / kap ** j for j, cj in enumerate(c)]
    bad = [j for j, v in enumerate(r) if v < 0]
    return len(bad) == 0, bad

def bisect(okfun, lo=Fr(3) + Fr(1, 10**9), hi=Fr(400), tol=Fr(1, 10**10)):
    ok, bad = okfun(hi)
    if ok: return None, []
    ok, bad = okfun(lo)
    if not ok: return Fr(3), bad
    while hi - lo > tol:
        mid = Fr(round(((lo + hi) / 2) * 2**40), 2**40)
        ok, bad = okfun(mid)
        if ok: lo = mid
        else: hi = mid
    return lo, okfun(hi)[1]

def scan(name, mk, nmax):
    res = []
    for n in range(1, nmax + 1):
        f, k = mk(n)
        d, bad = bisect(lambda D: ok_factors(f, k, D))
        if d is not None:
            res.append((d, n, bad))
    res.sort()
    print(f"{name}: Dcrit = {float(res[0][0]):.10f} set by level {res[0][1]} spins {res[0][2]}; next: "
          f"{[(r[1], round(float(r[0]), 4)) for r in res[1:4]]}; Dcrit(n={nmax}) = "
          f"{[round(float(r[0]),4) for r in res if r[1]==nmax]}", flush=True)
    return res

if __name__ == '__main__':
    scan('VS stripped (S=0,a=0)', lambda n: closed_facs(n, 0, 0, 0), NMAX)
    scan('closed bosonic tachyon (S=-1,a=1,m2=-1)', lambda n: closed_facs(n, 1, -1, -1), NMAX)
    # dilaton
    res = []
    for n in range(1, min(NMAX, 40) + 1):
        d, bad = bisect(lambda D: dilaton_ok(n, D))
        if d is not None: res.append((d, n, bad))
    res.sort()
    print(f"four-dilaton: Dcrit = {float(res[0][0]):.12f} level {res[0][1]} spins {res[0][2]}; next {[(r[1], round(float(r[0]),4)) for r in res[1:4]]}", flush=True)
    for S in (-1, 0):
        for a in (Fr(-1, 2), Fr(0), Fr(1, 2), Fr(1), Fr(3, 2)):
            m2 = (S - 3 * a) / 4
            scan(f'S={S} a={a} (m2={m2})', lambda n, a=a, m2=m2, S=S: closed_facs(n, a, m2, S), min(NMAX, 40))
    # CHR closed lambda family, via own residue eq.(13): factors (lam*t + (1+lam)/2 + k - 1), massless, s_n=(n+lam-1)/lam
    for lam in (Fr(1, 16), Fr(1, 2), Fr(1), Fr(3, 2), Fr(2), Fr(3)):
        def mk(n, lam=lam):
            sn = Fr(n) + lam - 1
            sn = sn / lam
            kap = q(sn / 2)
            L = q(lam)
            f = [(L * (-kap) + (1 + L) / 2 + kk - 1, L * kap) for kk in range(1, n)]
            return f + f, kap
        scan(f'CHR closed lambda={lam}', mk, min(NMAX, 40))
        def mkp(n, lam=lam):
            sn = (Fr(n) + lam - 1) / lam
            kap = q(sn / 2); L = q(lam)
            return [(L * (-kap) + (1 + L) / 2 + kk - 1, L * kap) for kk in range(1, n)], kap
        scan(f'CHR planar lambda={lam}', mkp, min(NMAX, 40))
