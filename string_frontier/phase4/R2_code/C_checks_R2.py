"""
R2 independent checks of C1 (Bousso-Polchinski), C2 (DESI exponential-quintessence tension), C3 (dark dimension).
"""
import numpy as np
import mpmath as mp
from math import gamma, pi, log10, sqrt

# ---------------- C1: BP hep-th/0004134 eqs. (2.24)-(2.26) ----------------
def omega(Jm1):          # area of unit S^{J-1}:  2 pi^{J/2}/Gamma(J/2)
    J = Jm1 + 1
    return 2*mp.pi**(mp.mpf(J)/2)/mp.gamma(mp.mpf(J)/2)

def log10_q(J, lam_bare, dlam, D):
    # equal charges: q^J = omega_{J-1}/D |2 lam|^{J/2-1} dlam      (2.24)
    v = mp.log10(omega(J-1)) - mp.log10(D) + (mp.mpf(J)/2-1)*mp.log10(abs(2*lam_bare)) + mp.log10(dlam)
    return v/J

print("C1: BP equal-charge estimate at J=100, Delta lambda = 1e-120")
for lab, lb in [("lambda_bare=1 (phase-1 choice)", 1), ("|2 lambda_bare|=1", mp.mpf(1)/2),
                ("(sqrt2 kappa)^4 lambda_bare=1 with kappa=1 -> lambda_bare=1/4", mp.mpf(1)/4)]:
    print(f"   {lab:60s}: log10 q = {mp.nstr(log10_q(100, lb, mp.mpf(10)**-120, 1),5)} (D=1), "
          f"{mp.nstr(log10_q(100, lb, mp.mpf(10)**-120, mp.mpf(2)**100),5)} (D=2^J)")

# brute-force shell counts vs (2.24) expectation  N_exp = omega_{J-1} r^{J-1} dr / prod q
def shell_count(q, r2, dr2):
    J = len(q)
    nmax = [int(np.floor(np.sqrt(r2+dr2)/qi)) for qi in q]
    # meet-in-the-middle over first J-1 coords, last coordinate solved analytically
    grids = np.meshgrid(*[np.arange(-n, n+1) for n in nmax[:-1]], indexing='ij')
    part = sum((g.astype(float)*qi)**2 for g, qi in zip(grids, q[:-1])).ravel()
    part = part[part < r2+dr2]
    lo = np.clip(r2 - part, 0, None); hi = r2 + dr2 - part
    ql = q[-1]
    # count integers n with lo < (n ql)^2 < hi
    a = np.sqrt(lo)/ql; b = np.sqrt(hi)/ql
    # n in (a, b) and (-b, -a), plus n=0 when lo==0 (lo<0 case)
    cnt_pos = np.floor(np.nextafter(b, 0)) - np.floor(a)
    cnt_pos = np.clip(cnt_pos, 0, None)
    zero = (r2 - part < 0) & (hi > 0)       # n=0 allowed when part itself in shell
    total = 2*cnt_pos.sum() + zero.sum()
    return int(total)

rng = np.random.default_rng(1)
print("C1: brute-force lattice-point counts in the BP shell vs eq. (2.24) volume estimate")
for J in [3, 4, 5]:
    for trial in range(3):
        q = np.sqrt(np.array([2, 3, 5, 7, 11][:J], float))*rng.uniform(0.8, 1.2, J)
        r2 = {3: 4.0e4, 4: 3.0e3, 5: 6.0e2}[J]
        dr2 = {3: 40.0, 4: 3.0, 5: 0.6}[J]
        N = shell_count(q, r2, dr2)
        r = np.sqrt(r2); dr = dr2/(2*r)       # Delta r = Delta(r^2)/(2r)  (2.23 with r^2 = 2|lam|)
        Nexp = float(omega(J-1))*r**(J-1)*dr/np.prod(q)
        print(f"   J={J} q={np.round(q,3)}  count={N}  expected={Nexp:.1f}  ratio={N/Nexp:.3f}")

# ---------------- C2: DESI exponential quintessence 2504.04226 ----------------
lam0, up, dn = 0.698, 0.173, 0.202
print("\nC2: lambda = 0.698 +0.173 -0.202 (2504.04226, flat, DESI DR2+CMB+DESY5)")
for name, b in [("TCC sqrt(2/3)", sqrt(2/3)), ("strong sqrt(2)", sqrt(2)), ("OOSV-acceleration sqrt(2)", sqrt(2))][:2]:
    z_lin = (b - lam0)/up
    # alternative: likelihood ~ Gaussian in y = lambda^2 (observables depend on lambda^2 at leading order)
    y0, yu = lam0**2, (lam0+up)**2
    z_y = (b**2 - y0)/(yu - y0)
    print(f"   {name}: Gaussian-in-lambda (upper sigma) z = {z_lin:.2f};  Gaussian-in-lambda^2 z = {z_y:.2f}")
print(f"   lambda=0 with lower sigma (linear): z = {lam0/dn:.2f}  (paper quotes ~3.5 sigma from posterior, 99.95% CL)")
y_lo = (lam0-dn)**2
print(f"   lambda=0 if Gaussian in lambda^2 (lower sigma): z = {(lam0**2)/(lam0**2-y_lo):.2f}  -> this parametrisation is NOT what the posterior shows")

# ---------------- C3: dark dimension 2205.12293 ----------------
hbarc_eVm = 1.973269804e-7
L14 = 2.31e-3                          # eV, Lambda^{1/4} as quoted in 2205.12293 (text after 3.2)
lenL = hbarc_eVm/L14
print(f"\nC3: Lambda^(-1/4) = {lenL*1e6:.1f} micron (paper quotes ~88 micron)")
Mred, Mnon = 2.435e18, 1.2209e19
for lam in [1e-1, 1e-3]:
    l = lam*lenL; m = L14/lam*1e-9      # GeV
    for lab, MP in [("reduced", Mred), ("non-reduced", Mnon)]:
        print(f"   lambda={lam:g}: l={l*1e6:.3f} micron, m=1/l={m*1e9*1e3:.3g} meV, Mhat=m^(1/3) M^(2/3) [{lab}] = {m**(1/3)*MP**(2/3):.3g} GeV;"
              f" with m=2pi/l: {(2*pi*m)**(1/3)*MP**(2/3):.3g} GeV")
