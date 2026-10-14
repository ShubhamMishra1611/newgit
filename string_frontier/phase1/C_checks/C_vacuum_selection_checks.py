"""Front C quick checks (Phase 1). Run: python3 -I C_checks/C_vacuum_selection_checks.py (from phase1/)
Every number printed here is quoted in C_vacuum_selection.md as "verified by computation".
Sources for formulas are given inline (arXiv ID + eq. number as read in this session).
"""
import itertools, math
import numpy as np
from mpmath import mp, mpf, gamma, pi, log10

mp.dps = 30
print("=" * 70)
print("1. Bousso-Polchinski spacing, hep-th/0004134 eq. (2.26)")
# Delta_lambda_min = D * prod(q_i) / (omega_{J-1} |2 lambda_bare|^{J/2-1}),
# omega_{J-1} = 2 pi^{J/2} / Gamma(J/2).  Units: BP's (sqrt2 kappa4)^4 lambda_bare ~ 1.
def omega(Jm1):
    J = Jm1 + 1
    return 2 * pi ** (mpf(J) / 2) / gamma(mpf(J) / 2)

for J in (4, 6, 20, 100, 300):
    lam_bare = mpf(1)
    target = mpf(10) ** -120
    for label, logD in (("D=1", 0), ("D=2^J", J * math.log10(2))):
        # uniform q: q^J = target * omega * |2lam|^{J/2-1} / D
        logq = (log10(target) + log10(omega(J - 1)) + (mpf(J) / 2 - 1) * log10(2 * lam_bare) - logD) / J
        print(f"  J={J:4d} {label:6s}: uniform charge needed q = 10^{float(logq):7.3f}  (q^(1/2)=10^{float(logq)/2:6.3f})")

print("\n   Brute-force check of the shell-volume estimate (2.24) at small J")
rng = np.random.default_rng(1)
for J, R2, width in ((3, 400.0, 2.0), (4, 300.0, 2.0), (5, 150.0, 2.0)):
    q = 1.0 + rng.random(J)  # incommensurate O(1) charges
    nmax = [int(math.sqrt(R2 + width) / qi) + 1 for qi in q]
    grids = np.meshgrid(*[np.arange(-m, m + 1) for m in nmax], indexing="ij")
    s = sum((g.astype(float) * qi) ** 2 for g, qi in zip(grids, q))
    count = int(np.count_nonzero((s >= R2) & (s < R2 + width)))
    # shell in r: r in [sqrt(R2), sqrt(R2+width)]; volume = omega_{J-1}/J * (r2^J - r1^J)
    r1, r2 = math.sqrt(R2), math.sqrt(R2 + width)
    vol = float(omega(J - 1)) / J * (r2 ** J - r1 ** J)
    pred = vol / np.prod(q)
    print(f"  J={J}: lattice points in shell = {count:6d}, volume estimate = {pred:9.1f}, ratio = {count/pred:.3f}")

print("=" * 70)
print("2. Rigid CY toy (Denef-Douglas hep-th/0404116 sec. 5.1): W = A tau + B, A,B in Z[i],")
print("   tadpole 0 < Im(conj(A) B) <= L; SL(2,Z) reps a2=0, 0<=b1<a1, a1*b2 = n  ->  #orbits(n) = sigma_1(n)")
Lmax = 3000
sig = np.zeros(Lmax + 1, dtype=np.int64)
for d in range(1, Lmax + 1):
    sig[d::d] += d
cum = np.cumsum(sig)
for L in (50, 150, 500, 1000, 3000):
    exact = int(cum[L])
    print(f"  L={L:5d}: exact N(L)={exact:9d}  pi^2 L^2/12={math.pi**2*L*L/12:12.1f}  ratio={exact/(math.pi**2*L*L/12):.4f}"
          f"   continuum AD(1.5)/DD(3.15) eval = pi^2 L^2/6 -> ratio {exact/(math.pi**2*L*L/6):.4f}")

# brute-force check that the orbit count equals sum sigma_1 for small L, via explicit reduction
def reduce_to_fd(t):
    # map tau to the standard PSL(2,Z) fundamental domain
    for _ in range(1000):
        t = complex(t.real - math.floor(t.real + 0.5), t.imag)
        if abs(t) < 1 - 1e-13:
            t = -1 / t
        else:
            return t
    return t

L = 150
taus = []
for a1 in range(1, L + 1):
    for b2 in range(1, L // a1 + 1):
        for b1 in range(0, a1):
            A = complex(a1, 0)
            B = complex(b1, b2)
            tau = -(B / A).conjugate()  # DD (5.5): conj(tau) = -B/A
            taus.append(reduce_to_fd(tau))
taus = np.array(taus)
print(f"  L={L}: enumerated {len(taus)} vacua (should equal sum sigma_1 = {int(cum[L])})")
print(f"  all Im(tau) > 0 after reduction: {bool(np.all(taus.imag > 0))}")
for y in (1.5, 2.0, 4.0):
    frac = np.mean(taus.imag > y)
    print(f"   fraction with Im tau > {y}: {frac:.4f}  vs uniform hyperbolic prediction 3/(pi y) = {3/(math.pi*y):.4f}")
n_i = int(np.sum(np.abs(taus - 1j) < 1e-9)); n_rho = int(np.sum(np.abs(taus - complex(0.5, math.sqrt(3)/2)) < 1e-9)
                                                     + np.sum(np.abs(taus - complex(-0.5, math.sqrt(3)/2)) < 1e-9))
print(f"   vacua exactly at tau=i: {n_i}, at tau=rho: {n_rho} (out of {len(taus)})")

print("=" * 70)
print("3. Exponential quintessence vs swampland slopes (d=4)")
lam, up, dn = 0.698, 0.173, 0.202   # arXiv:2504.04226 abstract, DESI DR2 combo, Omega_K=0
bounds = {"TCC average  2/sqrt((d-1)(d-2)) [1909.11063]": 2 / math.sqrt(6),
          "strong dS    2/sqrt(d-2)       [2101.11617]": 2 / math.sqrt(2)}
for k, b in bounds.items():
    print(f"  {k}: {b:.4f}; DESI fit lambda=0.698(+0.173): (bound - lambda)/sigma_up = {(b-lam)/up:+.2f}")
for l in (lam, 2 / math.sqrt(6), math.sqrt(2)):
    print(f"   attractor w_inf = -1 + lambda^2/3 (1806.08362 eq.1.3) at lambda={l:.3f}: w={-1+l*l/3:+.3f}")

print("=" * 70)
print("4. Dark dimension numerics (2205.12293 eqs. (2.1),(3.2)-(3.5))")
Mpl = 2.435e18          # reduced Planck mass, GeV
L14 = 2.31e-12          # Lambda^{1/4} = 2.31 meV (2205.12293 text near (3.2))
hbarc_m = 1.97327e-16   # GeV*m
print(f"  Lambda^-1/4 = {hbarc_m/L14*1e6:.1f} micron  (paper: ~88 micron)")
print(f"  Lambda/Mpl^4 = {(L14/Mpl)**4:.2e}")
for lamDD in (1e-1, 1e-2, 1e-3):
    m = L14 / lamDD
    Mhat = m ** (1 / 3) * Mpl ** (2 / 3)   # eq. (2.1) with n=1
    print(f"  lambda={lamDD:.0e}: m={m*1e9:.3g} eV, l=1/m={hbarc_m/m*1e6:.3g} micron, species scale Mhat={Mhat:.2e} GeV")

print("=" * 70)
print("5. Size of the Ashok-Douglas prefactor (2 pi L)^K / K!  (hep-th/0307049 eq. (1.5)), log10")
from math import lgamma, log, pi as PI
for K, L in ((100, 100), (300, 300), (300, 1000), (300, 3000), (600, 1000)):
    val = (K * log(2 * PI * L) - lgamma(K + 1)) / log(10)
    print(f"  K=b3={K:4d}, L={L:5d}: log10[(2 pi L)^K/K!] = {val:8.1f}")
