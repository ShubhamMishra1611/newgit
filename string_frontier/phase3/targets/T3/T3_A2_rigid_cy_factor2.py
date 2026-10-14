"""T3 (A2): the factor of 2 between exact and continuum flux-vacuum counts in the
rigid-CY / "T^2 with fixed dilaton" toy.

Run:  python3 -I T3_A2_rigid_cy_factor2.py > T3_A2_rigid_cy_factor2.out

Sources read in this session (texts saved in sources/):
  * Ashok-Douglas hep-th/0307049, sec. 3.2 "Example of T^2":
      W = A Z + B, A,B in Z+iZ; DW=0 <-> Zbar = -B/A; tadpole Im(A* B) = L.
      "use SL(2,Z) invariance to set a2 = 0 ... bring b1 into 0 <= b1 < a1 ...
       We can furthermore take a1 > 0, taking into account a1 < 0 by multiplying
       this result by 2. ... Nvac(L) = 2 sigma(L)"   (no eq. number printed)
      and the continuum volume eq. (3.11):
       int d^2A d^2B delta(L - Im A*B) = pi L int d^2rho/(Im rho)^2  -> "the volume is
       pi^2 L/3, which agrees with the L derivative of the previous computation".
      Also (2.23): vol(F) = int_F d^2tau / (4 (Im tau)^2) = pi/12.
  * Denef-Douglas hep-th/0404116, sec. 5.1 "Flux vacua on rigid Calabi-Yau":
      same toy, tadpole (5.4) Im A*B <= L, (5.5) taubar = -B/A; "there are 240 vacua
      at tau = 2i" (L=150); "a disc of sufficiently large area A will contain
      approximately 2 pi A L^2 vacua" (fig. 7, discs around tau = 2i, L = 150).

What this script does
  1. Orbit count modulo the FULL SL(2,Z) (Hermite normal form) = sum sigma_1  (Phase 1).
  2. Brute-force count of all flux pairs (A,B) whose vacuum tau lies in the standard
     fundamental domain F of PSL(2,Z) (half-open boundary convention).  This is what
     the continuum density integrated over F counts.  It equals 2*sum sigma_1 + O(L)
     because N and -N (the central element -1 of SL(2,Z)) give the same tau.
  3. Reproduces DD's "240 vacua at tau = 2i" (L=150) and shows it counts N and -N
     separately (orbits mod full SL(2,Z): 120).
  4. Reproduces DD fig. 7: counts in discs around tau=2i vs 2 pi A L^2 with A the
     proper area in the metric g = 1/(4 y^2) of AD (2.23).
  5. Numerically checks int_F d^2tau/y^2 = pi/3 and the continuum total pi^2 L^2/6.
"""
import math
import numpy as np
from scipy import integrate

print("=" * 78)
print("1. Exact SL(2,Z)-orbit count (HNF: a2=0, a1>0, 0<=b1<a1, a1*b2=n) = sum sigma_1")
Lmax = 3000
sig = np.zeros(Lmax + 1, dtype=np.int64)
for d in range(1, Lmax + 1):
    sig[d::d] += d
cum = np.cumsum(sig)
for L in (50, 150, 500, 1000, 3000):
    s = int(cum[L])
    print(f"  L={L:5d}: sum sigma = {s:9d}; /(pi^2L^2/12) = {s/(math.pi**2*L*L/12):.4f};"
          f"  AD convention 2*sum sigma = {2*s:9d}; /(pi^2L^2/6) = {2*s/(math.pi**2*L*L/6):.4f}")

print("=" * 78)
print("2. Brute force: all flux pairs (A,B) in Z[i]^2 with 0 < Im(conj(A)B) <= L and")
print("   vacuum tau = -conj(B/A) inside the half-open fundamental domain F.")


def in_F(t, eps=1e-12):
    """half-open standard fundamental domain: -1/2 <= x < 1/2, |t|>=1, and on the arc
    |t|=1 only x <= 0 is kept."""
    x, y = t.real, t.imag
    if y <= 0:
        return False
    if x < -0.5 - eps or x >= 0.5 - eps:
        return False
    r2 = x * x + y * y
    if r2 < 1 - eps:
        return False
    if abs(r2 - 1) <= eps and x > eps:
        return False
    return True


def brute_count(L):
    """count (A,B) with tau(A,B) in F, 0<det<=L.  For tau in F: Im tau >= sqrt3/2 and
    det = |A|^2 Im tau, so |A|^2 <= 2L/sqrt3; |B| = |A||tau| with |Re tau|<=1/2,
    Im tau <= det/|A|^2 <= L."""
    count = 0
    at_i = 0
    at_rho = 0
    amax = int(math.isqrt(int(2 * L / math.sqrt(3)) + 1)) + 1
    for a1 in range(-amax, amax + 1):
        for a2 in range(-amax, amax + 1):
            nA = a1 * a1 + a2 * a2
            if nA == 0 or nA > 2 * L / math.sqrt(3) + 1e-9:
                continue
            A = complex(a1, a2)
            # tau = -conj(B/A)  <=>  B = -A*conj(tau);  |B|^2 = |A|^2 |tau|^2 <= |A|^2 (1/4 + (L/|A|^2)^2)
            bmax = int(math.sqrt(nA * (0.25 + (L / nA) ** 2))) + 1
            for b1 in range(-bmax, bmax + 1):
                for b2 in range(-bmax, bmax + 1):
                    det = a1 * b2 - a2 * b1  # Im(conj(A) B)
                    if det <= 0 or det > L:
                        continue
                    t = -(complex(b1, b2) / A).conjugate()
                    if in_F(t):
                        count += 1
                        if abs(t - 1j) < 1e-9:
                            at_i += 1
                        if abs(t - complex(-0.5, math.sqrt(3) / 2)) < 1e-9:
                            at_rho += 1
    return count, at_i, at_rho


for L in (10, 20, 40, 80):
    c, ci, cr = brute_count(L)
    s = int(cum[L])
    # orbits through tau=i (stabilizer order 4 in SL(2,Z)) and tau=rho (order 6):
    # each such orbit has 4 resp. 6 flux representatives with tau exactly at the point,
    # a generic orbit has 2 (N and -N).  So c = 2*s + 2*(#orbits at i) + 4*(#orbits at rho).
    n_orb_i = ci // 4
    n_orb_rho = cr // 6
    pred = 2 * s + 2 * n_orb_i + 4 * n_orb_rho
    print(f"  L={L:3d}: #(A,B) with tau in F = {c:6d}; sum sigma = {s:5d}; 2*sum sigma = {2*s:6d};"
          f" fluxes at tau=i: {ci}, at rho: {cr};  2*sum sigma + elliptic corr = {pred}  match={pred==c}")

print("=" * 78)
print("3. Vacua exactly at tau = 2i for L = 150 (DD sec. 5.1 quotes 240)")
L = 150
nA = sum(1 for x in range(-20, 21) for y in range(-20, 21) if 0 < 2 * (x * x + y * y) <= L)
print(f"  fluxes (A, B = 2i A) with 0 < det = 2|A|^2 <= 150: {nA}   (DD: 240)")
print(f"  modulo the full SL(2,Z) (stabilizer of 2i is {{+1,-1}}): {nA//2}")

# cross-check via the HNF enumeration of Phase 1 (one rep per full-SL(2,Z) orbit)
def reduce_to_F(t):
    for _ in range(10000):
        t = complex(t.real - math.floor(t.real + 0.5), t.imag)
        if abs(t) < 1 - 1e-13:
            t = -1 / t
        else:
            return t
    return t

taus = []
for a1 in range(1, L + 1):
    for b2 in range(1, L // a1 + 1):
        for b1 in range(0, a1):
            taus.append(reduce_to_F(-(complex(b1, b2) / complex(a1, 0)).conjugate()))
taus = np.array(taus)
print(f"  HNF orbit enumeration (L=150): {len(taus)} orbits, at tau=2i: {int(np.sum(np.abs(taus-2j)<1e-9))}")

print("=" * 78)
print("4. DD fig. 7: number of vacua in discs |tau - 2i| < R (coordinate radius), L = 150")
print("   estimate 2 pi A L^2 with A = proper area in g = 1/(4y^2) (AD (2.23)).")


def proper_area_disc(R, y0=2.0):
    # int over disc of dx dy /(4 y^2)
    f = lambda r, th: r / (4.0 * (y0 + r * math.sin(th)) ** 2)
    val, _ = integrate.dblquad(lambda r, th: f(r, th), 0, 2 * math.pi, 0, R)
    return val

d = np.abs(taus - 2j)
print("    R    | estimate 2piAL^2 | orbits mod SL(2,Z) | x2 (N and -N counted) ")
for R in (0.05, 0.10, 0.12, 0.15, 0.20, 0.25, 0.30):
    A = proper_area_disc(R)
    est = 2 * math.pi * A * L * L
    n_orb = int(np.sum(d < R))
    print(f"   {R:.2f}  |   {est:10.1f}     |   {n_orb:8d}         |  {2*n_orb:8d}")

print("=" * 78)
print("5. Continuum integral checks")
# int_F dx dy / y^2  (hyperbolic area) = pi/3
val, _ = integrate.dblquad(lambda y, x: 1.0 / y ** 2, -0.5, 0.5, lambda x: math.sqrt(1 - x * x), lambda x: np.inf)
print(f"  int_F d^2tau / y^2 = {val:.10f}  (pi/3 = {math.pi/3:.10f});  with 1/(4y^2): {val/4:.10f} (pi/12)")
for L in (150, 1000, 3000):
    cont = 2 * math.pi * (val / 4) * L * L  # DD: 2 pi A L^2 integrated over F
    ad = math.pi * L * L / 2 * val           # AD (3.11) integrated: int_0^L pi L' dL' * int d^2rho/y^2
    print(f"  L={L}: DD 2 pi A(F) L^2 = {cont:.1f};  AD (3.11) integrated = {ad:.1f};  pi^2L^2/6 = {math.pi**2*L*L/6:.1f};"
          f"  exact orbits = {int(cum[L])};  2*orbits = {2*int(cum[L])}")
print("\nVERDICT: AD's exact count explicitly multiplies by 2 for a1<0 (i.e. counts N and -N as distinct),")
print("and that is what the continuum density integrated over F reproduces.  Modulo the full SL(2,Z),")
print("which contains -1: (A,B) -> (-A,-B) at fixed tau, the count is sum sigma ~ pi^2 L^2/12.")
