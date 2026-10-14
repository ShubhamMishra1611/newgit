"""RT3 independent brute force for T3 claim 6 (AD hep-th/0307049 sec 3.2 / DD hep-th/0404116 sec 5.1).

Toy: W = A tau + B, A=a1+i a2, B=b1+i b2 in Z[i];  DW=0 <-> conj(tau) = -B/A, i.e. tau = -conj(B)/conj(A);
tadpole L' = Im(conj(A) B) = a1 b2 - a2 b1, 0 < L' <= L.
Half-open fundamental domain F: -1/2 <= Re tau < 1/2, |tau| >= 1, and Re tau <= 0 if |tau| = 1.
All arithmetic exact (Fractions).
"""
import math
from fractions import Fraction as Fr
from collections import Counter


def sigma(n):
    return sum(d for d in range(1, n + 1) if n % d == 0)


def tau_of(a1, a2, b1, b2):
    # tau = -conj(B)/conj(A) = -(b1 - i b2)/(a1 - i a2) = -(b1 - i b2)(a1 + i a2)/|A|^2
    nA = a1 * a1 + a2 * a2
    re = Fr(-(b1 * a1 + b2 * a2), nA)
    im = Fr(-(b1 * a2 - b2 * a1), nA)
    return re, im


def in_F(re, im):
    if im <= 0:
        return False
    if not (Fr(-1, 2) <= re < Fr(1, 2)):
        return False
    r2 = re * re + im * im
    if r2 < 1:
        return False
    if r2 == 1 and re > 0:
        return False
    return True


def fluxes_in_F(L):
    """all fluxes with 0<L'<=L and tau in F; returns list of (flux, tau, L')"""
    out = []
    # bounds: Im tau >= sqrt3/2 => |A|^2 <= 2L/sqrt3 ; |B|^2 = |tau|^2 |A|^2 <= (1/4 + L^2/|A|^4)|A|^2
    Amax2 = int(2 * L / math.sqrt(3)) + 1
    am = int(math.isqrt(Amax2)) + 1
    for a1 in range(-am, am + 1):
        for a2 in range(-am, am + 1):
            nA = a1 * a1 + a2 * a2
            if nA == 0 or nA > Amax2:
                continue
            Bmax2 = nA / 4 + L * L / nA + 1
            bm = int(math.isqrt(int(Bmax2))) + 1
            for b1 in range(-bm, bm + 1):
                for b2 in range(-bm, bm + 1):
                    Lp = a1 * b2 - a2 * b1
                    if Lp <= 0 or Lp > L:
                        continue
                    re, im = tau_of(a1, a2, b1, b2)
                    if in_F(re, im):
                        out.append(((a1, a2, b1, b2), (re, im), Lp))
    return out


print("T3 claim: #(fluxes with tau in F, L'<=L) = 2*sum sigma + 2*#(orbits at i) + 4*#(orbits at rho)")
for L in (10, 20, 40, 60, 80):
    fl = fluxes_in_F(L)
    n = len(fl)
    ssig = sum(sigma(k) for k in range(1, L + 1))
    at_i = sum(1 for f in fl if f[1] == (0, 1))
    at_rho = sum(1 for f in fl if f[1][0] == Fr(-1, 2) and f[1][1] * f[1][1] == Fr(3, 4))
    # orbits at i: SL(2,Z) stabiliser of i has order 4 and acts freely on fluxes -> orbits = at_i/4
    pred = 2 * ssig + 2 * (at_i // 4) + 4 * (at_rho // 6)
    print(f"  L={L:3d}: brute #fluxes={n:6d};  2*sum_sigma={2*ssig:6d};  fluxes at tau=i: {at_i} (orbits {at_i//4});"
          f" at rho: {at_rho};  2sum_sigma+2#orb(i)+4#orb(rho) = {pred}  -> {'EXACT' if pred == n else 'MISMATCH'};"
          f"  pi^2 L^2/6 = {math.pi**2*L*L/6:.0f}")

# DD's 240 at tau = 2i, L = 150
L = 150
cnt_all = 0
for a1 in range(-13, 14):
    for a2 in range(-13, 14):
        if 0 < 2 * (a1 * a1 + a2 * a2) <= L:
            cnt_all += 1
print(f"\nDD tau=2i, L=150: #A with 0<2|A|^2<=150 (B=2iA) = {cnt_all}; SL(2,Z)-orbits (only +-1 fix 2i) = {cnt_all//2}")
# DD's stated recipe: a2=0, 0<=b1<a1, a1 b2 <= L ; map tau to F and count those landing at 2i
def to_F(re, im):
    # standard reduction
    while True:
        # shift
        k = math.floor(re + Fr(1, 2))
        re = re - k
        r2 = re * re + im * im
        if r2 < 1:
            re, im = -re / r2, im / r2
            continue
        if re == Fr(1, 2):
            re = Fr(-1, 2)
        if r2 == 1 and re > 0:
            re = -re
        return re, im
landed = Counter()
nrep = 0
for a1 in range(1, L + 1):
    for b2 in range(1, L // a1 + 1):
        for b1 in range(0, a1):
            nrep += 1
            re, im = tau_of(a1, 0, b1, b2)
            landed[to_F(re, im)] += 1
print(f"DD recipe (a2=0, 0<=b1<a1, a1 b2<=150): {nrep} representatives (= sum_{{L'<=150}} sigma(L') = {sum(sigma(k) for k in range(1,151))});"
      f" landing at tau=2i: {landed[(Fr(0), Fr(2))]}")
# DD fig 7: discs of coordinate radius R around 2i, counting all fluxes (N and -N separate)
fl = fluxes_in_F(150)
print(f"all fluxes with tau in F at L=150: {len(fl)} (2*sum sigma = {2*sum(sigma(k) for k in range(1,151))})")
for R in (0.15, 0.20, 0.25, 0.30):
    c = sum(1 for f in fl if (float(f[1][0]) ** 2 + (float(f[1][1]) - 2) ** 2) < R * R)
    # DD estimate 2 pi A L^2 with A = int d^2tau /(4 y^2) over the disc (numerical)
    import numpy as np
    n = 2000
    xs = np.linspace(-R, R, n); ys = np.linspace(2 - R, 2 + R, n)
    X, Y = np.meshgrid(xs, ys)
    m = X ** 2 + (Y - 2) ** 2 < R * R
    A = np.sum(1 / (4 * Y[m] ** 2)) * (xs[1] - xs[0]) * (ys[1] - ys[0])
    print(f"  R={R:.2f}: brute count (N,-N separate) = {c};  half = {c//2};  2 pi A L^2 = {2*math.pi*A*L*L:.0f}")

print("\nboundary check with exact arithmetic (|tau-2i|^2 < R^2 vs <= R^2):")
for R in (Fr(15, 100), Fr(20, 100), Fr(25, 100), Fr(30, 100)):
    lt = sum(1 for f in fl if f[1][0] ** 2 + (f[1][1] - 2) ** 2 < R * R)
    le = sum(1 for f in fl if f[1][0] ** 2 + (f[1][1] - 2) ** 2 <= R * R)
    print(f"  R={float(R):.2f}: strict {lt}, closed {le}")
