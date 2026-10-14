"""RT3 check of T3 claim 2 at large J (J=100, 300): is the face-corrected density right, and
how far is BP eq. (2.26) with D=2^J from the TRUE number of orthant points (distinct values for
incommensurate charges)?  Independent method: exact lattice counting with generating functions.

For equal charges q the orthant count is the coefficient of x^N in f(x)^J, f(x)=sum_{n>=0} x^{n^2}.
Incommensurate charges of the same magnitude have the same orthant-point count up to the smearing
of each integer N into a cloud; so the true distinct-value density per unit lambda near the shell is
   rho_true = [sum_{N in window} coeff_N / window width in N] * (2/q^2)      (lambda = lam_bare + q^2 N/2).
Computed in float64 (coefficients < 1e308) with numpy convolution and binary exponentiation.
Units: BP, lambda_bare = -1 (|2 lambda_bare| = 2), as in T3's part (e) and A1e.
"""
import math
import numpy as np


def polypow(f, J, deg):
    res = np.zeros(deg + 1); res[0] = 1.0
    base = f.copy()
    while J:
        if J & 1:
            res = np.convolve(res, base)[: deg + 1]
        J >>= 1
        if J:
            base = np.convolve(base, base)[: deg + 1]
    return res


def log10_face_density(J, q, S):
    # 2^-J sum_k omega_{k-1} C(J,k) q^-k S^{(k-2)/2}  (per unit Lambda)
    terms = []
    for k in range(1, J + 1):
        lom = math.log(2) + (k / 2) * math.log(math.pi) - math.lgamma(k / 2)
        terms.append(lom + math.lgamma(J + 1) - math.lgamma(k + 1) - math.lgamma(J - k + 1) - k * math.log(q)
                     + (k - 2) / 2 * math.log(S))
    m = max(terms)
    return (m + math.log(sum(math.exp(t - m) for t in terms)) - J * math.log(2)) / math.log(10), int(np.argmax(terms)) + 1


def log10_bp(J, q, S):
    lom = math.log(2) + (J / 2) * math.log(math.pi) - math.lgamma(J / 2)
    return (lom + (J / 2 - 1) * math.log(S) - J * math.log(2) - J * math.log(q)) / math.log(10)


S = 2.0
print(" J   log10 q  n_max  | log10 rho: BP(D=2^J) | face formula (k*) | exact lattice (window) | exact/face | exact/BP")
for J, lq in ((100, -1.6), (100, -1.43), (100, -1.73), (300, -0.87), (300, -1.171), (24, -1.0)):
    q = 10 ** lq
    N0 = S / q ** 2               # shell in units of q^2
    deg = int(N0 + 40)
    nm = int(math.isqrt(deg))
    f = np.zeros(deg + 1); f[[k * k for k in range(nm + 1)]] = 1.0
    F = polypow(f, J, deg)
    # window: +-w integers around N0 (smooth out arithmetic fluctuations of r_J(N))
    w = 8
    lo, hi = int(math.floor(N0)) - w + 1, int(math.floor(N0)) + w
    dens_N = F[lo:hi + 1].sum() / (hi - lo + 1)
    l_exact = math.log10(dens_N) + math.log10(2 / q ** 2)
    l_face, kstar = log10_face_density(J, q, S)
    l_bp = log10_bp(J, q, S)
    print(f"{J:4d} {lq:7.3f} {math.sqrt(N0):6.1f} | {l_bp:10.2f}           | {l_face:8.2f} ({kstar:3d})    | {l_exact:10.2f}             |"
          f" 10^{l_exact-l_face:+.2f}   | 10^{l_exact-l_bp:+.2f}")
print("\nNB: J=300, q=10^-0.87 and 10^-1.171 do not occur in BP (hep-th/0004134); BP's J=100 example uses D 'small' (q=10^-1.6).")
