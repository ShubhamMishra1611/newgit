"""T3 (A1e): face-corrected density of DISTINCT Lambda values for BP's own large-J examples.
rho(0) = 2^{-J} sum_k omega_{k-1} C(J,k) q^{-k} X^{(k-2)/2},  X = R^2 = 2|lambda_bare|, uniform typical charge q
(our smooth orthant formula, bp_lib.smooth_orthant_density_Lambda, with all q_i = q), vs BP (2.26) with D = 2^J
(the k = J term alone).  BP's own caveat (text after (2.26)): (2.24)-(2.26) break down if q > J^{-1/2}|2 lambda_bare|^{1/2},
i.e. n_max = R/q < sqrt(J).  Run: python3 -I T3_A1e_bp_faces.py > T3_A1e_bp_faces.out
"""
from math import lgamma, log, log10, pi, sqrt

def log10_terms(J, q, X):
    out = []
    for k in range(1, J + 1):
        lw = log(2) + (k / 2) * log(pi) - lgamma(k / 2)
        lC = lgamma(J + 1) - lgamma(k + 1) - lgamma(J - k + 1)
        out.append((lw + lC - k * log(q) + (k - 2) / 2 * log(X) - J * log(2)) / log(10))
    return out

print(" J    log10 q  n_max=R/q  sqrt(J) | log10 rho: BP(2.26),D=2^J | face-corrected | dominant k | gain")
for J, lq in ((100, -1.6), (100, -1.43), (100, -1.73), (300, -0.87), (300, -1.171)):
    q = 10 ** lq
    X = 2.0
    t = log10_terms(J, q, X)
    m = max(t)
    tot = m + log10(sum(10 ** (x - m) for x in t))
    print(f"{J:4d}  {lq:7.3f}  {sqrt(X)/q:8.1f}  {sqrt(J):6.1f} |        {t[-1]:8.2f}            |   {tot:8.2f}     |   {t.index(m)+1:4d}     | 10^{tot-t[-1]:.2f}")
print("(lambda_bare = 1 in BP's units; '# distinct vacua with |lambda| < 1e-120' = 2 rho 1e-120.)")
