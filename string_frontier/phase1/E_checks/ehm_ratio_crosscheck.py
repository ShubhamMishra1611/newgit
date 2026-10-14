"""Independent float64 cross-check of a21/a20 for the planar CHR family (see ehm_test_planar_family.py):
partial sums to N=4e6 with explicit power-law tail correction."""
import numpy as np
from scipy.special import gammaln, digamma
for lam in [0.25, 0.5, 0.9, 1.0, 1.1, 2.0, 4.0]:
    a = (1+lam)/2; b = (1+3*lam)/2
    n = np.arange(1, 4_000_001, dtype=np.float64)
    logR0 = gammaln(a+n-1) + gammaln(b) - gammaln(a) - gammaln(b+n-1)
    R0 = np.exp(logR0); mu = (n+lam-1)/lam
    T1 = R0*lam*(digamma(a+n-1) - digamma(a))
    t20 = R0/mu**3; t21 = T1/mu**2
    # tail estimate: fit last-decade terms to C n^-p (log n)^q is overkill; use integral of last term * n/(p-1)
    def tot(x, p):
        return x.sum() + x[-1]*n[-1]/(p-1)
    a20 = tot(t20, lam+3); a21 = tot(t21, lam+2)
    a20_half = t20[:2_000_000].sum() + t20[1_999_999]*2e6/(lam+2); a21_half = t21[:2_000_000].sum() + t21[1_999_999]*2e6/(lam+1)
    print("lam=%-5s a20=%.10f a21=%.10f ratio=%.8f   (N=2e6 ratio=%.8f)" % (lam, a20, a21, a21/a20, a21_half/a20_half))
