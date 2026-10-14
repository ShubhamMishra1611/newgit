"""RT3 independent BP (Bousso-Polchinski) discretuum tools.

Model (BP hep-th/0004134 eq. 2.21): lambda = lambda_bare + 1/2 sum q_i^2 n_i^2.
Write X = sum c_i n_i^2 with c_i = q_i^2, shell value S = 2|lambda_bare|; Lambda = (X - S)/2.
Distinct Lambda values <-> orthant points n_i >= 0 (generic c_i).

Written from scratch by RT3 (did not import or copy T3's bp_lib).
Enumeration: split coordinates into two blocks, enumerate each block's orthant
points with partial sum <= Xmax by recursive numpy expansion, sort block B, and
count/list pairs with a+b in a window via searchsorted.  Float64 sums; anything
used as a 'minimum' is re-checked with exact Fractions.
"""
from fractions import Fraction
import math
import numpy as np


def block_sums(c, Xmax):
    """All orthant vectors in len(c) dims with sum c_i n_i^2 <= Xmax.
    Returns (sums float64 unsorted, list-of-int-arrays n) -- n stored as int16 matrix."""
    sums = np.zeros(1)
    ns = np.zeros((1, 0), dtype=np.int16)
    for ci in c:
        nmax = int(math.isqrt(int(Xmax / ci))) + 1
        new_s, new_n = [], []
        for k in range(nmax + 1):
            add = ci * k * k
            if add > Xmax:
                break
            m = sums + add <= Xmax
            if not m.any():
                break
            new_s.append(sums[m] + add)
            new_n.append(np.hstack([ns[m], np.full((m.sum(), 1), k, dtype=np.int16)]))
        sums = np.concatenate(new_s)
        ns = np.vstack(new_n)
    return sums, ns


class Shell:
    def __init__(self, c, Xmax, split=None):
        self.c = np.asarray(c, dtype=np.float64)
        J = len(self.c)
        h = J // 2 if split is None else split
        self.ca, self.cb = self.c[:h], self.c[h:]
        a, na = block_sums(self.ca, Xmax)
        b, nb = block_sums(self.cb, Xmax)
        oa = np.argsort(a)
        ob = np.argsort(b)
        self.a, self.na = a[oa], na[oa]
        self.b, self.nb = b[ob], nb[ob]

    def count(self, lo, hi):
        a = self.a[self.a < hi]
        return int((np.searchsorted(self.b, hi - a) - np.searchsorted(self.b, lo - a)).sum())

    def listing(self, lo, hi, with_n=False):
        a = self.a[self.a < hi]
        s = np.searchsorted(self.b, lo - a)
        e = np.searchsorted(self.b, hi - a)
        cnt = e - s
        keep = np.nonzero(cnt)[0]
        ia = np.repeat(keep, cnt[keep])
        st = np.repeat(s[keep], cnt[keep])
        off = np.arange(ia.size) - np.repeat(np.cumsum(cnt[keep]) - cnt[keep], cnt[keep])
        ib = st + off
        X = self.a[ia] + self.b[ib]
        o = np.argsort(X)
        if with_n:
            return X[o], ia[o], ib[o]
        return X[o]

    def nvec(self, ia, ib):
        return [int(x) for x in self.na[ia]] + [int(x) for x in self.nb[ib]]


def exact_X(c, n):
    return sum(Fraction(float(ci)) * (k * k) for ci, k in zip(c, n))


def esym(x):
    e = np.zeros(len(x) + 1)
    e[0] = 1.0
    for v in x:
        e[1:] = e[1:] + v * e[:-1]
    return e


def smooth_count(c, X):
    """face-corrected smooth orthant count  2^-J sum_k V_k e_k(1/sqrt c) X^{k/2}"""
    w = 1 / np.sqrt(np.asarray(c))
    e = esym(w)
    J = len(c)
    X = np.asarray(X, dtype=float)
    out = 0.0
    for k in range(J + 1):
        Vk = math.exp((k / 2) * math.log(math.pi) - math.lgamma(k / 2 + 1))
        out = out + Vk * e[k] * X ** (k / 2)
    return out / 2 ** J


def smooth_rho(c, S, Lam=0.0):
    """d/dLambda of smooth_count at X = S + 2 Lambda"""
    w = 1 / np.sqrt(np.asarray(c))
    e = esym(w)
    J = len(c)
    X = S + 2 * Lam
    tot = 0.0
    for k in range(1, J + 1):
        om = 2 * math.pi ** (k / 2) / math.gamma(k / 2)
        tot += om * e[k] * X ** ((k - 2) / 2)
    return tot / 2 ** J


def bp_rho(c, S):
    """BP eq. 2.26 inverse with D = 2^J (density of distinct values per unit lambda)."""
    J = len(c)
    om = 2 * math.pi ** (J / 2) / math.gamma(J / 2)
    return om * S ** (J / 2 - 1) / (2 ** J * np.prod(np.sqrt(c)))
