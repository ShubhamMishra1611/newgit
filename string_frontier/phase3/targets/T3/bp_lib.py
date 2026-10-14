"""Shared helpers for the Bousso-Polchinski (BP) discretuum computations (T3 / A1).

BP, hep-th/0004134, eq. (2.21):  lambda = lambda_bare + (1/2) sum_i n_i^2 q_i^2,  lambda_bare < 0.
Write R^2 = 2|lambda_bare| and Q(n) = sum_i q_i^2 n_i^2, so  Lambda = (Q - R^2)/2.

Because n_i -> -n_i leaves Lambda unchanged (BP's degeneracy D = 2^J, text after (2.29)),
we enumerate the orthant n_i >= 0: each orthant point is one DISTINCT value of Lambda
when the q_i^2 are rationally independent (true for generic float64 charges in our range).

Meet-in-the-middle (MITM): split n = (n^(1), n^(2)); enumerate all half-vectors with
Q_h <= Qmax, sort one half, and for each value of the other half binary-search the
partner(s).  This is an EXACT enumeration of the shell, not sampling.

Charges are stored as float64 q_i^2.  Those float64 numbers are taken as the *exact*
definition of the model (dyadic rationals), so exact recomputation of any Lambda is
possible with Python Fractions; float sums are only used for search, and the smallest
|Lambda| values are re-verified exactly.
"""
from fractions import Fraction
import math
import numpy as np


def half_enumerate(q2, Qmax):
    """All orthant vectors n >= 0 (len(q2) coords) with sum q2_i n_i^2 <= Qmax.
    Returns (vals sorted, codes) where codes is a mixed-radix int64 encoding of n."""
    q2 = np.asarray(q2, dtype=np.float64).astype(np.longdouble)   # exact embedding of float64 charges
    Qmax = np.longdouble(Qmax)
    radices = [int(math.floor(math.sqrt(float(Qmax) / float(c)))) + 1 for c in q2]
    vals = np.zeros(1, dtype=np.longdouble)
    codes = np.zeros(1, dtype=np.int64)
    mult = 1
    for c, rad in zip(q2, radices):
        order = np.argsort(vals, kind="stable")
        vals = vals[order]
        codes = codes[order]
        parts_v, parts_c = [], []
        for n in range(rad):
            add = c * np.longdouble(n * n)
            k = np.searchsorted(vals, Qmax - add, side="right")
            if k == 0:
                break
            parts_v.append(vals[:k] + add)
            parts_c.append(codes[:k] + n * mult)
        vals = np.concatenate(parts_v)
        codes = np.concatenate(parts_c)
        mult *= rad
        if mult > 2 ** 62:
            raise OverflowError("code overflow")
    order = np.argsort(vals, kind="stable")
    return vals[order], codes[order], radices


def decode(code, radices):
    n = []
    for r in radices:
        n.append(int(code % r))
        code //= r
    return n


def exact_Q(n, q2):
    return sum(Fraction(float(c)) * (k * k) for k, c in zip(n, q2))


def split_q(q2):
    J = len(q2)
    J1 = J // 2
    return np.asarray(q2[:J1]), np.asarray(q2[J1:])


class BPModel:
    """Enumerate the BP orthant lattice near the shell Q = R^2."""

    def __init__(self, q2, R2, Wq):
        """q2: array of q_i^2; R2 = 2|lambda_bare|; Wq: half-width of the Q-window kept
        (we need half vectors up to Qmax = R2 + Wq)."""
        self.q2 = np.asarray(q2, dtype=np.float64)
        self.R2 = np.longdouble(R2)
        self.Wq = np.longdouble(Wq)
        qa, qb = split_q(self.q2)
        Qmax = self.R2 + self.Wq
        self.A, self.Acode, self.Arad = half_enumerate(qa, Qmax)
        self.B, self.Bcode, self.Brad = half_enumerate(qb, Qmax)
        self.qa, self.qb = qa, qb

    def count_Q_window(self, lo, hi):
        """# orthant points with lo <= Q < hi  (exact up to float rounding at edges)."""
        A, B = self.A, self.B
        msk = A < hi
        a = A[msk]
        c = np.searchsorted(B, hi - a, side="left") - np.searchsorted(B, lo - a, side="left")
        return int(c.sum())

    def counts_many(self, edges):
        """counts in [edges[k], edges[k+1]) for a sorted array of Q edges (vectorized per edge)."""
        A, B = self.A, self.B
        out = []
        for lo, hi in zip(edges[:-1], edges[1:]):
            out.append(self.count_Q_window(lo, hi))
        return np.array(out)

    def list_Q_window(self, lo, hi, max_points=5_000_000):
        """Return sorted Q values of all orthant points with lo <= Q < hi (plus codes)."""
        A, B = self.A, self.B
        msk = A < hi
        idxA = np.nonzero(msk)[0]
        a = A[msk]
        s = np.searchsorted(B, lo - a, side="left")
        e = np.searchsorted(B, hi - a, side="left")
        cnt = e - s
        tot = int(cnt.sum())
        if tot > max_points:
            raise MemoryError(f"{tot} points in window")
        keep = cnt > 0
        idxA, s, cnt = idxA[keep], s[keep], cnt[keep]
        rep_a = np.repeat(idxA, cnt)
        # offsets within each block
        starts = np.repeat(s, cnt)
        offs = np.arange(tot) - np.repeat(np.cumsum(cnt) - cnt, cnt)
        idxB = starts + offs
        Q = A[rep_a] + B[idxB]
        order = np.argsort(Q)
        return Q[order], rep_a[order], idxB[order]

    def closest(self, target, k=50):
        """k candidates with smallest |Q - target| (float search), then exact re-evaluation.
        Returns list of (|Q-target| exact as float, n-vector) sorted exactly."""
        A, B = self.A, self.B
        a = A[A <= target + self.Wq]
        idx = np.searchsorted(B, target - a)
        best = np.full(a.shape, np.inf)
        bidx = np.zeros(a.shape, dtype=np.int64)
        for off in (-1, 0):
            j = np.clip(idx + off, 0, len(B) - 1)
            d = np.abs(a + B[j] - target)
            better = d < best
            best[better] = d[better]
            bidx[better] = j[better]
        # also the second neighbours so that k candidates are not all from distinct a only
        kk = min(k, len(best))
        top = np.argpartition(best, kk - 1)[:kk]
        cands = []
        tgt = Fraction(float(target)) if not isinstance(target, Fraction) else target
        for t in top:
            na = decode(int(self.Acode[t]), self.Arad)
            nb = decode(int(self.Bcode[bidx[t]]), self.Brad)
            n = na + nb
            Qex = exact_Q(n, self.q2)
            cands.append((abs(float(Qex - tgt)), float(Qex - tgt), n, float(best[t])))
        cands.sort(key=lambda x: x[0])
        return cands


def elem_sym(x):
    """elementary symmetric polynomials e_0..e_J of the entries of x."""
    e = np.zeros(len(x) + 1)
    e[0] = 1.0
    for v in x:
        e[1:] = e[1:] + v * e[:-1]
    return e


def smooth_orthant_count(q2, X):
    """Smooth (Euler-Maclaurin type) prediction for #{n in Z_{>=0}^J : Q(n) <= X}:
       N_orth(X) ~ 2^{-J} sum_{U subset [J]} Vol_U(X),  Vol_U = V_{|U|} X^{|U|/2} / prod_{i in U} q_i
    (derived in the report: Moebius inversion of the full-lattice counts on each coordinate face).
    """
    q = np.sqrt(np.asarray(q2))
    J = len(q)
    e = elem_sym(1.0 / q)
    X = np.asarray(X, dtype=np.float64)
    tot = np.zeros_like(X)
    for k in range(0, J + 1):
        Vk = math.pi ** (k / 2) / math.gamma(k / 2 + 1)
        tot = tot + Vk * e[k] * X ** (k / 2)
    return tot / 2 ** J


def smooth_orthant_density_Lambda(q2, R2, Lam=0.0):
    """d N_orth / d Lambda at Lambda (X = R2 + 2 Lambda)."""
    q = np.sqrt(np.asarray(q2))
    J = len(q)
    e = elem_sym(1.0 / q)
    X = R2 + 2 * Lam
    tot = 0.0
    for k in range(1, J + 1):
        omega = 2 * math.pi ** (k / 2) / math.gamma(k / 2)   # area of unit S^{k-1}
        tot += omega * e[k] * X ** ((k - 2) / 2)
    return tot / 2 ** J


def bp_leading_density(q2, R2):
    """BP (2.26) with D = 2^J: 1/Delta lambda = omega_{J-1} |2 lambda_bare|^{J/2-1} / (D prod q)."""
    q = np.sqrt(np.asarray(q2))
    J = len(q)
    omega = 2 * math.pi ** (J / 2) / math.gamma(J / 2)
    return omega * R2 ** (J / 2 - 1) / (2 ** J * np.prod(q))
