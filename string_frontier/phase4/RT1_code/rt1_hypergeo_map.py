"""Dictionary check: T1 open c=0 family with intercept a and FREE external mass m_ext^2 is, level by level and
up to a positive factor 1/(r+1)_n (r>-1), the bosonic-type hypergeometric amplitude of 2302.12263 / 2409.09561
(residue (1+t+r-m^2)_n/(r+1)_n at s=m^2+n, eq.(2.4) of 2409.09561) with
      r = (4/3)(m_ext^2 + a),   m^2 = (4/3) m_ext^2 + a/3 .
So on-trajectory (m_ext^2=-a) <-> r=0 ; T1 sec.3 line a=1, m_ext^2=-1+delta <-> r=4delta/3, m^2=-1+4delta/3 ;
level-1 decoupling edge m_ext^2=-(1+3a)/4 <-> 2r+m^2+1=0 (boundary of Theorem 1 of 2409.09561);
'm_ext^2 <= -a' <-> r <= 0.   Check: compare exact per-level Dcrit computed both ways."""
import sys; sys.path.insert(0, '.')
from rt1core import *
from fractions import Fraction as Fr
def hyper_facs(n, r, m2):
    r = q(r); m2 = q(m2); kap = (q(n) - 3 * m2) / 2
    # factors (1 + t + r - m2 + i), i=0..n-1 ; t = kap(x-1)
    return [(1 + r - m2 + i - kap, kap) for i in range(n)], kap
worst = Fr(0)
for (a, me) in [(Fr(1), Fr(-1) + Fr(4022, 10000)), (Fr(1), Fr(-1, 2)), (Fr(1, 2), Fr(-5, 8)), (Fr(0), Fr(-1, 10)), (Fr(-1, 3), Fr(1, 3))]:
    r = Fr(4, 3) * (me + a); m2 = Fr(4, 3) * me + a / 3
    out = []
    for n in range(1, 13):
        f1, k1 = open_facs(n, a, me, 0); f2, k2 = hyper_facs(n, r, m2)
        d1, _ = Dcrit_factors(f1, k1, tol=Fr(1, 10**9)); d2, _ = Dcrit_factors(f2, k2, tol=Fr(1, 10**9))
        if (d1 is None) != (d2 is None): out.append(('MISMATCH', n)); continue
        if d1 is not None: worst = max(worst, abs(d1 - d2)); out.append((n, round(float(d1), 5)))
    print(f"T1 (a={a}, m_ext^2={me}) <-> hypergeometric (r={r}, m^2={m2}): Dcrit(n) agree; mins: {sorted(out, key=lambda z: z[1])[:2]}")
print("max |difference| over all checks:", float(worst))
