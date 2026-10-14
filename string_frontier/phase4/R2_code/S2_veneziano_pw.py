"""
R2 independent check of S2: partial-wave (Gegenbauer) positivity of the bosonic Veneziano amplitude
  A = Gamma(-a(s))Gamma(-a(t))/Gamma(-a(s)-a(t)),  a(x) = 1 + x  (alpha'=1), external tachyons m^2 = -1.
Residue at a(s)=n:  R_n = (a(t)+1)_n / n!  (polynomial of degree n in t).
Kinematics (all outgoing/incoming 2->2, equal masses m^2=-1):  s+t+u = 4m^2 = -4,
  t = -(s - 4 m^2)(1 - x)/2 = -(s+4)(1-x)/2,  x = cos(theta); at level n, s = n-1.
Gegenbauer basis C_j^{lam}(x), lam = (D-3)/2 (standard normalisation, C_1 = 2 lam x).
Monomials -> Gegenbauer via the exact expansion
  x^k = k!/2^k * sum_{m=0}^{floor(k/2)} (k-2m+lam)/(m! (lam)_{k-m+1}) C_{k-2m}^{lam}(x),
which is verified numerically below (independently of the phase-3 script, which I read only afterwards).
"""
import sympy as sp
D, x = sp.symbols('D x')
lam = (D - 3)/2

def mono_to_geg(k):
    out = {}
    for m in range(k//2 + 1):
        j = k - 2*m
        out[j] = out.get(j, 0) + sp.factorial(k)/2**k*(j + lam)/(sp.factorial(m)*sp.rf(lam, k - m + 1))
    return out

# numeric verification of the expansion formula
import mpmath as mp
for kk in range(0, 7):
    for Dv in [4, 7.5, 26]:
        lv = (Dv-3)/2; xv = 0.37
        lhs = xv**kk
        rhs = sum(float(c.subs(D, Dv))*float(mp.gegenbauer(j, lv, xv)) for j, c in mono_to_geg(kk).items())
        assert abs(lhs-rhs) < 1e-12, (kk, Dv, lhs, rhs)
print("monomial->Gegenbauer expansion verified for k<=6")

def level_coeffs(n):
    s = n - 1
    t = -(s + 4)*(1 - x)/2
    R = sp.expand(sp.rf(t + 2, n)/sp.factorial(n))
    P = sp.Poly(R, x)
    c = {}
    for (k,), a in P.terms():
        for j, g in mono_to_geg(k).items():
            c[j] = c.get(j, 0) + a*g
    return {j: sp.factor(sp.simplify(v)) for j, v in c.items()}

def max_D(cs):
    # largest D > 3 such that all c_j(D) >= 0 on (3, Dmax]
    crit = []
    for j, cj in cs.items():
        num, den = sp.fraction(sp.together(cj))
        for r in sp.Poly(num, D).nroots(n=30) if sp.Poly(num, D).degree() > 0 else []:
            if abs(sp.im(r)) < 1e-20 and sp.re(r) > 3.0001:
                crit.append(float(sp.re(r)))
    # test positivity just above 3 and find first sign change
    cand = sorted(set(crit))
    for r in cand:
        # if some c_j becomes negative right after r, r is the bound
        if any(float(cj.subs(D, r + 1e-6)) < -1e-30 for cj in cs.values()):
            return r
    return float('inf')

if __name__ == "__main__":
    for n in range(1, 9):
        cs = level_coeffs(n)
        if n <= 3:
            for j in sorted(cs, reverse=True):
                print(f"  n={n} j={j}: c = {cs[j]}")
        print(f"level n={n}: all c_j >= 0 for 3 < D <= {max_D(cs):.6f}")
    # check positivity at D=26 and D=10 numerically for n<=20, and D=26.0001 failure
    for Dv in [4, 10, 26, 26.01]:
        neg = []
        for n in range(1, 21):
            for j, cj in level_coeffs(n).items():
                v = float(cj.subs(D, Dv))
                if v < -1e-14: neg.append((n, j, v))
        print(f"D={Dv}: negative coefficients up to n=20: {neg[:4]}{' ...' if len(neg)>4 else ''}")
