"""Claim 2: what is the correct positivity condition for a pole BELOW the 2-particle threshold?
Explicit Feynman-rule check in a manifestly unitary Lagrangian (mostly-minus metric):
  L = |d phi|^2 - m^2|phi|^2 - 1/4 F^2 + 1/2 M^2 A^2 - i e A^mu (phi* d_mu phi - phi d_mu phi*)   (Proca, real e)
      + spin-2: g h_{mu nu} d^mu phi* d^nu phi (Fierz-Pauli massive h)          + scalar: -g0 sigma |phi|^2
s-channel phi(p1) phi*(p2) -> phi(p3) phi*(p4).  Writing A ~ Res/(M^2 - s):
   scalar : Res = g0^2
   vector : Res = -e^2 (p1-p2).(p3-p4) = e^2 (t-u)
   spin 2 : Res = g^2 [ (q.Pt.q')^2 - (q.Pt.q)(q'.Pt.q')/(D-1) ] ,  Pt = -eta + k k/M^2
We evaluate these with EXPLICIT complex on-shell momenta in D=5 (p imaginary when M^2<4m^2) and compare with
   Res_j = g_j^2 (s-4m^2)^j * C_j^{nu}(x) / C_j^{nu} leading-normalisation,  x = 1 + 2t/(s-4m^2).
Conclusion printed: c_1 has the sign of (M^2-4m^2) for a healthy (real e) vector -> the correct condition is
rho_j = c_j/kappa^j >= 0, NOT c_j >= 0.  Also: the Veneziano level-1 residue IS exactly this Proca residue.
"""
import sympy as sp
Dd = 5
m2, M2, th = sp.Rational(1, 2), sp.Rational(3, 2), sp.Rational(7, 10)   # M^2 < 4 m^2 = 2 : below threshold
def run(m2, M2, th):
    E = sp.sqrt(M2) / 2
    p = sp.sqrt(E**2 - m2)            # imaginary below threshold
    n1 = [1, 0, 0, 0]; n2 = [sp.cos(th), sp.sin(th), 0, 0]
    def vec(E, p, n): return [E] + [p * c for c in n]
    p1 = vec(E, p, n1); p2 = vec(E, -p, n1); p3 = vec(E, p, n2); p4 = vec(E, -p, n2)
    eta = sp.diag(1, *([-1] * (Dd - 1)))
    dot = lambda a, b: sp.expand(sum(eta[i, i] * a[i] * b[i] for i in range(Dd)))
    sub = lambda a, b: [a[i] - b[i] for i in range(Dd)]
    add = lambda a, b: [a[i] + b[i] for i in range(Dd)]
    s = dot(add(p1, p2), add(p1, p2)); t = dot(sub(p1, p3), sub(p1, p3)); u = dot(sub(p1, p4), sub(p1, p4))
    q = sub(p1, p2); qp = sub(p3, p4); k = add(p1, p2)
    Pt = sp.Matrix(Dd, Dd, lambda i, j: -eta[i, j] + (eta[i, i] * k[i]) * (eta[j, j] * k[j]) / M2)  # lower-index projector
    bil = lambda a, b: sp.expand(sum(a[i] * Pt[i, j] * b[j] for i in range(Dd) for j in range(Dd)))
    resV = sp.simplify(-dot(q, qp))
    resV_proj = sp.simplify(bil(q, qp))          # polarisation-sum form, must equal resV
    res2 = sp.simplify(bil(q, qp)**2 - bil(q, q) * bil(qp, qp) / (Dd - 1))
    x = sp.simplify(1 + 2 * t / (s - 4 * m2))
    return dict(s=sp.nsimplify(s), t=sp.N(t), u=sp.N(u), x=sp.N(x), resV=sp.N(resV), resV_proj=sp.N(resV_proj),
                predV=sp.N((s - 4 * m2) * x), res2=sp.N(res2), pred2=sp.N((s - 4 * m2)**2 * (x**2 - sp.Rational(1, Dd - 1))))
for (mm, MM) in [(sp.Rational(1, 2), sp.Rational(3, 2)), (sp.Rational(1, 2), sp.Rational(7, 2))]:
    r = run(mm, MM, th)
    tag = 'BELOW' if MM < 4 * mm else 'ABOVE'
    print(f"m^2={mm}, M^2={MM} ({tag} threshold 4m^2={4*mm}):")
    for kk, v in r.items():
        print(f"    {kk:10s} = {v}")
    print(f"    => spin-1 residue = e^2 (M^2-4m^2) x : coefficient of x has sign {sp.sign(MM-4*mm)} for a HEALTHY real-e Proca field")
print()
print("Veneziano-type level-1 residue (al_t+1) with t=kappa(x-1): = (1-a)/2 + kappa x, kappa=(s_1-4m^2)/2;")
print("its spin-1 part kappa*x = (1/2)(M^2-4m^2) x is EXACTLY the Proca residue with e^2=1/2.  Hence c_{1,1}<0 for")
print("m^2>1/3 (2201.12331 eq.(13)) is the healthy sign of a real-coupled vector below threshold, not a ghost.")
