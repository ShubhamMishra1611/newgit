"""Independent symbolic check (monomial -> Gegenbauer formula, NOT the recurrence) of low-level rho_j
for the open c=0 family on trajectory (m^2=-a) and the c=1 family; D symbolic.
x^k = k!/2^k sum_m (nu+k-2m)/(m! (nu)_{k-m+1}) C_{k-2m}(x)."""
import sympy as sp
a, D, x = sp.symbols('a D x')
nu = (D-3)/2
def gegen_coeffs(poly_x):
    P = sp.Poly(sp.expand(poly_x), x)
    deg = P.degree()
    c = [0]*(deg+1)
    for (k,), rk in P.terms():
        for m in range(k//2+1):
            j = k-2*m
            c[j] += rk*sp.factorial(k)/2**k*(nu+j)/(sp.factorial(m)*sp.rf(nu, k-m+1))
    return [sp.factor(sp.simplify(cj)) for cj in c]
def residue(n, cc, m2):
    kap = (n - a - 4*m2)/2
    t = kap*(x-1)
    R = sp.prod([a + i + t for i in range(1, n-cc+1)])/sp.factorial(n)
    return sp.expand(R), kap
out = []
for cc in [0,1]:
    for n in range(1, 6):
        R, kap = residue(n, cc, -a)
        cs = gegen_coeffs(R)
        rhos = [sp.factor(sp.simplify(cs[j]/kap**j)) for j in range(len(cs))]
        print(f'c={cc} level {n}: kappa={sp.factor(kap)}')
        for j, r in enumerate(rhos):
            print(f'   rho_{j} =', r)
