"""
R2: all-orders test of the CHR planar lambda-family against the master equation of
Huang-Wan-Wang-Zhou 2610.02124 eq. (45) (derived from N=4 SUSY + 6pt factorization + scalar parity
+ a discrete R-symmetry):
  (x1x3 - x2x4) B(x1,x2)B(x3,x4) + x2x4 B(x1,x2+x3)B(x3,x1+x4) = x1x3 B(x2,x1+x4)B(x2+x3,x4),
  x4 = -(x1+x2+x3),  A(s,u) = -B(s,u)/(su),  B(s,0)=1.
For the CHR family A_lam(s,u) = -rho/(su) + Sat_lam(s,u) we set B = -su A/rho = 1 - su Sat/rho.
We scan rho to minimise the residual: if no rho gives ~0 residual, the family member is excluded
(given the parity assumption).  lam=1, rho=1 must give ~0 (Veneziano).
"""
import mpmath as mp
mp.mp.dps = 25

def Sat(lam, s, u):
    lam = mp.mpf(lam)
    f = lambda k: lam*(3*lam-1)*mp.gamma(k+lam-lam*s)*mp.gamma(k+lam-lam*u)/((3*lam-1+2*k)*mp.factorial(k)*mp.gamma(k+2*lam-lam*s-lam*u))
    return mp.nsum(f, [0, mp.inf])

def B(lam, rho, s, u):
    if s == 0 or u == 0:
        return mp.mpf(1)
    return 1 - s*u*Sat(lam, s, u)/rho

def residual(lam, rho, x1, x2, x3):
    x4 = -(x1+x2+x3)
    lhs = (x1*x3 - x2*x4)*B(lam,rho,x1,x2)*B(lam,rho,x3,x4) + x2*x4*B(lam,rho,x1,x2+x3)*B(lam,rho,x3,x1+x4)
    rhs = x1*x3*B(lam,rho,x2,x1+x4)*B(lam,rho,x2+x3,x4)
    return lhs - rhs

pts = [(mp.mpf('0.21'), mp.mpf('-0.13'), mp.mpf('0.17')),
       (mp.mpf('0.31'), mp.mpf('0.07'), mp.mpf('-0.23')),
       (mp.mpf('-0.11'), mp.mpf('0.26'), mp.mpf('0.09'))]

if __name__ == "__main__":
    # sanity: lam=1, rho=1 -> B = Gamma(1-s)Gamma(1-u)/Gamma(1-s-u)
    s,u = mp.mpf('0.2'), mp.mpf('-0.35')
    print("lam=1 B vs Veneziano B:", mp.nstr(B(1,1,s,u),15), mp.nstr(mp.gamma(1-s)*mp.gamma(1-u)/mp.gamma(1-s-u),15))
    for lam in ['1', '0.5', '2']:
        print(f"lambda = {lam}")
        for rho in ['0.3','0.6','0.86','1.0','1.15','1.5','3']:
            r = [residual(lam, mp.mpf(rho), *p) for p in pts]
            print(f"   rho={rho:>5}: residuals =", [mp.nstr(x, 4) for x in r])

# residual is exactly quadratic in y=1/rho with zero constant term: r(y) = c1 y + c2 y^2.
# The nontrivial root rho* = -c2/c1 must be the same at every kinematic point for a solution to exist.
def rho_star(lam, p):
    r1 = residual(lam, mp.mpf(1), *p)        # y=1
    r2 = residual(lam, mp.mpf(2), *p)        # y=1/2
    # c1 + c2 = r1 ; c1/2 + c2/4 = r2  -> c2 = 2(r1 - 2 r2) ... solve
    c2 = (r1 - 2*r2)*2
    c1 = r1 - c2
    return -c2/c1
if __name__ == "__main__":
    print("point-wise root rho* of the master-equation residual (must be point-independent):")
    for lam in ['1','0.5','0.9','2']:
        print("  lambda =", lam, [mp.nstr(rho_star(lam,p),8) for p in pts])
