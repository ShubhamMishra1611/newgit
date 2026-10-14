"""
R2 independent checks of S3 (one-loop sigma-model beta functions), conventions of Polchinski (3.7.14):
  beta^G_mn   = a' R_mn + 2a' D_m D_n Phi - (a'/4) H_mpq H_n^pq
  beta^Phi    = (D-26)/6 - (a'/2) D^2 Phi + a' (D Phi)^2 - (a'/24) H^2
  S_spacetime ~ int sqrt(G) e^{-2Phi} [ -2(D-26)/(3a') + R + 4 (D Phi)^2 - H^2/12 ]   (Polchinski 3.7.20)
"""
import sympy as sp
import numpy as np
from scipy import integrate

def ricci(g, coords):
    n = len(coords); ginv = g.inv()
    Gam = [[[sum(ginv[a, d]*(sp.diff(g[d, b], coords[c]) + sp.diff(g[d, c], coords[b]) - sp.diff(g[b, c], coords[d]))
                 for d in range(n))/2 for c in range(n)] for b in range(n)] for a in range(n)]
    Ric = sp.zeros(n, n)
    for b in range(n):
        for c in range(n):
            Ric[b, c] = sp.simplify(sum(sp.diff(Gam[a][b][c], coords[a]) - sp.diff(Gam[a][b][a], coords[c])
                        + sum(Gam[a][a][d]*Gam[d][b][c] - Gam[a][c][d]*Gam[d][b][a] for d in range(n)) for a in range(n)))
    return Ric, Gam

# (1a) RNC on S^2: exact normal-coordinate metric d rho^2 + r^2 sin^2(rho/r) dphi^2 ; angular coefficient vs
#      prediction from g_ij = delta_ij - (1/3) R_ikjl x^k x^l with R_ikjl = (delta_ij delta_kl - delta_il delta_kj)/r^2
rho, r = sp.symbols('rho r', positive=True)
exact = sp.series(r**2*sp.sin(rho/r)**2, rho, 0, 6).removeO()
pred = rho**2 - rho**4/(3*r**2)          # angular part: rho^2 (1 - rho^2/(3 r^2))
print("(1a) S^2 RNC: exact - predicted angular metric through O(rho^4):", sp.simplify(exact - pred))

# (1b) coincident 2D lattice propagator  G(m^2) = int_{-pi}^{pi} d^2k/(2pi)^2 1/(4 sin^2(k1/2)+4 sin^2(k2/2)+m^2)
def Glat(m2):
    # integrate k2 analytically: int dk2/(2pi) 1/(A + 4 sin^2(k2/2)) = 1/sqrt(A(A+4))
    f = lambda k1: 1/np.sqrt((4*np.sin(k1/2)**2 + m2)*(4*np.sin(k1/2)**2 + m2 + 4))
    val, _ = integrate.quad(f, -np.pi, np.pi, points=[0], limit=400, epsabs=1e-13)
    return val/(2*np.pi)
m2s = np.array([1e-3, 1e-4, 1e-5, 1e-6])
Gs = np.array([Glat(m) for m in m2s])
slope = np.polyfit(np.log(m2s), Gs, 1)[0]
print(f"(1b) dG/dlog m^2 = {slope:.6f}  vs -1/(4 pi) = {-1/(4*np.pi):.6f};  G + log(m^2)/(4pi) = {Gs[-1]+np.log(m2s[-1])/(4*np.pi):.6f} vs log(32)/(4pi) = {np.log(32)/(4*np.pi):.6f}")

# (1c) O(N) model: sphere S^{N-1} radius R has R_mn = (N-2)/R^2 G_mn ; mu dG/dmu = a' R_mn => mu d(R^2)/dmu = a'(N-2)
#      g = 2 pi a'/R^2 => beta(g) = -(N-2) g^2/(2 pi)
for N in [3, 4, 5]:
    th = sp.symbols(f'th1:{N}')
    R = sp.symbols('R', positive=True)
    # round metric on S^{N-1}
    g = sp.zeros(N-1, N-1); fac = R**2
    for i in range(N-1):
        g[i, i] = fac
        fac = fac*sp.sin(th[i])**2
    Ric, _ = ricci(g, th)
    k = sp.simplify(Ric[0, 0]/g[0, 0])
    a, mu = sp.symbols("alpha mu", positive=True)
    gcoup = 2*sp.pi*a/R**2
    # mu dR^2/dmu = a' * k * R^2
    beta = sp.simplify(sp.diff(gcoup, R)*(a*k*R**2)/(2*R))
    print(f"(1c) O({N}): R_mn/G_mn = {k};  beta(g) = {sp.simplify(beta.subs(a, gcoup*R**2/(2*sp.pi)))}   (Polyakov -(N-2) g^2/(2 pi) = {sp.simplify(-(N-2)*gcoup**2/(2*sp.pi)).subs(a, gcoup*R**2/(2*sp.pi))})")

# (2) cigar: K(dr^2 + tanh^2 r dth^2), Phi = -log cosh r
rr, tt, K = sp.symbols('r theta K', positive=True)
g = sp.diag(K, K*sp.tanh(rr)**2); co = (rr, tt)
Ric, Gam = ricci(g, co)
Phi = -sp.log(sp.cosh(rr))
def hess(Phi):
    return sp.Matrix(2, 2, lambda m, n: sp.diff(Phi, co[m], co[n]) - sum(Gam[a][m][n]*sp.diff(Phi, co[a]) for a in range(2)))
H = hess(Phi)
ginv = g.inv()
lap = sp.simplify(sum(ginv[m, n]*H[m, n] for m in range(2) for n in range(2)))
grad2 = sp.simplify(sum(ginv[m, n]*sp.diff(Phi, co[m])*sp.diff(Phi, co[n]) for m in range(2) for n in range(2)))
print("(2) cigar: R_mn + 2 D_m D_n Phi =", sp.simplify(Ric + 2*H), "; -lap/2 + (DPhi)^2 =", sp.simplify(-lap/2 + grad2),
      "; D=2 & beta^Phi=0 => K =", sp.solve(sp.Rational(2-26, 6) + sp.Symbol('alpha')*sp.simplify(-lap/2 + grad2), K))

# (3) S^3 with H = (h/R^3) vol: R_mn - (1/4) H_mpq H_n^pq = (2/R^2 - h^2/(2R^6)) G_mn  (eps_mpq eps_n^pq = 2 g_mn)
R, h = sp.symbols('R h', positive=True)
print("(3) S^3: h solving 2/R^2 - h^2/(2 R^6) = 0:", sp.solve(2/R**2 - h**2/(2*R**6), h), "; H^2 = 6 h^2/R^6 at h=2R^2:", sp.simplify(6*(2*R**2)**2/R**6))

# (4) dilaton EOM vs beta combination in D dims with H=0, check for a 2D conformally-flat metric e^{2w} delta, general Phi(u,v)
u, v, a = sp.symbols('u v alpha')
w = sp.Function('w')(u, v); P = sp.Function('P')(u, v); Dd = sp.Symbol('D')
g = sp.diag(sp.exp(2*w), sp.exp(2*w)); co = (u, v)
Ric, Gam = ricci(g, co); ginv = g.inv()
Rs = sp.simplify(sum(ginv[m, n]*Ric[m, n] for m in range(2) for n in range(2)))
Hs = sp.Matrix(2, 2, lambda m, n: sp.diff(P, co[m], co[n]) - sum(Gam[q][m][n]*sp.diff(P, co[q]) for q in range(2)))
lapP = sp.simplify(sum(ginv[m, n]*Hs[m, n] for m in range(2) for n in range(2)))
dP2 = sp.simplify(sum(ginv[m, n]*sp.diff(P, co[m])*sp.diff(P, co[n]) for m in range(2) for n in range(2)))
sqrtg = sp.exp(2*w)
L = sqrtg*sp.exp(-2*P)*(-2*(Dd-26)/(3*a) + Rs + 4*dP2)
EL = sp.simplify(sp.diff(L, P) - sum(sp.diff(sp.diff(L, sp.diff(P, c)), c) for c in co))
betaG_tr = a*Rs + 2*a*lapP
betaPhi = (Dd-26)/6 - a/2*lapP + a*dP2
ratio = sp.simplify(EL/(sqrtg*sp.exp(-2*P)*(betaG_tr - 4*betaPhi)/a))
print("(4) dilaton Euler-Lagrange / [sqrt(G) e^{-2Phi}(G^mn beta^G_mn - 4 beta^Phi)/a'] =", ratio)
