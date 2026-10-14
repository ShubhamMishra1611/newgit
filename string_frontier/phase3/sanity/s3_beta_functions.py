"""
Sanity check S3: leading-order (one-loop) sigma-model beta functions.

Target formulas (ESTABLISHED; Callan-Friedan-Martinec-Perry, Nucl. Phys. B262
(1985) 593 [pre-arXiv]; Polchinski Vol. 1, Sec. 3.7):
  beta^G_{mn}   = a' R_{mn} + 2a' D_m D_n Phi - (a'/4) H_{mpq} H_n^{pq}
  beta^B_{mn}   = -(a'/2) D^p H_{pmn} + a' D^p Phi H_{pmn}
  beta^Phi      = (D-26)/6 - (a'/2) D^2 Phi + a' (D Phi)^2 - (a'/24) H^2

What is computed here:
  (1) the background-field derivation of beta^G = a' R_{mn}:
      (1a) the Riemann-normal-coordinate coefficient g_ij = delta_ij - 1/3 R_ikjl x^k x^l
           checked symbolically on S^2 and H^2 (built from the exponential map);
      (1b) the coincident propagator <xi xi> = a' log(Lambda/mu), checked on a
           2D lattice: (1/V) sum_k 1/(khat^2+m^2) = -(1/4pi) log m^2 + const;
      (1c) combining: g_mn(mu) = g_mn - a' R_mn log(Lambda/mu)  =>  beta^G = a' R_mn,
           and cross-check against Polyakov's O(N) model result
           beta(g) = -(N-2) g^2/(2 pi)  (Polyakov, Phys. Lett. 59B (1975) 79).
  (2) exact backgrounds solve beta = 0 at this order (symbolic):
      flat linear dilaton (D = 26 - 6 a' V^2), Witten's 2D cigar, the S^3 WZW
      background with H flux.
  (3) the beta functions are the Euler-Lagrange equations of
      S = int sqrt(G) e^{-2Phi} [ -2(D-26)/(3a') + R + 4(DPhi)^2 - H^2/12 ]:
      dS/dPhi  is proportional to  (G^{mn} beta^G_mn - 4 beta^Phi),
      checked with sympy's Euler-Lagrange operator on a general 2D background.
"""
import itertools
import sympy as sp
import numpy as np

out = []
def log(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    out.append(s)

def christoffel(g, X):
    n = len(X)
    gi = g.inv()
    return [[[sp.simplify(sum(gi[a, d] * (sp.diff(g[d, b], X[c]) + sp.diff(g[d, c], X[b]) - sp.diff(g[b, c], X[d]))
                              for d in range(n)) / 2) for c in range(n)] for b in range(n)] for a in range(n)]

def riemann(g, X):
    """R^a_{bcd} = d_c Gam^a_{db} - d_d Gam^a_{cb} + Gam^a_{ce} Gam^e_{db} - Gam^a_{de} Gam^e_{cb}"""
    n = len(X)
    Gm = christoffel(g, X)
    R = [[[[sp.simplify(sp.diff(Gm[a][d][b], X[c]) - sp.diff(Gm[a][c][b], X[d])
                        + sum(Gm[a][c][e] * Gm[e][d][b] - Gm[a][d][e] * Gm[e][c][b] for e in range(n)))
            for d in range(n)] for c in range(n)] for b in range(n)] for a in range(n)]
    return R, Gm

def ricci(g, X):
    n = len(X)
    R, Gm = riemann(g, X)
    Ric = sp.Matrix(n, n, lambda b, d: sp.simplify(sum(R[a][b][a][d] for a in range(n))))
    return Ric, R, Gm

def hessian_scalar(f, g, X, Gm):
    n = len(X)
    return sp.Matrix(n, n, lambda a, b: sp.simplify(sp.diff(f, X[a], X[b]) - sum(Gm[c][a][b] * sp.diff(f, X[c]) for c in range(n))))

def laplacian(f, g, X):
    n = len(X)
    sg = sp.sqrt(sp.simplify(g.det()))
    gi = g.inv()
    return sp.simplify(sum(sp.diff(sg * gi[a, b] * sp.diff(f, X[b]), X[a]) for a in range(n) for b in range(n)) / sg)

# ---------------------------------------------------------------------------
log("=== (1a) Riemann normal coordinates: g_ij = delta_ij - (1/3) R_ikjl x^k x^l + O(x^3) ===")
rho, phi, r0 = sp.symbols('rho phi r0', positive=True)
x, y = sp.symbols('x y', real=True)
for name, fsin, K in (("S^2 radius r0", lambda u: r0 * sp.sin(u / r0), 1 / r0**2),
                      ("H^2 radius r0", lambda u: r0 * sp.sinh(u / r0), -1 / r0**2)):
    # geodesic polar coordinates = exponential map: ds^2 = drho^2 + f(rho)^2 dphi^2
    # Cartesian RNC x = rho cos phi, y = rho sin phi:
    # g_ij = rhohat_i rhohat_j + (f^2/rho^2)(delta_ij - rhohat_i rhohat_j)
    rr = sp.sqrt(x**2 + y**2)
    f2_over_rho2 = sp.series((fsin(rho) / rho)**2, rho, 0, 4).removeO().subs(rho, rr)
    X = [x, y]
    gRNC = sp.Matrix(2, 2, lambda i, j: sp.expand((X[i] * X[j] / rr**2) * (1 - f2_over_rho2)
                                                    + (1 if i == j else 0) * f2_over_rho2))
    gRNC = gRNC.applyfunc(lambda e: sp.simplify(e))
    # constant-curvature Riemann tensor (all indices down, flat at origin): R_ikjl = K (d_ij d_kl - d_il d_kj)
    d = lambda a, b: 1 if a == b else 0
    pred = sp.Matrix(2, 2, lambda i, j: d(i, j) - sp.Rational(1, 3) * sum(
        K * (d(i, j) * d(k, l) - d(i, l) * d(k, j)) * X[k] * X[l] for k in range(2) for l in range(2)))
    diff = (gRNC - pred).applyfunc(lambda e: sp.simplify(sp.expand(e)))
    # keep only terms up to second order
    diff2 = diff.applyfunc(lambda e: sp.simplify(sp.series(e.subs({x: sp.Symbol('eps') * x, y: sp.Symbol('eps') * y}),
                                                           sp.Symbol('eps'), 0, 3).removeO()))
    log(f"  {name}: g_RNC - [delta - R x x/3] through O(x^2) = {diff2.tolist()}")

# ---------------------------------------------------------------------------
log("\n=== (1b) coincident propagator on a 2D lattice: coefficient of -log m^2 ===")
def lattice_tadpole(L, m2):
    k = 2 * np.pi * np.arange(L) / L
    kh2 = (2 * np.sin(k / 2))**2
    K2 = kh2[:, None] + kh2[None, :]
    return np.mean(1.0 / (K2 + m2))
L = 4096
m2s = np.array([1e-2, 4e-3, 1e-3, 4e-4])
vals = np.array([lattice_tadpole(L, m2) for m2 in m2s])
slope = np.polyfit(np.log(m2s), vals, 1)[0]
log(f"  L={L}: tadpole vs m^2 = {list(zip(m2s.tolist(), np.round(vals, 6).tolist()))}")
log(f"  fitted d(tadpole)/d(log m^2) = {slope:.6f};   -1/(4 pi) = {-1/(4*np.pi):.6f}")
log(f"  constant: tadpole + log(m^2)/(4pi) = {np.round(vals + np.log(m2s)/(4*np.pi), 5).tolist()} vs log(32)/(4pi) = {np.log(32)/(4*np.pi):.5f}")
log("  => int d^2k/(2pi)^2 1/k^2 = (1/2pi) log(Lambda/mu); with S = (1/4pi a') int (d xi)^2 the propagator is")
log("     2 pi a'/k^2, so <xi^a xi^b> = delta^ab a' log(Lambda/mu).")

# ---------------------------------------------------------------------------
log("\n=== (1c) beta^G = a' R_mn, and Polyakov's O(N) beta function ===")
log("  second-order background-field action: (1/4pi a') int R_{m a b n} xi^a xi^b dX^m dX^n  (RNC, from (1a))")
log("  <...> -> (1/4pi a') a' log(Lambda/mu) R_{m a a n} dX dX = -(1/4pi a') a' log(Lambda/mu) R_mn dX dX")
log("  => G_mn(mu) = G_mn - a' R_mn log(Lambda/mu),  mu dG/dmu = a' R_mn.")
ap = sp.symbols("alpha'", positive=True)
for N in (3, 4, 5):
    n = N - 1
    th = sp.symbols(f'th0:{n}')
    r = sp.symbols('r', positive=True)
    # round S^{N-1} metric of radius r
    diag = []
    pref = r**2
    for i in range(n):
        diag.append(pref)
        pref = pref * sp.sin(th[i])**2
    g = sp.diag(*diag)
    Ric, _, _ = ricci(g, list(th))
    ratio = sp.simplify(Ric[0, 0] / g[0, 0])
    gcoup = sp.symbols('g', positive=True)
    # d r^2 / d log mu = a' * ratio * r^2 ; g = 2 pi a'/r^2
    dr2 = ap * ratio * r**2
    dg = sp.simplify(-2 * sp.pi * ap / r**4 * dr2)
    beta_g = sp.simplify(dg.subs(r, sp.sqrt(2 * sp.pi * ap / gcoup)))
    log(f"  O({N}) model, target S^{N-1}: R_mn = {sp.simplify(ratio*r**2)}/r^2 G_mn ;  with g = 2 pi a'/r^2:  "
        f"beta(g) = {beta_g}   (Polyakov: {sp.simplify(-(N-2)*gcoup**2/(2*sp.pi))})")

# ---------------------------------------------------------------------------
log("\n=== (2) exact backgrounds solve the one-loop beta equations ===")
Dd = sp.symbols('D')
# (2a) linear dilaton
V = sp.symbols('V', positive=True)
log(f"  linear dilaton (flat G, Phi = V.X): beta^G = 2a' d d Phi = 0;  beta^Phi = (D-26)/6 + a' V^2 = 0 "
    f"=> D = {sp.solve(sp.Eq((Dd-26)/6 + ap*V**2, 0), Dd)[0]}")
# (2b) Witten cigar
rr_ = sp.symbols('r', positive=True); th_ = sp.symbols('theta', real=True)
k = sp.symbols('k', positive=True)
g = sp.diag(k, k * sp.tanh(rr_)**2)
X = [rr_, th_]
Phi = -sp.log(sp.cosh(rr_))
Ric, _, Gm = ricci(g, X)
H = hessian_scalar(Phi, g, X, Gm)
betaG = (Ric + 2 * H).applyfunc(sp.simplify)
gi = g.inv()
dPhi2 = sp.simplify(sum(gi[a, b] * sp.diff(Phi, X[a]) * sp.diff(Phi, X[b]) for a in range(2) for b in range(2)))
lapPhi = laplacian(Phi, g, X)
betaPhi_noconst = sp.simplify(-sp.Rational(1, 2) * lapPhi + dPhi2)
log(f"  2D cigar ds^2 = k(dr^2 + tanh^2 r dth^2), Phi = -log cosh r:")
log(f"    R_mn + 2 D_m D_n Phi = {betaG.tolist()}")
log(f"    -(1/2) D^2 Phi + (D Phi)^2 = {betaPhi_noconst}  (constant)  => (D-26)/6 + a'/k = 0, D=2 => k = {sp.solve(sp.Eq(sp.Rational(2-26,6) + ap/k, 0), k)}")
# (2c) S^3 with H flux
ps, th2, ph2 = sp.symbols('psi theta2 phi2', real=True)
R3 = sp.symbols('R', positive=True)
hflux = sp.symbols('h', positive=True)
g3 = sp.diag(R3**2, R3**2 * sp.sin(ps)**2, R3**2 * sp.sin(ps)**2 * sp.sin(th2)**2)
X3 = [ps, th2, ph2]
Ric3, _, _ = ricci(g3, X3)
sqrtg = sp.sqrt(g3.det())
# H = h * (volume form)/R^3 ; H_{mnp} = h/R^3 * sqrt(g) eps_{mnp}
def Hdn(a, b, c):
    return hflux / R3**3 * sqrtg * sp.LeviCivita(a, b, c)
g3i = g3.inv()
HH = sp.Matrix(3, 3, lambda m, n: sp.simplify(sum(Hdn(m, p, q) * Hdn(n, pp, qq) * g3i[p, pp] * g3i[q, qq]
                                                  for p in range(3) for q in range(3) for pp in range(3) for qq in range(3))))
betaG3 = (Ric3 - HH / 4).applyfunc(sp.simplify)
hsol = sp.solve(sp.Eq(betaG3[0, 0], 0), hflux)
log(f"  S^3 radius R with H = (h/R^3) vol:  R_mn - H_mpq H_n^pq/4 = {[sp.factor(betaG3[i,i]) for i in range(3)]} (diag) => h = {hsol}")
log("    (H proportional to the covariantly-constant volume form => D^p H_pmn = 0, so beta^B = 0 with const Phi)")
H2 = sp.simplify(sum(HH[m, n] * g3i[m, n] for m in range(3) for n in range(3)))
log(f"    H^2 = {sp.simplify(H2.subs(hflux, hsol[0]))};  beta^Phi = (D-26)/6 - a' H^2/24 = (D-26)/6 - a'/R^2  (const, fixes the CFT c)")

# ---------------------------------------------------------------------------
log("\n=== (3) beta functions as Euler-Lagrange equations of the spacetime action (2D, H=0) ===")
u, v = sp.symbols('u v', real=True)
w = sp.Function('w')(u, v)        # conformal factor
P = sp.Function('P')(u, v)        # dilaton
g2 = sp.diag(sp.exp(2 * w), sp.exp(2 * w))
X2 = [u, v]
Ric2, _, Gm2 = ricci(g2, X2)
Rs = sp.simplify(sum(g2.inv()[a, b] * Ric2[a, b] for a in range(2) for b in range(2)))
dP2 = sp.simplify(sum(g2.inv()[a, b] * sp.diff(P, X2[a]) * sp.diff(P, X2[b]) for a in range(2) for b in range(2)))
lapP = laplacian(P, g2, X2)
cst = sp.symbols('Lambda0')   # stands for 2(D-26)/(3a')
Lag = sp.exp(2 * w) * sp.exp(-2 * P) * (-cst + Rs + 4 * dP2)
EL_P = sp.euler_equations(Lag, [P], X2)[0].lhs
# beta combination (in units of a'): trace(beta^G)/a' - 4 beta^Phi/a'  with (D-26)/(6a') = cst/4
HP = hessian_scalar(P, g2, X2, Gm2)
trbetaG = sp.simplify(Rs + 2 * lapP)
betaPhi = cst / 4 - sp.Rational(1, 2) * lapP + dP2
combo = sp.simplify(trbetaG - 4 * betaPhi)
ratio = sp.simplify(EL_P / (sp.exp(2 * w) * sp.exp(-2 * P) * combo))
log(f"  dL/dPhi (Euler-Lagrange) / [sqrt(G) e^(-2Phi) (G^mn beta^G_mn - 4 beta^Phi)/a'] = {ratio}")
log("  => the dilaton equation of motion IS the stated combination of one-loop beta functions.")
# w-variation (trace of metric equation): in 2D, R sqrt(g) is a total derivative
EL_w = sp.euler_equations(Lag, [w], X2)[0].lhs
combo_w = sp.simplify(EL_w / (sp.exp(2 * w) * sp.exp(-2 * P)))
log(f"  dL/dw / (sqrt(G)e^(-2Phi)) = {sp.simplify(combo_w)}")
log(f"  compare  -8*beta^Phi/a' = {sp.simplify(-8*betaPhi)}  -> difference = {sp.simplify(combo_w + 8*betaPhi)}")

with open(__file__.replace('.py', '_output.txt'), 'w') as f:
    f.write("\n".join(out) + "\n")
