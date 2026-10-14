"""T3 (B1): late-time exponential slopes from 10d -> 4d (dilaton, volume) scalings of type II sources.

Run:  python3 -I T3_B1_swampland_slopes.py > T3_B1_swampland_slopes.out

DERIVATION (string frame -> 4d Einstein frame), conventions of Hertzberg-Kachru-Taylor-Tegmark
(HKTT, 0711.2512; read this session):
  rho = Vol^{1/3}, tau = e^{-phi} sqrt(Vol) = e^{-phi} rho^{3/2}        (HKTT eq. 2.3)
  g^E_{mu nu} = tau^2 g_{mu nu}  (units M_P = 1)                          (HKTT eq. 2.4)
  canonical fields rhohat = sqrt(3/2) ln rho, tauhat = sqrt(2) ln tau      (HKTT eq. 2.5)
  We re-derive (2.5) below with sympy: reducing e^{-2phi} sqrt(G)(R + 4(d phi)^2) on a 6d space of
  metric e^{2B} g6 gives, in 4d Einstein frame, (1/2)R - (d ln tau)^2 - (3/4)(d ln rho)^2  [rho=e^{2B}].
A 10d string-frame term with string-coupling weight e^{-chi phi} (chi = Euler number of the
worldsheet: 2 sphere, 1 disk/RP2, 0 torus or RR kinetic term) and internal scaling rho^w
(w = (volume factor) - (number of inverse internal metrics)) gives
  V_E  propto  tau^{-4} e^{-chi phi} rho^w = rho^{w - 3 chi/2} tau^{chi - 4}.
  term                 chi   w          V_E                         source of the scaling
  R_6 (curvature)       2    3-1        rho^-1 tau^-2               HKTT (4.2) [geometric flux f]
  H_3                   2    3-3        rho^-3 tau^-2               HKTT (3.3)
  F_p (internal)        0    3-p        rho^{3-p} tau^-4            HKTT (3.3)
  Dp/Op on (p-3)-cycle  1    (p-3)/2    +-rho^{(p-6)/2} tau^-3      HKTT (4.1)  (O-planes negative)
  NS5 on 2-cycle        2    1          rho^-2 tau^-2               HKTT (4.4)
  Q / R non-geom. flux  2    -          +-rho^{1} tau^-2, +-rho^{3} tau^-2   HKTT (4.2)
  Casimir (1-loop)      0    3-5        +-rho^-2 tau^-4             our derivation: 10d energy density ~ L^-10
  10d 1-loop Lambda     0    3          rho^{3} tau^-4              our derivation (non-SUSY strings)
  tree-level dilaton    2    3          tau^-2                      our derivation (sphere-level constant,
       tadpole (formal)                                              e.g. non-critical strings) -- formal only
Writing V = A e^{-k.phihat} with phihat = (rhohat, tauhat):  k = (-a sqrt(2/3), -b/sqrt(2)) for V ~ rho^a tau^b.
Note the tauhat component of every k is (4-chi)/sqrt2 >= sqrt2: this is Shiu-Tonioni-Tran's
eq. (II.3) (2303.03418) gamma_dilaton^2 >= 4/(d-2) specialised to d = 4.

LATE-TIME SLOPE.  For V = sum_i A_i e^{-k_i.phi} with all A_i > 0 the late-time scaling
attractor has |grad V|/V -> distance from the origin to conv{k_i} (Calderon-Infante-Ruiz-
Valenzuela 2209.11821 'convex hull dS conjecture'; Shiu-Tonioni-Tran 2306.07327 convex-hull
criterion), and epsilon = gamma^2/2 in d=4 (acceleration iff gamma < sqrt2).
With terms of both signs we find all scaling fixed points of the autonomous system of
2303.03418 eqs. (A.1)-(A.2) (d=4: c_i = (1/2) sqrt(2/3) k_i, epsilon = 3 x^2) and keep the
linearly stable ones with V > 0.
"""
import itertools, math, json, os
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp

S23 = math.sqrt(2 / 3)


def hdr(s):
    print("\n" + "=" * 90 + "\n" + s + "\n" + "=" * 90, flush=True)


# ------------------------------------------------------------------------------------------
hdr("0. sympy minisuperspace check of the 4d Einstein-frame kinetic terms (HKTT eq. 2.5)")
# 10d string frame: ds^2 = -N^2 dt^2 + a^2 dx_3^2 + b^2 ds^2(T^6, unit volume), dilaton phi(t).
# L10 = sqrt(-G) e^{-2phi} (R10 + 4 (dphi)^2).  Remove second derivatives by parts, then rewrite in
# 4d Einstein-frame variables a_E = tau a, N_E = tau N, tau = e^{-phi} b^3, rho = b^2 and compare with
# 2 x [ -3 a_E adot_E^2/N_E + (a_E^3/N_E)( (dln tau)^2 + (3/4)(dln rho)^2 ) ],
# i.e. L4 = R/2 - (d ln tau)^2 - (3/4)(d ln rho)^2  <=>  tauhat = sqrt2 ln tau, rhohat = sqrt(3/2) ln rho.
tt = sp.symbols('t')
N, a, b, ph = [sp.Function(n)(tt) for n in ('N', 'a', 'b', 'phi')]
dim = 10
hs = [a] * 3 + [b] * 6
g = sp.diag(*([-N ** 2] + [h ** 2 for h in hs]))
coords = [tt] + list(sp.symbols('x1:10'))
ginv = g.inv()
def Gam(l, m, n_):
    return sum(ginv[l, s_] * (sp.diff(g[s_, m], coords[n_]) + sp.diff(g[s_, n_], coords[m]) - sp.diff(g[m, n_], coords[s_]))
               for s_ in range(dim)) / 2
Gm = [[[sp.simplify(Gam(l, m, n_)) for n_ in range(dim)] for m in range(dim)] for l in range(dim)]
def Ric(m, n_):
    r = 0
    for l in range(dim):
        r += sp.diff(Gm[l][m][n_], coords[l]) - sp.diff(Gm[l][m][l], coords[n_])
        for s_ in range(dim):
            r += Gm[l][l][s_] * Gm[s_][m][n_] - Gm[l][n_][s_] * Gm[s_][m][l]
    return r
Rs = sp.simplify(sum(ginv[m, m] * Ric(m, m) for m in range(dim)))
sqrtG = N * a ** 3 * b ** 6
L10 = sqrtG * sp.exp(-2 * ph) * (Rs - 4 * sp.diff(ph, tt) ** 2 / N ** 2)
L10 = sp.expand(L10)
# integrate by parts the second-derivative terms
L1 = L10
for f in (a, b):
    f2 = sp.diff(f, tt, 2)
    A_ = sp.expand(L1).coeff(f2)
    L1 = sp.expand(L1 - A_ * f2 - sp.diff(A_, tt) * sp.diff(f, tt))
assert not L1.has(sp.Derivative(a, (tt, 2))) and not L1.has(sp.Derivative(b, (tt, 2)))
# Einstein-frame variables
aE, NE, lt, lr = [sp.Function(n)(tt) for n in ('aE', 'NE', 'lntau', 'lnrho')]
sub = {b: sp.exp(lr / 2), ph: -lt + sp.Rational(3, 2) * lr, a: aE * sp.exp(-lt), N: NE * sp.exp(-lt)}
L1E = sp.simplify(sp.expand(L1.subs(sub).doit()))
target = 2 * (-3 * aE * sp.diff(aE, tt) ** 2 / NE + aE ** 3 / NE * (sp.diff(lt, tt) ** 2 + sp.Rational(3, 4) * sp.diff(lr, tt) ** 2))
diff_ = sp.simplify(sp.expand(L1E - target))
print("  L10 (after parts, Einstein-frame variables) - 2*L4(target) =", diff_)
print("  -> canonical fields tauhat = sqrt(2) ln tau, rhohat = sqrt(3/2) ln rho confirmed" if diff_ == 0 else "  MISMATCH")


def kvec(a, b):
    """V ~ rho^a tau^b -> k with V ~ e^{-k.(rhohat,tauhat)}"""
    return np.array([-a * S23, -b / math.sqrt(2)])


def scaling_from_chi_w(chi, w):
    return (w - 1.5 * chi, chi - 4)


# Catalogue: name -> (a, b, sign, theory, comment)
cat = {}
def add(name, chi, w, sign, theory, note):
    a, b = scaling_from_chi_w(chi, w)
    cat[name] = dict(a=a, b=b, k=kvec(a, b), sign=sign, theory=theory, note=note, chi=chi, w=w)

add("R6(<0 curv)", 2, 2, "+", "II", "negative internal curvature (V>0); same scaling as geometric flux / KK monopoles")
add("R6(>0 curv)", 2, 2, "-", "II", "positive internal curvature (V<0)")
add("H3", 2, 0, "+", "II", "NSNS 3-form flux")
for p in range(0, 7):
    add(f"F{p}", 0, 3 - p, "+", "IIA" if p % 2 == 0 else "IIB", f"internal RR {p}-form flux")
for p in range(4, 10):
    th = "IIA" if p % 2 == 0 else "IIB"
    if p == 9:
        th = "IIB/typeI"
    add(f"D{p}", 1, (p - 3) / 2, "+", th, f"(anti-)D{p} on a {p-3}-cycle, positive tension")
    add(f"O{p}", 1, (p - 3) / 2, "-", th, f"O{p}-plane, negative tension")
add("D3", 1, 0, "+", "IIB", "(anti-)D3 at points")
add("O3", 1, 0, "-", "IIB", "O3-plane")
add("NS5", 2, 1, "+", "II", "(anti-)NS5 on 2-cycle (HKTT 4.4)")
cat["Qflux"] = dict(a=1, b=-2, k=kvec(1, -2), sign="+-", theory="II(nongeom)", note="non-geometric Q flux (HKTT 4.2)", chi=2, w=None)
cat["Rflux"] = dict(a=3, b=-2, k=kvec(3, -2), sign="+-", theory="II(nongeom)", note="non-geometric R flux (HKTT 4.2)", chi=2, w=None)
add("Casimir", 0, -2, "+-", "any", "one-loop Casimir energy, sign depends on spectrum/BCs")
add("Lambda10(1-loop)", 0, 3, "+", "nonSUSY", "10d one-loop vacuum energy of a non-SUSY string (e.g. O(16)xO(16))")
add("dil.tadpole(tree)", 2, 3, "+", "formal", "sphere-level 10d constant (non-critical strings) -- formal")

print("\n  term                 V_E ~ rho^a tau^b      k = (k_rhohat, k_tauhat)     |k|     sign  theory")
for n, d in cat.items():
    print(f"  {n:20s} a={d['a']:+5.1f} b={d['b']:+3d}     ({d['k'][0]:+.4f}, {d['k'][1]:+.4f})       {np.linalg.norm(d['k']):.4f}   {d['sign']:3s}  {d['theory']}")
print(f"\n  min over all terms of k_tauhat = {min(d['k'][1] for d in cat.values()):.4f}  (sqrt2 = {math.sqrt(2):.4f}):"
      " every term has dilaton coupling >= sqrt2 [STT 2303.03418 eq. II.3 at d=4]")


# ------------------------------------------------------------------------------------------
def hull_distance(K):
    """distance from origin to convex hull of rows of K (2D), plus the closest point."""
    K = np.asarray(K, float)
    best = (np.inf, None)
    for i in range(len(K)):
        d = np.linalg.norm(K[i])
        if d < best[0]:
            best = (d, K[i].copy())
    for i, j in itertools.combinations(range(len(K)), 2):
        a, b = K[i], K[j]
        ab = b - a
        L2 = ab @ ab
        if L2 == 0:
            continue
        t = -(a @ ab) / L2
        if 0 < t < 1:
            p = a + t * ab
            d = np.linalg.norm(p)
            if d < best[0]:
                best = (d, p)
    # origin inside hull?  (2D test via angles)
    if len(K) >= 3:
        ang = np.sort(np.arctan2(K[:, 1], K[:, 0]))
        gaps = np.diff(np.concatenate([ang, [ang[0] + 2 * math.pi]]))
        if gaps.max() < math.pi:
            return 0.0, np.zeros(2)
    return best


def control_of_direction(theta):
    """direction of field motion theta (unit) in (rhohat,tauhat): does rho grow, does 10d g_s shrink?
    ln g_s = phi = -ln tau + (3/2) ln rho = -tauhat/sqrt2 + sqrt(3/2) rhohat."""
    dlnrho = theta[0] / math.sqrt(1.5)
    dphi10 = -theta[1] / math.sqrt(2) + math.sqrt(1.5) * theta[0]
    return dlnrho, dphi10


hdr("1. All-positive combinations: late-time slope = distance from origin to the convex hull")
pos_terms = {n: d for n, d in cat.items() if "+" in d["sign"]}
groups = {
    "IIA geometric (R6,H3,F0..F6 even,D4,D6,D8,NS5)": [n for n in pos_terms if n in ("R6(<0 curv)", "H3", "F0", "F2", "F4", "F6", "D4", "D6", "D8", "NS5")],
    "IIB geometric (R6,H3,F1,F3,F5,D3,D5,D7,D9,NS5)": [n for n in pos_terms if n in ("R6(<0 curv)", "H3", "F1", "F3", "F5", "D3", "D5", "D7", "D9", "NS5")],
    "IIA geometric + Casimir(+)": [n for n in pos_terms if n in ("R6(<0 curv)", "H3", "F0", "F2", "F4", "F6", "D4", "D6", "D8", "NS5", "Casimir")],
    "IIB geometric + Casimir(+)": [n for n in pos_terms if n in ("R6(<0 curv)", "H3", "F1", "F3", "F5", "D3", "D5", "D7", "D9", "NS5", "Casimir")],
    "+ non-SUSY 10d terms (Lambda10 1-loop, D9-like)": [n for n in pos_terms if n not in ("Qflux", "Rflux", "dil.tadpole(tree)")],
    "+ non-geometric Q,R fluxes (positive)": [n for n in pos_terms if n != "dil.tadpole(tree)"],
    "+ formal tree-level dilaton tadpole": list(pos_terms),
}
res1 = {}
for g, names in groups.items():
    K = np.array([cat[n]["k"] for n in names])
    d, p = hull_distance(K)
    # which terms realise the minimum (on the closest edge/vertex)
    phat = p / np.linalg.norm(p)
    on = [n for n in names if abs((cat[n]["k"] - p) @ phat) < 1e-9]   # terms on the supporting line
    theta = p / np.linalg.norm(p)
    dl, dp = control_of_direction(theta)
    res1[g] = (d, p.tolist(), on)
    print(f"  {g}:\n      min slope gamma = {d:.4f}   closest point {np.round(p,4)}  realised by {on}"
          f"\n      epsilon = gamma^2/2 = {d*d/2:.4f};  w_inf = -1 + gamma^2/3 = {-1+d*d/3:+.4f};"
          f"  along the attractor: d ln rho/d s = {dl:+.3f}, d ln g_s(10d)/d s = {dp:+.3f}")

print("\n  single-term slopes (each alone, positive):")
for n in ("R6(<0 curv)", "NS5", "H3", "D9", "D8", "D7", "D6", "F0", "F1", "Casimir", "Lambda10(1-loop)", "Qflux"):
    print(f"     {n:18s} |k| = {np.linalg.norm(cat[n]['k']):.4f}")

print("\n  pairwise hull distances among the low-lying positive terms (all-positive pairs):")
low = ["R6(<0 curv)", "NS5", "H3", "D9", "D8", "D7", "D6", "D5", "D4", "F0", "F1", "F2", "Casimir", "Lambda10(1-loop)", "Qflux"]
pairs = []
for a, b in itertools.combinations(low, 2):
    d, p = hull_distance(np.array([cat[a]["k"], cat[b]["k"]]))
    pairs.append((d, a, b))
pairs.sort()
for d, a, b in pairs[:12]:
    print(f"     {a:18s} + {b:18s}: gamma = {d:.4f}")


# ------------------------------------------------------------------------------------------
hdr("2. Mixed-sign combinations: stable scaling fixed points of the 2303.03418 autonomous system")
# variables: x (2-vector), w_i = y_i^2 (sign = sign of A_i);  d/dN x = -x(1-x^2) + sum c_i w_i,
# d/dN w_i = 2 (x^2 - c_i.x) w_i, constraint x^2 + sum w = 1, epsilon = 3 x^2 (d=4).
def c_of(k):
    return 0.5 * math.sqrt(2 / 3) * np.asarray(k)


def rhs(xw, C):
    x = xw[:2]
    w = xw[2:]
    x2 = x @ x
    dx = -x * (1 - x2) + C.T @ w
    dw = 2 * (x2 - C @ x) * w
    return np.concatenate([dx, dw])


def reduced_jac(x, w, C, elim):
    """Jacobian on the constraint surface: eliminate w[elim] = 1 - x^2 - sum others."""
    m = len(w)
    keep = [i for i in range(m) if i != elim]
    def f(z):
        xx = z[:2]
        ww = np.zeros(m)
        ww[keep] = z[2:]
        ww[elim] = 1 - xx @ xx - ww[keep].sum()
        full = rhs(np.concatenate([xx, ww]), C)
        return np.concatenate([full[:2], full[2:][keep]])
    z0 = np.concatenate([x, w[keep]])
    n = len(z0)
    J = np.zeros((n, n))
    h = 1e-7
    for j in range(n):
        e = np.zeros(n); e[j] = h
        J[:, j] = (f(z0 + e) - f(z0 - e)) / (2 * h)
    return J


def fixed_points(names, signs):
    """scaling fixed points with participating subset S (|S|=1,2) of the given terms."""
    K = np.array([cat[n]["k"] for n in names])
    C = np.array([c_of(k) for k in K])
    m = len(names)
    out = []
    for r in (1, 2):
        for S in itertools.combinations(range(m), r):
            CS = C[list(S)]
            if r == 1:
                c = CS[0]
                if c @ c >= 1:
                    continue
                x = c.copy()
                wS = np.array([1 - x @ x])
            else:
                M = CS
                if abs(np.linalg.det(M)) < 1e-12:
                    continue
                u = np.linalg.solve(M, np.ones(2))
                x = u / (u @ u)
                if x @ x >= 1:
                    continue
                try:
                    wS = np.linalg.solve(M.T, x * (1 - x @ x))
                except np.linalg.LinAlgError:
                    continue
            w = np.zeros(m)
            w[list(S)] = wS
            ok = all((w[i] > 1e-12 and signs[i] == "+") or (w[i] < -1e-12 and signs[i] == "-") for i in S)
            if not ok:
                continue
            # stability: eliminate one participating w
            J = reduced_jac(x, w, C, S[0])
            ev = np.linalg.eigvals(J)
            stable = np.all(ev.real < -1e-9)
            out.append(dict(S=[names[i] for i in S], x=x, w=w, eps=3 * x @ x, stable=bool(stable),
                            ev=np.round(np.sort(ev.real), 4).tolist()))
    return out


neg_terms = [n for n, d in cat.items() if "-" in d["sign"]]
posable = [n for n, d in cat.items() if "+" in d["sign"]]
print(f"  negative-capable terms: {neg_terms}")
results = []
checked = 0
for P in posable:
    for N in neg_terms:
        if P == N:
            continue
        for extra in [None] + [e for e in posable if e not in (P,)]:
            names = [P, N] + ([extra] if extra and extra != N else [])
            signs = ["+", "-"] + (["+"] if extra and extra != N else [])
            fps = fixed_points(names, signs)
            checked += 1
            for fp in fps:
                if fp["stable"] and fp["eps"] < 1 + 1e-9 and any(cat[s]["sign"] == "-" or s == N for s in fp["S"]):
                    results.append((fp["eps"], names, fp))
print(f"  checked {checked} combinations (1 positive + 1 negative [+1 positive]) ")
if results:
    results.sort(key=lambda t: t[0])
    print("  STABLE mixed-sign scaling attractors with V>0 and epsilon <= 1 (accelerating if < 1):")
    seen = set()
    for eps, names, fp in results[:25]:
        key = (tuple(fp["S"]),)
        if key in seen:
            continue
        seen.add(key)
        theta = fp["x"] / np.linalg.norm(fp["x"])
        dl, dp = control_of_direction(theta)
        print(f"     eps={eps:.4f} gamma*={math.sqrt(2*eps):.4f}  participating {fp['S']}  w={np.round(fp['w'][:len(fp['S'])+1],3)}"
              f"  dlnrho={dl:+.2f} dln gs={dp:+.2f}  (in set {names})")
else:
    print("  NONE: no stable mixed-sign scaling solution with V>0 and epsilon<=1 among these combinations.")

# also report the lowest-epsilon stable mixed-sign attractor regardless of acceleration
allstable = []
for P in posable:
    for N in neg_terms:
        names = [P, N]
        for fp in fixed_points(names, ["+", "-"]):
            if fp["stable"] and len(fp["S"]) == 2:
                allstable.append((fp["eps"], fp["S"], fp))
allstable.sort(key=lambda t: t[0])
print("\n  lowest-epsilon STABLE two-term (one +, one -) scaling attractors with V>0:")
for eps, S, fp in allstable[:8]:
    theta = fp["x"] / np.linalg.norm(fp["x"])
    dl, dp = control_of_direction(theta)
    print(f"     eps={eps:.4f}  gamma*={math.sqrt(2*eps):.4f}  {S}  dlnrho={dl:+.2f}  dln gs={dp:+.2f}")
print(f"  number of stable 2-term mixed attractors found: {len(allstable)}")

# numerical confirmation by integrating the ODE for the minimal all-positive IIA/IIB sets
hdr("3. ODE check of the convex-hull prediction (random initial data, 60 e-folds)")
rng = np.random.default_rng(3)
for g in ("IIA geometric (R6,H3,F0..F6 even,D4,D6,D8,NS5)", "IIB geometric (R6,H3,F1,F3,F5,D3,D5,D7,D9,NS5)",
          "+ non-geometric Q,R fluxes (positive)"):
    names = groups[g]
    C = np.array([c_of(cat[n]["k"]) for n in names])
    finals = []
    for trial in range(5):
        x0 = rng.normal(size=2) * 0.3
        w0 = rng.random(len(names)) + 0.05
        w0 *= (1 - x0 @ x0) / w0.sum()
        sol = solve_ivp(lambda N, z: rhs(z, C), (0, 60), np.concatenate([x0, w0]), rtol=1e-10, atol=1e-12)
        xf = sol.y[:2, -1]
        finals.append(3 * xf @ xf)
    d = res1[g][0]
    print(f"  {g[:45]:45s}: late-time epsilon from ODE = {np.round(finals,4)}  vs hull prediction gamma^2/2 = {d*d/2:.4f}")

hdr("4. Comparison with bounds and with DESI DR2 exponential quintessence (2504.04226)")
lam, up, dn = 0.698, 0.173, 0.202   # abstract of 2504.04226 (Omega_K = 0)
lamK, upK, dnK = 0.722, 0.182, 0.208
bounds = {"TCC 2/sqrt((d-1)(d-2)) [1909.11063]": 2 / math.sqrt(6),
          "strong/asymptotic 2/sqrt(d-2) [2101.11617]": math.sqrt(2)}
for g in groups:
    bounds[f"hull min: {g}"] = res1[g][0]
for k_, b in bounds.items():
    print(f"  {k_:62s} gamma = {b:.4f}:  (gamma - 0.698)/0.173 = {(b-lam)/up:+.2f} sigma"
          f"   [Omega_K free: (gamma-0.722)/0.182 = {(b-lamK)/upK:+.2f}]")
print("  NB: one-sided Gaussian reading of an asymmetric posterior; the DESI-era fit constrains the slope of the"
      " potential traversed since z~few, not the asymptotic (infinite-distance) slope.")

json.dump({g: res1[g] for g in res1}, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "T3_B1_swampland_slopes_data.json"), "w"), default=float)
