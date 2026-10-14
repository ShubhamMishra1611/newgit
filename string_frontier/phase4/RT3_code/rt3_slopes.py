"""RT3 independent check of T3 claims 8-9 (type II dilaton-volume exponent vectors, hull minima).

(1) Minisuperspace reduction (sympy) of the 10d string-frame action
      S = int sqrt(-G) e^{-2 phi} (R + 4 (d phi)^2)
    on ds^2 = -N^2 dt^2 + a^2 dx_3^2 + rho dy_6^2 (rho = string-frame L^2), then Weyl rescaling to 4d
    Einstein frame with tau = e^{-phi} rho^{3/2}:  check K = (3/4)(dln rho)^2 + (dln tau)^2 in units of M_P^2
    relative to (1/2) R  [HKTT 0711.2512 eq. 2.5].
(2) Exponent vectors: each source's Einstein-frame potential obtained by explicit substitution of the
    string-frame scalings (sqrt(G10) e^{-chi phi} x internal scaling) -> V_E = rho^a tau^b;
    k = -grad ln V in canonical (rho_hat, tau_hat).
(3) dist(0, conv{k_i}) exactly (sympy) for source sets; physical sub-sets; cone of weak coupling/large volume.
"""
import itertools, math
import sympy as sp
import numpy as np

# ---------------- (1) kinetic normalisation ----------------
t = sp.symbols('t')
N = sp.Function('N')(t); a = sp.Function('a')(t); r = sp.Function('r')(t); ph = sp.Function('ph')(t)
coords = sp.symbols('t x1 x2 x3 y1 y2 y3 y4 y5 y6')
g = sp.diag(-N**2, a**2, a**2, a**2, r, r, r, r, r, r)
ginv = g.inv()
D = 10
def d(f, mu):
    return sp.diff(f, t) if mu == 0 else 0
Gam = [[[sum(ginv[l, s] * (d(g[s, m], n) + d(g[s, n], m) - d(g[m, n], s)) for s in range(D)) / 2
         for n in range(D)] for m in range(D)] for l in range(D)]
def Ric(m, n):
    res = 0
    for l in range(D):
        res += d(Gam[l][m][n], l) - d(Gam[l][m][l], n)
        for s in range(D):
            res += Gam[l][l][s] * Gam[s][m][n] - Gam[l][n][s] * Gam[s][m][l]
    return res
Rs = sp.simplify(sum(ginv[m, m] * Ric(m, m) for m in range(D)))
sqrtg = N * a**3 * r**3
L10 = sqrtg * sp.exp(-2 * ph) * (Rs + 4 * ginv[0, 0] * sp.diff(ph, t)**2)
# Einstein frame variables: tau = e^{-ph} r^{3/2}; a = aE/tau, N = NE/tau
aE = sp.Function('aE')(t); NE = sp.Function('NE')(t); lr = sp.Function('lr')(t); lt = sp.Function('lt')(t)
subs = {r: sp.exp(lr), ph: sp.Rational(3, 2) * lr - lt}
L = L10.subs(subs).doit()
L = L.subs({a: aE * sp.exp(-lt), N: NE * sp.exp(-lt)}).doit()
L = sp.simplify(sp.expand(L))
# integrate second derivatives by parts: L = L0 + sum_q B_q qdd ; B_q must be velocity independent
qs = [aE, lr, lt, NE]
dd = {q: sp.diff(q, t, 2) for q in qs}
Leff = sp.expand(L)
for q in qs:
    B = sp.expand(Leff).coeff(dd[q])
    if B != 0:
        assert all(sp.diff(B, sp.diff(qq, t)) == 0 for qq in qs), "B depends on velocities"
        Leff = sp.expand(Leff - B * dd[q] - sp.diff(B, t) * sp.diff(q, t))   # B qdd ~ -dB/dt qdot
Leff = sp.simplify(Leff)
# expected 4d Einstein frame minisuperspace L/(norm) = -3 aE aEd^2/NE + aE^3/(2 NE) [ (3/2) lrd^2 + 2 ltd^2 ]
aEd, lrd, ltd = sp.diff(aE, t), sp.diff(lr, t), sp.diff(lt, t)
c_a = sp.simplify(Leff.coeff(aEd**2));
P = sp.Poly(sp.expand(Leff * NE / aE**3), aEd, lrd, ltd)
print("(1) reduced Einstein-frame minisuperspace Lagrangian (up to normalisation):")
print("    L =", sp.simplify(Leff))
coef = {m: sp.simplify(c) for m, c in zip(P.monoms(), P.coeffs())}
grav = coef.get((2, 0, 0)); kr = coef.get((0, 2, 0)); kt = coef.get((0, 0, 2))
print("    coefficient ratios: (lnrho)^2 / [aEdot^2 term]*(-6) ... ")
# normalise so that the gravity term matches -3 aE aEd^2/NE   (from int sqrt(-g) R/2)
norm = sp.simplify(grav * aE**2 / (-3)) if grav is not None else None
print("    normalisation factor:", norm)
# minisuperspace: a K-term  -c (d ln x)^2  in  sqrt(-g)[R/2 - c (d ln x)^2]  gives  + c aE^3 (ln x)dot^2 / NE
print("    K_rho coefficient c (expect 3/4):", sp.simplify(kr / norm), "  K_tau coefficient c (expect 1):", sp.simplify(kt / norm),
      "  cross terms:", {m: c for m, c in coef.items() if m not in [(2, 0, 0), (0, 2, 0), (0, 0, 2)]})

# ---------------- (2) exponent vectors ----------------
rho, tau = sp.symbols('rho tau', positive=True)
# string-frame data: (chi, internal scaling w of sqrt(g6)*f) with V_s4 = rho^3 * e^{-chi phi} * f(rho)
# Einstein frame: V_E = tau^{-4} * V_s4 (from g_s = tau^{-2} g_E => sqrt(g_s) = tau^{-4} sqrt(g_E))
def VE(chi, f_scaling):
    e_mphi = tau * rho**sp.Rational(-3, 2)            # e^{-phi}
    return sp.simplify(tau**-4 * rho**3 * e_mphi**chi * rho**f_scaling)
src = {}
src['R6'] = (2, -1)            # sqrt(g6) e^{-2phi} R6, R6 ~ 1/L^2 = rho^-1
src['H3'] = (2, -3)            # |H3|^2 ~ rho^-3
src['NS5'] = (2, sp.Rational(-2, 1))   # e^{-2phi} wraps 2 internal dims: rho^{2/2} * delta^4 ~ rho^-2 relative to rho^3 -> net rho^{1}; f = rho^-2
src['Q'] = (2, 1)              # non-geometric Q: |Q|^2 ~ rho^{+1}
src['Rng'] = (2, 3)            # non-geometric R: rho^{+3}
for p in range(0, 7):
    src[f'F{p}'] = (0, -p)     # |F_p|^2 ~ rho^-p
for p in range(3, 10):
    src[f'D{p}'] = (1, sp.Rational(p - 3, 2) - 3)   # DBI wraps (p-3) internal dims: rho^{(p-3)/2}, f = that / rho^3
src['Cas'] = (0, -5)           # one-loop Casimir: L^-10 = rho^-5 (T3's derivation)
src['L10'] = (0, 0)            # 10d one-loop cosmological constant (non-SUSY)
kvec = {}
for nme, (chi, w) in src.items():
    V = VE(chi, w)
    A = sp.simplify(rho * sp.diff(V, rho) / V); B = sp.simplify(tau * sp.diff(V, tau) / V)
    kvec[nme] = sp.Matrix([-A * sp.sqrt(sp.Rational(2, 3)), -B / sp.sqrt(2)])
print("\n(2) exponent vectors k = -grad ln V in canonical (rho_hat, tau_hat):")
for nme, k in kvec.items():
    print(f"   {nme:4s}: V ~ rho^{sp.simplify(-k[0]/sp.sqrt(sp.Rational(2,3)))} tau^{sp.simplify(-k[1]*sp.sqrt(2))}"
          f"  k = ({float(k[0]):+.4f}, {float(k[1]):+.4f})  |k| = {float(sp.sqrt(k.dot(k))):.4f}")

# ---------------- (3) hull minima ----------------
def seg_dist2(p, q):
    dvec = q - p
    tt = sp.Rational(0)
    den = dvec.dot(dvec)
    tt = sp.nsimplify(-(p.dot(dvec)) / den)
    tt = min(max(tt, 0), 1)
    x = p + tt * dvec
    return sp.nsimplify(sp.simplify(x.dot(x))), tt
def hull_min(names):
    best = None
    for n1 in names:
        d2 = sp.simplify(kvec[n1].dot(kvec[n1]))
        if best is None or d2 < best[0]:
            best = (d2, (n1,), 0)
    for n1, n2 in itertools.combinations(names, 2):
        d2, tt = seg_dist2(kvec[n1], kvec[n2])
        if d2 < best[0]:
            best = (d2, (n1, n2), tt)
    return best
sets = {
    'IIB geometric (T3 set: R6,H3,F1,F3,F5,D3,D5,D7,D9,NS5)': ['R6', 'H3', 'F1', 'F3', 'F5', 'D3', 'D5', 'D7', 'D9', 'NS5'],
    'IIB without D9 (IIB proper: C10 tadpole forbids net D9; positive D9-like tadpole only in non-SUSY open-string models, e.g. Sugimoto)': ['R6', 'H3', 'F1', 'F3', 'F5', 'D3', 'D5', 'D7', 'NS5'],
    'IIB without D9, + 10d one-loop Lambda (non-SUSY)': ['R6', 'H3', 'F1', 'F3', 'F5', 'D3', 'D5', 'D7', 'NS5', 'L10'],
    'IIA geometric (T3 set: R6,H3,F0..F6,D4,D6,D8,NS5)': ['R6', 'H3', 'F0', 'F2', 'F4', 'F6', 'D4', 'D6', 'D8', 'NS5'],
    'IIA without F0 and D8': ['R6', 'H3', 'F2', 'F4', 'F6', 'D4', 'D6', 'NS5'],
    'IIB + Q': ['R6', 'H3', 'F1', 'F3', 'F5', 'D3', 'D5', 'D7', 'D9', 'NS5', 'Q'],
    'fluxes + branes only, no curvature (IIA)': ['H3', 'F0', 'F2', 'F4', 'F6', 'D4', 'D6', 'D8', 'NS5'],
    'fluxes + branes only, no curvature (IIB)': ['H3', 'F1', 'F3', 'F5', 'D3', 'D5', 'D7', 'D9', 'NS5'],
}
print("\n(3) all-positive hull minima |grad V|/V >= dist(0, conv k):")
for lab, names in sets.items():
    d2, pair, tt = hull_min(names)
    print(f"   {lab}: min^2 = {d2} -> {math.sqrt(float(d2)):.4f}  at {pair} (t={tt})")

# pointwise check with random positive coefficients
rng = np.random.default_rng(0)
names = sets['IIB geometric (T3 set: R6,H3,F1,F3,F5,D3,D5,D7,D9,NS5)']
K = np.array([[float(kvec[n][0]), float(kvec[n][1])] for n in names])
mins = []
for _ in range(20000):
    A = np.exp(rng.normal(0, 3, len(names)))
    x = rng.normal(0, 5, 2)
    w = A * np.exp(-K @ x)
    gvec = -(w @ K) / w.sum()
    mins.append(np.linalg.norm(gvec))
print(f"   pointwise random test (IIB set, 2e4 points): min |grad V|/V = {min(mins):.4f} (bound {math.sqrt(18/7):.4f})")

# ---------------- mixed-sign collinear sets; weak-coupling / large-volume cone ----------------
print("\n(4) mixed-sign scaling solutions on lines carrying >=3 collinear source vectors")
print("    rule: V = e^{-f.phi} V_perp(u); a local transverse minimum with V>0 exists for suitable |Lambda_i| iff the two")
print("    extreme terms along the line are positive and a middle one negative (foot may lie outside: then the minimum is")
print("    only local/metastable, V_perp -> 0+ beyond a barrier). Two-term (+,-): the critical point with V>0 is a maximum.")
print("    flow direction = foot f; weak coupling & non-shrinking string-frame volume need f_rho>=0 and f_tau>=sqrt(3) f_rho")
can_neg = {'R6', 'Cas', 'L10', 'Q', 'Rng'} | {f'D{p}' for p in range(3, 10)}   # O-planes, positive curvature, ...
can_pos = set(kvec)                                                          # all can be positive (fluxes only positive)
names = list(kvec)
# collect lines
lines = {}
for n1, n2 in itertools.combinations(names, 2):
    p_, q_ = kvec[n1], kvec[n2]
    if sp.simplify((q_ - p_).dot(q_ - p_)) == 0:
        continue
    members = tuple(sorted(n for n in names if sp.simplify((kvec[n] - p_)[0] * (q_ - p_)[1] - (kvec[n] - p_)[1] * (q_ - p_)[0]) == 0))
    lines[members] = (p_, q_)
found = []
for members, (p_, q_) in lines.items():
    if len(members) < 3:
        continue
    dv = q_ - p_
    tt = -(p_.dot(dv)) / dv.dot(dv)
    f = p_ + tt * dv
    u = dv / sp.sqrt(dv.dot(dv))
    pos = {n: float((kvec[n] - f).dot(u)) for n in members}
    # merge degenerate points (same vector, e.g. F5 and Cas)
    svals = sorted(set(round(v, 12) for v in pos.values()))
    inside = (svals[0] < -1e-12 and svals[-1] > 1e-12)
    ext_lo = [n for n in members if abs(pos[n] - svals[0]) < 1e-9]
    ext_hi = [n for n in members if abs(pos[n] - svals[-1]) < 1e-9]
    mids = [n for n in members if svals[0] + 1e-9 < pos[n] < svals[-1] - 1e-9]
    neg_mid = [n for n in mids if n in can_neg]
    if not neg_mid:
        continue
    f2 = sp.nsimplify(sp.simplify(f.dot(f)))
    fr, ftau = float(f[0]), float(f[1])
    cone = (fr >= -1e-12) and (ftau >= math.sqrt(3) * fr - 1e-12)
    nongeo = bool(set(members) & {'Q', 'Rng'})
    found.append((float(f2), f2, members, ext_lo, neg_mid, ext_hi, fr, ftau, cone, nongeo, inside))
found.sort(key=lambda x: x[0])
for x in found[:16]:
    print(f"   |f|={math.sqrt(x[0]):.4f} (f^2={x[1]}) line {x[2]}: + {x[3]} / - {x[4]} / + {x[5]};"
          f" f=({x[6]:+.3f},{x[7]:+.3f}) cone={x[8]} nongeo={x[9]} foot_inside={x[10]}")
g_in = [x for x in found if x[8] and not x[9] and x[0] > 1e-12]
n_in = [x for x in found if x[8] and x[0] > 1e-12]
print(f"   min |f| inside cone, geometric: {math.sqrt(g_in[0][0]):.4f} {g_in[0][2]};  incl. non-geometric: {math.sqrt(n_in[0][0]):.4f} {n_in[0][2]}")
print("   lines through the origin with a valid (+,-,+) pattern (formal dS with flat direction):",
      [x[2] for x in found if x[0] < 1e-12])
