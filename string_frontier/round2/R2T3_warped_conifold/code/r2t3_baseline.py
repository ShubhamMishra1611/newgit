"""R2-T3 baseline checks.

(A) 2512.17995 (Hebecker-Lust-Schachner-Schreyer), Sec. 4 + App. E.1: for the warped Kahler potential
       K = K_k[T^A + Tbar^A + f^A(phi, phibar)] + K_cs(phi, phibar)            (4.17)/(E.1)
    with W = W(phi) independent of T, the N=1 potential collapses to
       V = e^K K^{I Jbar} D^(0)_I W conj(D^(0)_J W)                             (4.18)/(E.26)
    with K^{I Jbar} = [ (d_I dbar_J f^A) K_A + (K_cs)_{I Jbar} ]^{-1}            (E.10)
    and K^{A Jbar} = -(d_I f^A) K^{I Jbar}                                        (E.13)
    Two-moduli case: K^{M N} K_M K_N = 3 + (T/A)|d_z K_cs|^2                      (4.10)
                     K^{z zbar} = T / (T d d K_cs - 3 f_{z zbar})                 (4.12)
    We verify these to machine precision on random points, with generic f and K_cs, by building the
    full Kahler metric numerically (complex-step-free: exact symbolic derivatives with sympy).
    We also check the first-order expansion (4.14).

(B) BDGL 1809.06861 baselines: c'' = 2^(1/3)/I(0) (2.33), sqrt(g_s) M_min = (8/3) sqrt(pi c' c'') (3.3),
    g_s M^2 >= 64 pi c' c''/9 (Junghans 2201.03572 eq. 2.9), and critical points (3.2).

(C) Normalisation of the conifold-stability parameter beta (OUR ATTEMPT):
    c' (DST 0704.4001 / BDGL 2.18 quote 1.18) recomputed as 2^(2/3) * int_0^inf I(tau) l'(tau)/l'(inf) dtau,
    where l(tau) = (tau coth tau - 1)(sinh 2tau - 2tau)/(4 sinh^2 tau) is the KS dissolved-D3-charge profile
    (H3 ^ F3 = dF5 ~ l'(tau)), so that i chi_z ^ chibar_z has radial density l'(tau)/(pi l'(inf)).
    Then in the conventions of 2406.13751 [K (2.34), W = sqrt(2/pi) int G^Omega (2.48), V_antiD3 (B.10)]
    beta = 16 pi J / (N_W^2 g_s M^2),  J = int (I/I(0)) i chi^chibar = c'c''/(2 pi),  N_W^2 = 2/pi
         = 4 pi c' c'' / (g_s M^2)   (BDGL 3.2, BCPZ 2105.03370 eq. 3.9)  and NOT 4 pi^2 c'c''/(g_s M^2) (BKS 2.47).
"""
import numpy as np
import sympy as sp
from mpmath import mp, mpf, quad as mquad, coth, sinh, diff, sqrt as msqrt, pi as mpi

rng = np.random.default_rng(1)
print("=" * 78)
print("(A) 2512.17995 warped Kahler potential: no-scale identities")
print("=" * 78)


def check_model(nK, nC, Kk_kind, seed):
    """Build K(T, phi) with nK Kahler moduli, nC complex-structure moduli; return max relative errors."""
    r = np.random.default_rng(seed)
    T = sp.symbols(f't0:{nK}'); Tb = sp.symbols(f'tb0:{nK}')
    z = sp.symbols(f'z0:{nC}'); zb = sp.symbols(f'zb0:{nC}')
    # generic real f^A(z, zbar): polynomial + |z|^(2/3)-type (conifold-like) pieces
    fA = []
    for A in range(nK):
        expr = 0
        for i in range(nC):
            a1, a2, a3 = r.normal(size=3)
            expr += a1 * z[i] * zb[i] + a2 * (z[i] ** 2 + zb[i] ** 2) / 2 - abs(a3) * (z[i] * zb[i]) ** sp.Rational(1, 3)
        if nC > 1:
            b = r.normal()
            expr += b * (z[0] * zb[1] + z[1] * zb[0])
        fA.append(expr)
    calT = [T[A] + Tb[A] + fA[A] for A in range(nK)]
    if Kk_kind == 'single':
        Kk = -3 * sp.log(calT[0])
    else:  # swiss cheese V = tau_b^(3/2) - tau_s^(3/2) with tau = calT/2
        Kk = -2 * sp.log((calT[0] / 2) ** sp.Rational(3, 2) - sp.Rational(1, 10) * (calT[1] / 2) ** sp.Rational(3, 2))
    # generic K_cs: -log(Y) with Y a positive real function
    Y = 5
    for i in range(nC):
        c1 = 1 + abs(r.normal())
        Y += c1 * z[i] * zb[i] + 0.3 * r.normal() * (z[i] ** 2 + zb[i] ** 2) + 0.2 * (z[i] * zb[i]) * sp.log(z[i] * zb[i])
    Kcs = -sp.log(Y)
    K = Kk + Kcs
    hol = list(T) + list(z); ahol = list(Tb) + list(zb)
    n = len(hol)
    # superpotential depending on phi only
    coeffs = r.normal(size=(nC, 3)) + 1j * r.normal(size=(nC, 3))
    W = sp.Rational(1, 3) + sum((coeffs[i, 0] + coeffs[i, 1] * z[i] + coeffs[i, 2] * z[i] ** 2) * z[i] for i in range(nC)) \
        + sum(sp.I * 0.7 * z[i] * sp.log(z[i]) for i in range(nC))
    Wb = sp.conjugate(W).subs({sp.conjugate(z[i]): zb[i] for i in range(nC)})
    # point: choose values (treat holomorphic and antiholomorphic as independent, then set zb=conj(z))
    tv = [3 + 2 * abs(r.normal()) for _ in range(nK)]
    zv = [0.2 * (r.normal() + 1j * r.normal()) for _ in range(nC)]
    subs = {}
    for A in range(nK):
        th = r.normal()
        subs[T[A]] = tv[A] + 1j * th; subs[Tb[A]] = tv[A] - 1j * th
    for i in range(nC):
        subs[z[i]] = zv[i]; subs[zb[i]] = np.conj(zv[i])
    ev = lambda e: complex(sp.N(e.subs(subs), 30))
    Kval = ev(K).real
    dK = [ev(sp.diff(K, h)) for h in hol]
    Kmat = np.array([[ev(sp.diff(K, hol[a], ahol[b])) for b in range(n)] for a in range(n)])
    Kinv = np.linalg.inv(Kmat)  # Kinv[a,b] = K^{a bbar} in the convention (K_{a bbar})^{-1} transposed appropriately
    # Note: K_{a bbar} matrix M[a,b]; inverse metric K^{b abar} satisfies sum_b K^{a bbar} K_{c bbar} = delta;
    # we use G = inv(M).T so that sum_b G[a,b] M[c,b] = delta_ac
    G = np.linalg.inv(Kmat).T
    Wv = ev(W); dW = [ev(sp.diff(W, h)) for h in hol]
    DW = [dW[a] + dK[a] * Wv for a in range(n)]
    V_full = np.exp(Kval) * (sum(G[a, b] * DW[a] * np.conj(DW[b]) for a in range(n) for b in range(n)) - 3 * abs(Wv) ** 2)
    # 2512.17995 (4.18): V = e^K K^{I J} D0_I W conj(D0_J W), with D0 using K_cs only, K^{IJ} the cs block of full inverse
    dKcs = [ev(sp.diff(Kcs, zz)) for zz in z]
    D0 = [dW[nK + i] + dKcs[i] * Wv for i in range(nC)]
    Gcs = G[nK:, nK:]
    V_418 = np.exp(Kval) * sum(Gcs[i, j] * D0[i] * np.conj(D0[j]) for i in range(nC) for j in range(nC))
    # (E.10): K^{IJ} = [ (d_I dbar_J f^A) K_A + (K_cs)_{IJ} ]^{-1}
    KA = dK[:nK]
    E9 = np.array([[sum(ev(sp.diff(fA[A], z[i], zb[j])) * KA[A] for A in range(nK)) + ev(sp.diff(Kcs, z[i], zb[j]))
                    for j in range(nC)] for i in range(nC)])
    GE10 = np.linalg.inv(E9).T
    # (E.13): K^{A Jbar} = -(d_I f^A) K^{I Jbar}
    dfA = np.array([[ev(sp.diff(fA[A], z[i])) for i in range(nC)] for A in range(nK)])
    GE13 = -dfA @ Gcs
    err_418 = abs(V_full - V_418) / abs(V_full)
    err_E10 = np.max(np.abs(Gcs - GE10)) / np.max(np.abs(Gcs))
    err_E13 = np.max(np.abs(G[:nK, nK:] - GE13)) / np.max(np.abs(G[:nK, nK:]))
    # (E.24): K^{I Jbar} (dbar_J f^B) K_Bbar + K^{I Bbar} K_Bbar = 0
    dKb = [ev(sp.diff(K, h)) for h in ahol]
    dfAb = np.array([[ev(sp.diff(fA[A], zb[j])) for j in range(nC)] for A in range(nK)])
    E24 = np.array([sum(Gcs[i, j] * dfAb[B, j] * dKb[B] for j in range(nC) for B in range(nK))
                    + sum(G[nK + i, B] * dKb[B] for B in range(nK)) for i in range(nC)])
    err_E24 = np.max(np.abs(E24)) / max(1e-300, np.max(np.abs(Gcs)))
    out = dict(V_full=V_full.real, V_418=V_418.real, err_418=err_418, err_E10=err_E10, err_E13=err_E13, err_E24=err_E24)
    if nK == 1 and nC == 1:
        # (4.10) and (4.12)
        calTv = ev(calT[0]).real
        Kzz = ev(sp.diff(Kcs, z[0], zb[0])).real; fzz = ev(sp.diff(fA[0], z[0], zb[0])).real
        calA = calTv * Kzz - 3 * fzz
        lhs410 = sum(G[a, b] * dK[a] * np.conj(dK[b]) for a in range(n) for b in range(n)).real
        rhs410 = 3 + calTv / calA * abs(dKcs[0]) ** 2
        out['err_410'] = abs(lhs410 - rhs410) / abs(rhs410)
        out['err_412'] = abs(G[1, 1] - calTv / calA) / abs(calTv / calA)
    return out


for (nK, nC, kind, seed) in [(1, 1, 'single', 11), (1, 1, 'single', 12), (1, 2, 'single', 13), (2, 2, 'swiss', 14), (2, 1, 'swiss', 15)]:
    o = check_model(nK, nC, kind, seed)
    extra = '' if 'err_410' not in o else f"  err(4.10)={o['err_410']:.1e} err(4.12)={o['err_412']:.1e}"
    print(f"nK={nK} nC={nC} Kk={kind:6s}: V_full={o['V_full']:+.6e}  V(4.18)={o['V_418']:+.6e}  rel.err={o['err_418']:.1e}"
          f"  err(E.10)={o['err_E10']:.1e} err(E.13)={o['err_E13']:.1e} err(E.24)={o['err_E24']:.1e}{extra}")

# (4.14): first-order expansion in 1/(T+Tbar) for the two-moduli model
print("\n(4.14) first-order expansion check: V/[e^{K0} K_cs^{zz}|D0W|^2] - 1 vs -3(f - K_cs^{zz} f_zz)/(T+Tbar)")
for tt in [10, 30, 100, 300, 1000]:
    Kzz, fz, fzz = 0.7, 0.9, -0.4  # numbers: K_cs,zzbar ; f ; f_zzbar at the point
    calT = 2 * tt + fz
    exact = (calT ** -3 * calT / (calT * Kzz - 3 * fzz)) / ((2 * tt) ** -3 / Kzz)
    first = 1 - 3 * (fz - fzz / Kzz) / (2 * tt)
    print(f"  T+Tbar={2*tt:6d}: exact ratio={exact:.8f}  first-order={first:.8f}  diff*(T+Tbar)^2={(exact-first)*(2*tt)**2:.4f}")

print("\n" + "=" * 78)
print("(B) BDGL 1809.06861 baselines")
print("=" * 78)
mp.dps = 30
Iint = lambda x: (x * coth(x) - 1) / sinh(x) ** 2 * (sinh(2 * x) - 2 * x) ** (mpf(1) / 3)
I = lambda t: mquad(Iint, [t, t + 5, t + 20, mp.inf])
I0 = I(mpf('1e-30'))
cpp = 2 ** (mpf(1) / 3) / I0
cp_lit = mpf('1.18')
print(f"I(0) = {float(I0):.6f};  c'' = 2^(1/3)/I(0) = {float(cpp):.5f}   [BDGL (2.33): ~1.75]")
print(f"sqrt(g_s) M_min = (8/3) sqrt(pi c' c'') = {float(8/3*msqrt(mpi*cp_lit*cpp)):.4f}   [BDGL (3.3): 6.8]")
print(f"g_s M^2 >= 64 pi c'c''/9 = {float(64*mpi*cp_lit*cpp/9):.3f} (c'=1.18)  [Junghans (2.9): 46.1 with c''=1.75 -> {float(64*mpi*cp_lit*1.75/9):.2f}]")
# critical points of V(L) = e^{4L/3}(L^2 + beta), compare with BDGL (3.2): L = -3/4 +- sqrt(9/16 - beta)
from scipy.optimize import brentq
for gM2 in [160.0, 60.0, 46.0]:
    beta = float(4 * mpi * cp_lit * cpp / gM2)
    dV = lambda L: np.exp(4 * L / 3) * ((4 / 3) * (L ** 2 + beta) + 2 * L)
    Ls = np.linspace(-6, 3, 90001); s = np.sign(dV(Ls)); roots = [brentq(dV, Ls[k], Ls[k + 1]) for k in range(len(Ls) - 1) if s[k] != s[k + 1]]
    disc = 9 / 16 - beta
    pred = [] if disc < 0 else [-0.75 - np.sqrt(disc), -0.75 + np.sqrt(disc)]
    print(f"  g_sM^2={gM2:6.1f} beta={beta:.4f}: numerical critical L = {[round(x,6) for x in roots]}  BDGL(3.2): {[round(x,6) for x in pred]}")

print("\n" + "=" * 78)
print("(C) c' from the KS geometry and the normalisation of beta (OUR ATTEMPT)")
print("=" * 78)
ell = lambda t: (t * coth(t) - 1) * (sinh(2 * t) - 2 * t) / (4 * sinh(t) ** 2)
lp = lambda t: diff(ell, t)
lp_inf = lp(mpf(40))
cp_calc = 2 ** (mpf(2) / 3) * mquad(lambda t: I(t) * lp(t) / lp_inf, [mpf('1e-12'), 1, 3, 8, 20, 60])
J = cp_calc * cpp / (2 * mpi)
print(f"l'(inf) = {float(lp_inf):.6f}")
print(f"c' = 2^(2/3) int_0^inf I(tau) l'(tau)/l'(inf) dtau = {float(cp_calc):.5f}   [DST 0704.4001 via BDGL (2.18): c' ~ 1.18]")
print(f"J = int (I/I(0)) i chi^chibar = c'c''/(2pi) = {float(J):.5f}")
NW2 = 2 / mpi  # 2406.13751 (2.48)
print(f"beta(2406 conventions) = 16 pi J/(N_W^2 g_s M^2) = {float(16*mpi*J/NW2):.4f}/(g_sM^2);  4 pi c'c'' = {float(4*mpi*cp_calc*cpp):.4f};"
      f"  4 pi^2 c'c'' (BKS 2.47) = {float(4*mpi**2*cp_calc*cpp):.4f}")
print(f"=> bound g_s M^2 >= (16/9)*4 pi c'c'' = {float(64*mpi*cp_calc*cpp/9):.3f}  (with computed c'={float(cp_calc):.4f});"
      f"  BKS normalisation would give {float(64*mpi**2*cp_calc*cpp/9):.1f}")
# sanity: the density l'(tau)/(pi l'(inf)) integrates to (1/pi) * tau_UV + const  => log coefficient 1/pi
for tU in [10, 20, 40]:
    tot = mquad(lambda t: lp(t) / lp_inf, [mpf('1e-12'), 1, 5, tU])
    print(f"  int_0^{tU} l'/l'(inf) dtau = {float(tot):.4f}  (= tau_UV - 1 + ...: {tU-1})")
