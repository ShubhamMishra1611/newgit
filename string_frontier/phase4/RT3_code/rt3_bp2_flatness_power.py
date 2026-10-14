"""RT3 independent re-derivation of T3 claims 1, 3, 4, 5 (BP flatness, Poisson, min|Lambda|,
near-commensurate charges), with NEW seeds / radii, and a POWER test of the flatness statistic.

Flatness statistic (re-implemented from T3's description): N(Lambda in [0,eps]) and
N(Lambda in [-eps,0]) on a log grid eps_k from ~3 mean spacings to 1e-3 S; power-law fit p on
N >= 10; chi^2/dof of (n-n_pred)/sqrt(n_pred) on the disjoint shells eps_{k-1}<|Lambda|<=eps_k
(both signs), n_pred from the face-corrected smooth count.
"""
import sys, os, math, time, json
import numpy as np
from scipy import stats
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rt3_bp_core import Shell, smooth_count, smooth_rho, bp_rho, exact_X

OUT = {}


def flatness(c, S, nlev=200_000, ngrid=22, tag="", rho_ref=None, verbose=True):
    rho = smooth_rho(c, S) if rho_ref is None else rho_ref
    epsmax = 1e-3 * S
    sh = Shell(c, S + 2 * epsmax + 1e-9)
    eps = np.geomspace(3.0 / rho, epsmax, ngrid)
    Np, Nm, Pp, Pm = [], [], [], []
    for e in eps:
        Np.append(sh.count(S, S + 2 * e)); Nm.append(sh.count(S - 2 * e, S))
        Pp.append(float(smooth_count(c, S + 2 * e) - smooth_count(c, S)))
        Pm.append(float(smooth_count(c, S) - smooth_count(c, S - 2 * e)))
    Np, Nm, Pp, Pm = map(np.array, (Np, Nm, Pp, Pm))
    sel = Np >= 10
    p = np.polyfit(np.log(eps[sel]), np.log(Np[sel]), 1)[0] if sel.sum() > 2 else float("nan")
    # disjoint shells
    n = np.concatenate([np.diff(np.r_[0, Np]), np.diff(np.r_[0, Nm])])
    npred = np.concatenate([np.diff(np.r_[0, Pp]), np.diff(np.r_[0, Pm])])
    pull = (n - npred) / np.sqrt(npred)
    chi2 = float(np.sum(pull ** 2)); dof = len(n)
    pval = float(stats.chi2.sf(chi2, dof))
    # spacing statistics
    W = 0.5 * nlev / rho
    X = sh.listing(S - 2 * W, S + 2 * W)
    L = (X - S) / 2
    x = rho * L
    s = np.diff(x)
    s = s / s.mean()
    ks = stats.kstest(s, "expon")
    # number variance
    rng = np.random.default_rng(1)
    sig = {}
    for Lw in (1, 10, 100, 1000):
        if Lw < 0.2 * x[-1] - 0.2 * x[0]:
            st = rng.uniform(x[0], x[-1] - Lw, 20000)
            cnt = np.searchsorted(x, st + Lw) - np.searchsorted(x, st)
            sig[Lw] = float(cnt.var())
    res = dict(J=len(c), S=S, rho_smooth=rho, rho_bp=bp_rho(c, S), p=p, chi2_dof=chi2 / dof, chi2_p=pval,
               maxpull=float(np.abs(pull).max()), nlev=len(L), dens_ratio=(len(L) - 1) / (L[-1] - L[0]) / rho,
               var_s=float(s.var()), P_s01=float(np.mean(s < 0.1)), ks_p=float(ks.pvalue), maxgap=float(s.max()),
               sigma2=sig, Nmin=int(Np[0]), smallest_eps=float(eps[0]))
    if verbose:
        print(f"  {tag} J={len(c)} S={S:.3f} n_max={math.sqrt(S/c.min()):.0f}: rho_smooth={rho:.4e}  rho_BP(D=2^J)={res['rho_bp']:.4e}"
              f" (smooth/BP={rho/res['rho_bp']:.3f}); measured/smooth over {len(L)} levels = {res['dens_ratio']:.4f}")
        print(f"     fit p = {p:.4f}; chi2/dof = {chi2:.1f}/{dof} = {chi2/dof:.2f} (p={pval:.3g}); max|pull| = {res['maxpull']:.2f}")
        print(f"     spacings: var s = {s.var():.4f}, P(s<0.1) = {np.mean(s<0.1):.4f}, KS p = {ks.pvalue:.3f}, max s = {s.max():.1f}"
              f" (ln M = {math.log(len(s)):.1f});  Sigma^2: " + ", ".join(f"L={k}:{v:.1f}" for k, v in sig.items()), flush=True)
    return res


if __name__ == "__main__":
    T0 = time.time()
    print("== (1,3) flatness + spacing statistics, NEW seeds (RT3) ==")
    gen = []
    for J, R, seed in ((4, 800.0, 101), (6, 200.0, 102), (8, 70.0, 103), (10, 40.0, 104), (12, 26.0, 105)):
        rng = np.random.default_rng(seed)
        c = (1 + rng.random(J)) ** 2
        S = R * R + 0.3183
        gen.append(flatness(c, S, tag="generic"))
    OUT["generic"] = gen

    print("\n== POWER test: same statistic on near-commensurate charges c_i = 1 + delta u_i, J=8, S=400.3711 ==")
    S = 400.3711
    pw = []
    for dS in (120.0, 12.0, 4.0, 1.2):
        delta = dS / S
        rng = np.random.default_rng(300)
        c = 1.0 + delta * rng.random(8)
        try:
            r = flatness(c, S, nlev=100_000, tag=f"deltaS={dS:g}")
            pw.append(r)
        except Exception as ex:  # windows may be empty
            print("   failed:", ex)
    OUT["power"] = pw

    print("\n== (5) near-commensurate: gaps and min|Lambda| near Lambda=0 (J=8, S=400.3711) ==")
    nc = []
    for dS in (120.0, 12.0, 1.2, 0.12, 0.0):
        delta = dS / S
        rng = np.random.default_rng(300)
        c = 1.0 + delta * rng.random(8) if dS > 0 else np.ones(8)
        rho_gen = smooth_rho(1.0 + 1e-3 * np.random.default_rng(300).random(8), S)  # generic reference density
        sh = Shell(c, S + 2.0)
        X = sh.listing(S - 2.0, S + 2.0)
        L = (X - S) / 2
        Lu = np.unique(np.round(L, 11))
        near = int(np.sum(np.abs(L) < 50 / rho_gen))
        gaps = np.diff(Lu)
        mn = float(np.min(np.abs(Lu)))
        print(f"  delta*S={dS:7.3g}: levels in |Lam|<1: {len(L):,} (distinct {len(Lu):,}); generic expectation {2*rho_gen:.3g};"
              f" N(|Lam|<50/rho)={near} (expect 100); max gap = {gaps.max()*rho_gen:.3g} generic mean spacings;"
              f" min|Lam| = {mn:.3e} (generic Poisson mean {1/(2*rho_gen):.1e})", flush=True)
        nc.append(dict(dS=dS, nlev=len(L), ndist=len(Lu), near=near, maxgap_sp=float(gaps.max() * rho_gen), minL=mn))
    OUT["nearcomm"] = nc

    print("\n== (4) min|Lambda| statistics, RT3 seeds, exact re-check ==")
    zs = []
    for J, R, seed in ((4, 400.0, 401), (6, 80.0, 402), (8, 40.0, 403), (10, 24.0, 404), (12, 16.0, 405), (16, 10.0, 406)):
        rng = np.random.default_rng(seed)
        c = (1 + rng.random(J)) ** 2
        S0 = R * R
        sh = Shell(c, S0 + 1.0 + 1e-3)
        z = []
        for k in range(60):
            S = float(S0 + rng.random())
            rho = smooth_rho(c, S)
            w = 25.0 / rho
            X, ia, ib = sh.listing(S - w, S + w, with_n=True)
            j = int(np.argmin(np.abs(X - S)))
            # exact re-check of the 3 nearest candidates
            order = np.argsort(np.abs(X - S))[:3]
            ex = min(abs(float(exact_X(c, sh.nvec(ia[o], ib[o])) - __import__('fractions').Fraction(S))) for o in order)
            z.append(2 * rho * ex / 2)
        z = np.array(z)
        ks = stats.kstest(z, "expon")
        print(f"  J={J:2d} R={R:5.1f}: 60 shells, <z> = {z.mean():.3f} +- {z.std()/math.sqrt(len(z)):.3f}, var = {z.var():.3f}, KS p = {ks.pvalue:.3f}", flush=True)
        zs.extend(z.tolist())
    zs = np.array(zs)
    ks = stats.kstest(zs, "expon")
    print(f"  pooled N={len(zs)}: <z> = {zs.mean():.3f} +- {zs.std()/math.sqrt(len(zs)):.3f}, var = {zs.var():.3f}, KS p = {ks.pvalue:.3f}")
    OUT["minz"] = dict(mean=float(zs.mean()), se=float(zs.std() / math.sqrt(len(zs))), ks_p=float(ks.pvalue))

    json.dump(OUT, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "rt3_bp2_flatness_power.json"), "w"), indent=1, default=float)
    print(f"\nruntime {time.time()-T0:.0f}s")
