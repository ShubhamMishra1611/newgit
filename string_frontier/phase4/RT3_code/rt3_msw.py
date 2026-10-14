"""RT3 independent re-derivation of T3 claim 7 (Martel-Shapiro-Weinberg astro-ph/9701099 + Planck 2018).

MSW model (re-derived here from the physics of their sec. 2-3, not from their closed forms):
  core volume V with overdensity delta at recombination, compensating shell U, s = V/U.
  collapse threshold (MSW eq. 4): delta_c = (729 rho_V / 500 rhobar)^(1/3)
  infall fraction (MSW eq. 7, rho_V/rhobar dropped in numerator):  F(delta) = s delta/(delta_c + s delta)
  collapsed fraction (eq. 9): F(rho_V) = int_{delta_c}^inf N(delta) F(delta) d delta, N Gaussian (eq. 21)
  flat prior => P_obs(rho_V) prop. to F(rho_V);  everything depends only on x = rho_V/(sigma^3 rhobar).
Checks: <rho_V>/(sigma^3 rhobar) (eq. 25: 2.579, 3.010, 3.224 for s=0,1,inf); quantiles (text after eq. 30:
0.0590, 1.49, 11.1 at P=0.05,0.5,0.95 for s=1); Table 1 (H0=70, n=1, Omega_B h^2 = 0.015, BBKS eq. 41,
Gamma eq. 43, COBE delta_H eq. 49 with a=-0.95, b=-0.169, f(lambda0) eq. 48, z_rec = 1000).
Planck 2018 update: sigma_0(R_G) with (i) CAMB (own call; inputs from 1807.06209 Table 2,
TT,TE,EE+lowE+lensing) and (ii) Eisenstein-Hu 1998 no-wiggle transfer (own code) normalised to sigma8.
"""
import math, sys
import numpy as np
from scipy import integrate, optimize, special

# ---------------- MSW probability distribution ----------------
def Fcol(x, s):
    """collapsed fraction (unnormalised), x = rho_V/(sigma^3 rhobar); u = delta/sigma"""
    uc = (729.0 * x / 500.0) ** (1.0 / 3.0)
    if np.isinf(s):
        return 0.5 * math.sqrt(2 * math.pi) * special.erfc(uc / math.sqrt(2))
    if s == 0:   # limit s->0, divided by s
        return math.exp(-uc * uc / 2) / uc if uc > 0 else np.inf
    f = lambda u: math.exp(-u * u / 2) * s * u / (uc + s * u)
    return integrate.quad(f, uc, np.inf, epsabs=0, epsrel=1e-11, limit=200)[0]


def make_cdf(s):
    # integrate in t = x^(1/3) to resolve small x
    ts = np.concatenate([np.linspace(0, 1, 401)[1:], np.linspace(1, 6, 1001)[1:]])
    xs = ts ** 3
    fx = np.array([Fcol(x, s) for x in xs])
    # dP = F(x) dx = F 3 t^2 dt
    g = fx * 3 * ts ** 2
    g0 = np.r_[0.0 if s != 0 else 3 * (729 / 500) ** (-1 / 3) * 0.0, g]   # integrand -> 0 at t=0 (s=0: ~ t)
    t0 = np.r_[0.0, ts]
    cum = integrate.cumulative_trapezoid(g0, t0, initial=0)
    # tail beyond t=6 negligible (exp(-uc^2/2), uc = 1.13 t)
    norm = cum[-1]
    mean = integrate.trapezoid(g0 * t0 ** 3, t0) / norm
    return t0 ** 3, cum / norm, mean


def P_le(x, s, _cache={}):
    if s not in _cache:
        _cache[s] = make_cdf(s)
    X, C, _ = _cache[s]
    return float(np.interp(x ** (1 / 3), X ** (1 / 3), C))


def quantile(p, s):
    X, C, _ = make_cdf(s) if True else None
    return float(np.interp(p, C, X ** (1 / 3))) ** 3


print("== MSW eq. 25 / quantiles ==")
for s in (0, 1, np.inf):
    X, C, mean = make_cdf(s)
    q = [float(np.interp(p, C, X ** (1 / 3))) ** 3 for p in (0.05, 0.5, 0.95)]
    print(f"  s={s}: <rho_V>/sigma^3rhobar = {mean:.4f};  quantiles 5/50/95% = {q[0]:.5f} {q[1]:.4f} {q[2]:.3f}")
print("  MSW: 2.579 (s=0), 3.010 (s=1), 3.224 (s=inf); s=1 quantiles 0.0590, 1.49, 11.1")

# ---------------- MSW Table 1 ----------------
def f_growth(lam0):
    Om = 1 - lam0
    if lam0 == 0:
        return 1.0
    I = integrate.quad(lambda w: w ** (-1 / 6) * (1 + w) ** (-1.5), 0, lam0 / Om)[0]
    return 6 * lam0 ** (5 / 6) / (5 * Om ** (1 / 3)) / I


def T_bbks(q):
    return np.log(1 + 2.34 * q) / (2.34 * q) * (1 + 3.89 * q + (16.1 * q) ** 2 + (5.46 * q) ** 3 + (6.71 * q) ** 4) ** -0.25


def sigma_rec_msw(lam0, RG, h=0.7, obh2=0.015, n=1.0, zrec=1000.0):
    Om = 1 - lam0
    ob = obh2 / h ** 2
    Gam = Om * h * math.exp(-ob - ob / Om)
    dH = 1.94e-5 * Om ** (-0.785 - 0.05 * math.log(Om)) * math.exp(-0.95 * (n - 1) - 0.169 * (n - 1) ** 2)
    c100 = 2997.9
    # sigma^2 = (c100 Gam)^(n+3) dH^2 A^-2 int q^(n+2) T^2 W^2(k RG) dq,  k = q h Gam [1/Mpc]
    integrand = lambda lq: (np.exp(lq) ** (n + 3)) * T_bbks(np.exp(lq)) ** 2 * np.exp(-(np.exp(lq) * h * Gam * RG) ** 2)
    K = integrate.quad(integrand, math.log(1e-7), math.log(1e4), limit=400)[0]
    return (c100 * Gam) ** ((n + 3) / 2) * dH * math.sqrt(K) * f_growth(lam0) / (1 + zrec)


print("\n== MSW Table 1 (H0=70, n=s=1) ==")
tab = {0.1: (0.0067, 0.0005, 0.0042, 0.0019), 0.2: (0.0063, 0.0013, 0.0040, 0.0045), 0.3: (0.0059, 0.0025, 0.0038, 0.0084),
       0.4: (0.0054, 0.0049, 0.0036, 0.015), 0.5: (0.0048, 0.0097, 0.0032, 0.027), 0.6: (0.0041, 0.021, 0.0029, 0.054),
       0.7: (0.0033, 0.054, 0.0024, 0.12), 0.8: (0.0023, 0.19, 0.0017, 0.35), 0.9: (0.0011, 0.90, 0.0008, 0.98)}
print(" lam0 | RG=1: sigma (MSW)   P (MSW)  P(from MSW's printed sigma) | RG=2: sigma (MSW)  P (MSW)  P(printed sigma)")
for lam0, (s1, p1, s2, p2) in tab.items():
    row = []
    for RG, sm, pm in ((1.0, s1, p1), (2.0, s2, p2)):
        sg = sigma_rec_msw(lam0, RG)
        x = (lam0 / (1 - lam0)) / (sg * 1001) ** 3
        xp = (lam0 / (1 - lam0)) / (sm * 1001) ** 3
        row.append(f"{sg:.4f} ({sm})  {P_le(x,1):.4f} ({pm})  {P_le(xp,1):.4f}")
    print(f"  {lam0:.1f} | " + " | ".join(row))

# ---------------- Planck 2018 update ----------------
H0, obh2, och2, ns, lnAs, tau, mnu = 67.36, 0.02237, 0.1200, 0.9649, 3.044, 0.0544, 0.06
h = H0 / 100
OmL, Om = 0.6847, 0.3153         # Table 2 derived values
g0 = 1 / f_growth(OmL)            # radiation neglected, as in MSW eq. 48
print(f"\n== Planck 2018: g0 = D(a=1)/a = {g0:.5f};  rho_obs/rho_m0 = {OmL/Om:.4f};  "
      f"MSW L* filter formula R_G = h^-1/3 (obh2/0.015)^-1/3 = {h**(-1/3)*(obh2/0.015)**(-1/3):.4f} Mpc")

kk = np.logspace(-5, 2.0, 4000)  # 1/Mpc (Gaussian filter kills k>20/Mpc for R_G>=0.5)


def sig_R(Pk, k, R, kind="G"):
    x = k * R
    W = np.exp(-x * x / 2) if kind == "G" else 3 * (np.sin(x) - x * np.cos(x)) / x ** 3
    return math.sqrt(integrate.simpson(Pk * W * W * k ** 3, x=np.log(k)) / (2 * math.pi ** 2))


sig = {}
try:
    sys.path.insert(0, "/tmp/claude-0/-home-user-newgit/c6fe2909-bd06-5dcd-8945-17832f3fe595/scratchpad/pylib")
    import camb
    pars = camb.CAMBparams()
    pars.set_cosmology(H0=H0, ombh2=obh2, omch2=och2, mnu=mnu, num_massive_neutrinos=1, tau=tau)
    pars.InitPower.set_params(As=math.exp(lnAs) * 1e-10, ns=ns)
    pars.set_matter_power(redshifts=[0.0], kmax=100.0, nonlinear=False)
    res = camb.get_results(pars)
    s8 = res.get_sigma8_0()
    PK = res.get_matter_power_interpolator(nonlinear=False, hubble_units=False, k_hunit=False,
                                           var1="delta_nonu", var2="delta_nonu")
    Pcb = PK.P(0.0, kk)
    PKt = res.get_matter_power_interpolator(nonlinear=False, hubble_units=False, k_hunit=False,
                                            var1="delta_tot", var2="delta_tot")
    Ptot = PKt.P(0.0, kk)
    print(f"  CAMB {camb.__version__}: sigma8(total) = {s8:.4f} (Planck 0.8111);  sigma8 from my integral (tot) = "
          f"{sig_R(Ptot, kk, 8 / h, 'TH'):.4f}, (cb) = {sig_R(Pcb, kk, 8 / h, 'TH'):.4f}")
    sig["CAMB_cb"] = {R: sig_R(Pcb, kk, R) for R in (0.5, 1.0, 2.0, 3.0)}
    sig["CAMB_tot"] = {R: sig_R(Ptot, kk, R) for R in (0.5, 1.0, 2.0, 3.0)}
except Exception as ex:
    print("  CAMB unavailable:", ex)

# Eisenstein-Hu 1998 no-wiggle transfer (eqs. 26-31), normalised to sigma8 = 0.8111
Tcmb = 2.7255
th = Tcmb / 2.7
omh2 = obh2 + och2 + mnu / 93.14
fb = obh2 / omh2
s_eh = 44.5 * math.log(9.83 / omh2) / math.sqrt(1 + 10 * obh2 ** 0.75)
aG = 1 - 0.328 * math.log(431 * omh2) * fb + 0.38 * math.log(22.3 * omh2) * fb ** 2
Geff = (omh2 / h) * (aG + (1 - aG) / (1 + (0.43 * kk * s_eh) ** 4))   # = Omega_m h (...)
q = kk / h * th ** 2 / Geff
L0 = np.log(2 * math.e + 1.8 * q)
C0 = 14.2 + 731 / (1 + 62.5 * q)
T0 = L0 / (L0 + C0 * q * q)
Peh = kk ** ns * T0 ** 2
Peh *= (0.8111 / sig_R(Peh, kk, 8 / h, "TH")) ** 2
sig["EH98nw"] = {R: sig_R(Peh, kk, R) for R in (0.5, 1.0, 2.0, 3.0)}
# MSW's own BBKS (eq. 41) with Gamma (eq. 43), normalised to sigma8
ob = obh2 / h ** 2
Gam = Om * h * math.exp(-ob - ob / Om)
Pbb = kk ** ns * T_bbks(kk / (h * Gam)) ** 2
Pbb *= (0.8111 / sig_R(Pbb, kk, 8 / h, "TH")) ** 2
sig["BBKS_s8"] = {R: sig_R(Pbb, kk, R) for R in (0.5, 1.0, 2.0, 3.0)}

print("\n R_G[Mpc]  M_f(total matter)[Msun] | sigma_0 by transfer fn | x_obs = rho_obs/sigma^3rhobar | P(<=rho_obs) s=1 (s=0, s=inf) | median rho_V/rho_obs")
rhocrit = 2.775e11 * h * h
for R in (0.5, 1.0, 2.0, 3.0):
    Mf = (2 * math.pi) ** 1.5 * Om * rhocrit * R ** 3
    for key in sig:
        s0 = sig[key][R]
        x = (OmL / Om) / (s0 / g0) ** 3
        P1, P0, Pi = P_le(x, 1), P_le(x, 0), P_le(x, np.inf)
        med = quantile(0.5, 1) / x
        print(f"  {R:4.1f}  {Mf:9.2e} | {key:9s} {s0:.4f} | {x:.4f} | {P1:.4f} ({P0:.4f}, {Pi:.4f}) | {med:.1f}")
print("T3: R_G=0.5/1/2/3: sigma0 2.792/2.106/1.510/1.208, x 0.049/0.114/0.309/0.603, P 0.043/0.086/0.184/0.293, median 30.4/13.0/4.8/2.5")

# ---------------- prior tilt rho_V^alpha (T3: alpha=-0.05/0/+0.05 -> P = 0.102/0.086/0.072) ----------------
ts = np.concatenate([np.linspace(0, 1, 2001)[1:], np.linspace(1, 6, 2001)[1:]])
xs = ts ** 3
Fx = np.array([Fcol(x, 1) for x in xs])
xobs = (OmL / Om) / (sig["CAMB_cb"][1.0] / g0) ** 3
for alpha in (-0.05, -0.01, 0.0, 0.01, 0.05):
    gg = Fx * xs ** alpha * 3 * ts ** 2
    cum = integrate.cumulative_trapezoid(np.r_[0, gg], np.r_[0, ts], initial=0)
    print(f"  prior rho^alpha, alpha={alpha:+.2f}: P(<=rho_obs) = {np.interp(xobs ** (1/3), np.r_[0, ts], cum) / cum[-1]:.4f}")
