"""T3 (A3): Martel-Shapiro-Weinberg (MSW) anthropic distribution of rho_V.
Reproduce MSW (astro-ph/9701099 = ApJ 492, 29 (1998)), then update with Planck 2018,
then vary the prior (power law; discrete Poisson 'landscape' prior from A1; gappy comb).

Run:  PYTHONPATH=<dir containing camb> python3 T3_A3_msw_anthropic.py > T3_A3_msw_anthropic.out
(CAMB 2.0.4 is used for the Planck linear power spectrum; if absent, only the BBKS/MSW
transfer function is used.)

MSW equations used (equation numbers as printed in the ApJ version read this session,
sources/MSW_astro-ph_9701099_ApJ.txt):
  (4)  collapse iff delta >= (729 rho_V / 500 rhobar)^{1/3}          [k=0, Weinberg 1987]
  (7)  F(delta,rho_V) = s delta [1 + c]/[c + s delta],  c = (729 rho_V/500 rhobar)^{1/3}, s = V/U
       (MSW then drop c in the numerator, text before (9))
  (9)  F(rho_V) = int_c^inf d delta N(delta) s delta/(c + s delta)
  (3)  P_obs(rho_V) = F(rho_V)/int_0^inf F            (flat prior, rho_V >= 0)
  (13) I_0(s) = int_0^1 x^2/(x+s) dx;  (17) I_n(s) = int_0^1 x^{3n+2}/(x+s) dx
  (21) N(delta) = (1/sigma) sqrt(2/pi) exp(-delta^2/2sigma^2)  (positive fluctuations)
  (22) <delta^N> = 2^{N/2} sigma^N Gamma((N+1)/2)/sqrt(pi)
  (16),(24),(25) <rho_V> = sigma^3 rhobar x {2.579, 3.010, 3.224} for s = {0, 1, inf}
  (29) closed form for P(<= rho_V)  [checked numerically against (9)]
  text: for s=1, P(<=rho_V) = 0.05, 0.5, 0.95 at rho_V = 0.0590, 1.49, 11.1 sigma^3 rhobar
  (33) sigma^2 = (1/2pi^2) int P(k) W^2(kR) k^2 dk,  (36) Gaussian window W = exp(-u^2/2)
  (39) P(k,z) = 2 pi^2 (c/H0)^{3+n} delta_H^2 k^n T^2(q) A^{-2}(z,0)
  (41) T(q) = ln(1+2.34q)/(2.34q) [1+3.89q+(16.1q)^2+(5.46q)^3+(6.71q)^4]^{-1/4}
  (42) q = k/(h Gamma Mpc^-1),  (43) Gamma = Omega0 h exp(-Omega_B0 - Omega_B0/Omega0)
  (44),(48) f(lambda0) = (6 lambda0^{5/6}/(5 Omega0^{1/3})) [int_0^{lambda0/Omega0} dw w^{-1/6}(1+w)^{-3/2}]^{-1}
  (49) delta_H = 1.94e-5 Omega0^{-0.785-0.05 ln Omega0} exp[a(n-1)+b(n-1)^2]  (a=-0.95,b=-0.169)
  (50)-(52) sigma(z_rec) = sigma_0 f(lambda0)/(1+z_rec);  sigma^3 rhobar = rho_0 [f sigma_0]^3
  Table 1 (H0=70, n=s=1, R_G=1 and 2 Mpc): sigma and P(<= rho_V*) for lambda0 = 0.1..0.9.
  Omega_B0 = 0.015 h^-2 (MSW text after (43)).
"""
import math, sys, os, json
import numpy as np
from scipy import integrate, special, optimize

np.set_printoptions(linewidth=150)
OUT = {}


def hdr(s):
    print("\n" + "=" * 90 + "\n" + s + "\n" + "=" * 90, flush=True)


# ---------------------------------------------------------------- MSW section 3 (analytic)
def I_n(n, s):
    if s == 0:
        return 1.0 / (3 * n + 2)
    if math.isinf(s):
        return 0.0  # handled via s*I_n
    return integrate.quad(lambda x: x ** (3 * n + 2) / (x + s), 0, 1, epsabs=1e-14, epsrel=1e-13)[0]


def sI_n(n, s):
    """s*I_n(s), finite for s->inf: -> 1/(3n+3)."""
    if math.isinf(s):
        return 1.0 / (3 * n + 3)
    return s * I_n(n, s)


def mean_rhoV_over_sig3rho(s):
    # (16)+(22): <rho_V> = (500/729) rhobar <d^6>/<d^3> I_1/I_0
    ratio = 2 ** 1.5 * special.gamma(3.5) / special.gamma(2.0)
    if s == 0:
        r = (1 / 5) / (1 / 2)
    elif math.isinf(s):
        r = (1 / 6) / (1 / 3)
    else:
        r = I_n(1, s) / I_n(0, s)
    return 500 / 729 * ratio * r


def F_unnorm(y, s):
    """F(rho_V)/s for y = rho_V/(sigma^3 rhobar), sigma=1 units: c = (729 y/500)^{1/3}.
    F/s = int_c^inf N(d) d/(c + s d) dd  (s=0 limit: int N d / c)."""
    c = (729.0 * y / 500.0) ** (1 / 3)
    Nd = lambda d: math.sqrt(2 / math.pi) * math.exp(-d * d / 2)
    if s == 0:
        g = lambda d: Nd(d) * d / c
    elif math.isinf(s):
        g = lambda d: Nd(d)          # s*d/(c+s*d) -> 1 ; we return F itself (not /s) here
    else:
        g = lambda d: Nd(d) * d / (c + s * d)
    return integrate.quad(g, c, np.inf, epsabs=1e-15, epsrel=1e-12)[0]


def norm_const(s):
    # int_0^inf F dy  = (500/243) <d^3> s I_0(s)   (eq. 12, sigma=1, rhobar=1); divide by s for F/s
    d3 = 2 ** 1.5 * special.gamma(2.0) / math.sqrt(math.pi)
    if s == 0:
        return 500 / 243 * d3 * I_n(0, 0)
    if math.isinf(s):
        return 500 / 243 * d3 * sI_n(0, s)
    return 500 / 243 * d3 * I_n(0, s)


class MSWdist:
    """P_obs(y) dy with y = rho_V/(sigma^3 rhobar), for shape parameter s and prior y^alpha."""

    _cache = {}

    def __init__(self, s=1.0, alpha=0.0, ygrid=None):
        self.s = s
        self.alpha = alpha
        self.y = np.geomspace(1e-9, 400, 6000) if ygrid is None else ygrid
        if s not in MSWdist._cache:
            MSWdist._cache[s] = np.array([F_unnorm(v, s) for v in self.y])
        self.F = MSWdist._cache[s]
        w = self.F * self.y ** alpha
        # cumulative via trapezoid in log y (integrand w*y)
        ly = np.log(self.y)
        dens = w * self.y
        cum = np.concatenate([[0.0], np.cumsum(0.5 * (dens[1:] + dens[:-1]) * np.diff(ly))])
        # add analytic piece below y_min: F ~ const there -> int_0^ymin y^alpha dy * F(ymin)
        cum = cum + w[0] * self.y[0] / (alpha + 1)
        self.cum = cum / cum[-1]

    def P_le(self, y):
        return float(np.interp(np.log(y), np.log(self.y), self.cum))

    def quantile(self, p):
        return float(np.exp(np.interp(p, self.cum, np.log(self.y))))


def P_le_closed_form(y, s):
    """MSW eq. (29) as decoded: P(<=rho_V) = 1-(1+b)e^{-b} + 1/(2 I_0(s)) int_b^inf e^{-x}
       [ -2s sqrt(b x) + b + 2 s^2 x ln( (1/s) sqrt(b/x) + 1 ) ] dx,  b = (1/2)(729 y/500)^{2/3}."""
    b = 0.5 * (729 * y / 500) ** (2 / 3)
    I0 = I_n(0, s)
    g = lambda x: math.exp(-x) * (-2 * s * math.sqrt(b * x) + b + 2 * s * s * x * math.log(math.sqrt(b / x) / s + 1))
    val = integrate.quad(g, b, np.inf, epsabs=1e-14, epsrel=1e-12)[0]
    return 1 - (1 + b) * math.exp(-b) + val / (2 * I0)


hdr("1. MSW sec. 3 reproduction (flat prior)")
for s in (0.0, 1.0, math.inf):
    print(f"  s={s}: <rho_V>/(sigma^3 rhobar) = {mean_rhoV_over_sig3rho(s):.3f}   (MSW eq. 25: s=0 2.579, s=1 3.010, s=inf 3.224)")
print(f"  I_1(1)/I_0(1) = {I_n(1,1)/I_n(0,1):.4f}  (MSW eq. 20: 0.467);  closed form 47/60-ln2 over 1/2... "
      f"I0(1) = {I_n(0,1):.6f} vs 1/2-1+ln2 = {0.5-1+math.log(2):.6f}")
D1 = MSWdist(1.0)
Dinf = MSWdist(math.inf)
D0 = MSWdist(0.0)
for p in (0.05, 0.5, 0.95):
    print(f"  s=1: P(<=rho_V) = {p:.2f} at rho_V/(sigma^3 rhobar) = {D1.quantile(p):.4g}   "
          f"(MSW text: 0.0590, 1.49, 11.1)   [s=0: {D0.quantile(p):.4g}, s=inf: {Dinf.quantile(p):.4g}]")
print(f"  median/<rho_V> (s=1) = {D1.quantile(0.5)/mean_rhoV_over_sig3rho(1.0):.3f}  (MSW: 0.49)")
for y in (0.01, 0.0646, 0.5, 3.0):
    print(f"  check eq.(29) closed form vs direct integration of (9): y={y}: {P_le_closed_form(y,1.0):.6f} vs {D1.P_le(y):.6f}")

# ---------------------------------------------------------------- MSW Table 1 from their sigma values
hdr("2. MSW Table 1: P(<= rho_V*) from the tabulated sigma (z_rec = 1000, s = 1)")
table1 = [  # lambda0, rhoV*/rho0, sigma(R_G=1), P(R_G=1), sigma(R_G=2), P(R_G=2)
    (0.1, 0.11, 0.0067, 0.0005, 0.0042, 0.0019),
    (0.2, 0.25, 0.0063, 0.0013, 0.0040, 0.0045),
    (0.3, 0.43, 0.0059, 0.0025, 0.0038, 0.0084),
    (0.4, 0.67, 0.0054, 0.0049, 0.0036, 0.015),
    (0.5, 1.00, 0.0048, 0.0097, 0.0032, 0.027),
    (0.6, 1.50, 0.0041, 0.021, 0.0029, 0.054),
    (0.7, 2.33, 0.0033, 0.054, 0.0024, 0.12),
    (0.8, 4.00, 0.0023, 0.19, 0.0017, 0.35),
    (0.9, 9.00, 0.0011, 0.90, 0.0008, 0.98),
]
zrec = 1000.0
print("  lambda0 | P(R_G=1) ours  MSW   | P(R_G=2) ours  MSW")
for l0, r, s1, p1, s2, p2 in table1:
    y1 = r / (s1 ** 3 * (1 + zrec) ** 3)
    y2 = r / (s2 ** 3 * (1 + zrec) ** 3)
    print(f"   {l0:.1f}    |   {D1.P_le(y1):.4f}   {p1:.4f}  |   {D1.P_le(y2):.4f}   {p2:.4f}")


# ---------------------------------------------------------------- MSW sec. 4: sigma from COBE + BBKS
def growth_f_MSW(l0):
    Om = 1 - l0
    if l0 == 0:
        return 1.0
    Iw = integrate.quad(lambda w: w ** (-1 / 6) * (1 + w) ** (-1.5), 0, l0 / Om)[0]
    return 6 * l0 ** (5 / 6) / (5 * Om ** (1 / 3)) / Iw


def T_bbks(q):
    return np.log(1 + 2.34 * q) / (2.34 * q) * (1 + 3.89 * q + (16.1 * q) ** 2 + (5.46 * q) ** 3 + (6.71 * q) ** 4) ** -0.25


def sigma0_MSW(l0, h, RG, n=1.0, OmB_h2=0.015, deltaH=None, Gamma=None):
    """linear sigma today (z=0, with Lambda growth) in MSW eqs (33),(36),(39)-(43),(49)."""
    Om = 1 - l0
    OmB = OmB_h2 / h ** 2
    if Gamma is None:
        Gamma = Om * h * math.exp(-OmB - OmB / Om)
    if deltaH is None:
        a, b = -0.95, -0.169
        deltaH = 1.94e-5 * Om ** (-0.785 - 0.05 * math.log(Om)) * math.exp(a * (n - 1) + b * (n - 1) ** 2)
    cH = 2997.9 / h  # Mpc
    lk = np.linspace(math.log(1e-5), math.log(50.0 / RG), 6000)
    k = np.exp(lk)
    q = k / (h * Gamma)
    Delta2 = deltaH ** 2 * (cH * k) ** (3 + n) * T_bbks(q) ** 2   # = k^3 P/(2 pi^2)
    integrand = Delta2 * np.exp(-(k * RG) ** 2)
    return math.sqrt(np.trapezoid(integrand, lk))


hdr("3. MSW sec. 4 reproduction: sigma(z_rec=1000) from COBE-normalized BBKS CDM (h=0.7, n=1)")
print("  lambda0 | f(lambda0) | sigma_rec R_G=1: ours  MSW | R_G=2: ours  MSW | P(<=rho*) R_G=1 ours MSW | R_G=2 ours MSW")
rows3 = []
for l0, r, s1, p1, s2, p2 in table1:
    f = growth_f_MSW(l0)
    sr1 = sigma0_MSW(l0, 0.7, 1.0) * f / (1 + zrec)
    sr2 = sigma0_MSW(l0, 0.7, 2.0) * f / (1 + zrec)
    P1 = D1.P_le(r / (sr1 ** 3 * (1 + zrec) ** 3))
    P2 = D1.P_le(r / (sr2 ** 3 * (1 + zrec) ** 3))
    rows3.append((l0, f, sr1, sr2, P1, P2))
    print(f"   {l0:.1f}    |  {f:.4f}    |      {sr1:.5f}  {s1:.4f}  |  {sr2:.5f}  {s2:.4f} |"
          f"      {P1:.4f}  {p1:.4f}      |   {P2:.4f}  {p2:.4f}")
OUT["MSW_repro"] = rows3


# ---------------------------------------------------------------- Planck 2018 update
hdr("4. Planck 2018 update (TT,TE,EE+lowE+lensing, 1807.06209 Table 2)")
P18 = dict(ombh2=0.02237, omch2=0.1200, H0=67.36, ns=0.9649, lnAs=3.044, tau=0.0544,
           sigma8=0.8111, Om=0.3153, OL=0.6847, mnu=0.06)
print("  inputs:", P18)
h = P18["H0"] / 100
rho_ratio_obs = P18["OL"] / P18["Om"]   # rho_V*/rho_m0
print(f"  rho_V*/rho_m0 = Omega_L/Omega_m = {rho_ratio_obs:.4f}")


def growth_g0(Om, OL):
    """D(a=1)/a_early for flat LCDM (radiation neglected after recombination): MSW 1/f."""
    # D(a) = (5/2) Om H(a) int_0^a da'/(a' H)^3   (H in units of H0)
    Hf = lambda a: math.sqrt(Om / a ** 3 + OL)
    I = integrate.quad(lambda a: 1 / (a * Hf(a)) ** 3, 0, 1, epsabs=1e-13)[0]
    return 2.5 * Om * Hf(1.0) * I


g0 = growth_g0(P18["Om"], P18["OL"])
print(f"  linear growth suppression today g0 = D(1)/a = {g0:.5f}  -> MSW f = 1/g0 = {1/g0:.5f};"
      f"  MSW eq.(48) with lambda0=0.6847: {growth_f_MSW(P18['OL']):.5f}")

camb_ok = True
try:
    import camb
except Exception as e:
    camb_ok = False
    print("  CAMB not available:", e)

sig_tab = {}
RGs = [0.5, 1.0, 2.0, 3.0]
if camb_ok:
    pars = camb.set_params(H0=P18["H0"], ombh2=P18["ombh2"], omch2=P18["omch2"], mnu=P18["mnu"],
                           tau=P18["tau"], As=math.exp(P18["lnAs"]) * 1e-10, ns=P18["ns"],
                           WantTransfer=True, kmax=200.0)
    pars.set_matter_power(redshifts=[0.0], kmax=200.0, nonlinear=False)
    res = camb.get_results(pars)
    s8 = res.get_sigma8_0()
    print(f"  CAMB {camb.__version__}: sigma8(z=0) from A_s = {s8:.4f}  (Planck table: 0.8111);"
          f" Omega_m(CAMB) = {res.get_Omega('cdm')+res.get_Omega('baryon')+res.get_Omega('nu'):.4f}")
    for var, lab in ((("delta_nonu", "delta_nonu"), "CDM+baryons"), (("delta_tot", "delta_tot"), "total matter")):
        PK = res.get_matter_power_interpolator(nonlinear=False, hubble_units=False, k_hunit=False,
                                               var1=var[0], var2=var[1])
        lk = np.linspace(math.log(1e-5), math.log(150.0), 20000)
        k = np.exp(lk)
        Pk = PK.P(0.0, k)
        for RG in RGs:
            integrand = k ** 3 * Pk / (2 * math.pi ** 2) * np.exp(-(k * RG) ** 2)
            sig = math.sqrt(np.trapezoid(integrand, lk))
            sig_tab[(lab, RG)] = sig
        # top-hat 8/h check
        R8 = 8 / h
        u = k * R8
        W = 3 * (np.sin(u) - u * np.cos(u)) / u ** 3
        s8chk = math.sqrt(np.trapezoid(k ** 3 * Pk / (2 * math.pi ** 2) * W ** 2, lk))
        print(f"    [{lab}] sigma(8/h Mpc, top-hat) from our integration = {s8chk:.4f}")

# BBKS (MSW eq. 41-43) shape with Planck parameters, normalized to sigma8
OmB = P18["ombh2"] / h ** 2
Gam = P18["Om"] * h * math.exp(-OmB - OmB / P18["Om"])
lk = np.linspace(math.log(1e-5), math.log(150.0), 20000)
k = np.exp(lk)
Delta2_shape = (k / h) ** (3 + P18["ns"]) * T_bbks(k / (h * Gam)) ** 2
u = k * 8 / h
W = 3 * (np.sin(u) - u * np.cos(u)) / u ** 3
norm = P18["sigma8"] ** 2 / np.trapezoid(Delta2_shape * W ** 2, lk)
for RG in RGs:
    sig_tab[("BBKS(MSW eq.41) sigma8-normalized", RG)] = math.sqrt(norm * np.trapezoid(Delta2_shape * np.exp(-(k * RG) ** 2), lk))

rho_crit = 2.775e11 * h * h   # Msun / Mpc^3
rho_m0 = P18["Om"] * rho_crit
print("\n  R_G[Mpc]  M_f=(2pi)^1.5 rho_m0 R_G^3 [Msun] | transfer        | sigma_0(R_G)  A=sigma0/g0 | x*=rho_obs/(s^3 rb) |"
      " P(<=rho_obs) s=1  (s=0, s=inf) | median rho_V/rho_obs (s=1) | mean/rho_obs")
OUT["planck"] = []
for (lab, RG), sig in sorted(sig_tab.items(), key=lambda t: (t[0][1], t[0][0])):
    A = sig / g0                    # = (1+z) sigma at early matter-dominated times = MSW's f*sigma_0
    sig3rho_over_rhom0 = A ** 3     # sigma^3 rhobar / rho_m0
    xstar = rho_ratio_obs / sig3rho_over_rhom0
    P1 = D1.P_le(xstar); P0 = D0.P_le(xstar); Pi = Dinf.P_le(xstar)
    med = D1.quantile(0.5) / xstar
    mean = mean_rhoV_over_sig3rho(1.0) / xstar
    M = (2 * math.pi) ** 1.5 * rho_m0 * RG ** 3
    print(f"   {RG:4.1f}     {M:10.3e}                     | {lab[:15]:15s} |  {sig:.4f}      {A:.4f}    |   {xstar:.4f}           |"
          f"   {P1:.4f}   ({P0:.4f}, {Pi:.4f})       |    {med:8.2f}                 | {mean:7.2f}")
    OUT["planck"].append(dict(transfer=lab, RG=RG, M=M, sigma0=sig, A=A, xstar=xstar, P_s1=P1, P_s0=P0, P_sinf=Pi, median_ratio=med))
RG_msw = h ** (-1 / 3) * (P18["ombh2"] / 0.015) ** (-1 / 3)
print(f"\n  MSW's L*-galaxy filter R_G ~ h^(-1/3) (Omega_B0 h^2/0.015)^(-1/3) Mpc with Planck values = {RG_msw:.3f} Mpc;"
      f"  MSW's 'volume per L* galaxy' R_G ~ 2 h^-1 Mpc = {2/h:.2f} Mpc")

# ---------------------------------------------------------------- priors
hdr("5. Prior dependence at the Planck baseline (CAMB CDM+baryons, R_G = 1 Mpc, s = 1)")
key = ("CDM+baryons", 1.0) if camb_ok else ("BBKS(MSW eq.41) sigma8-normalized", 1.0)
A = sig_tab[key] / g0
xstar = rho_ratio_obs / A ** 3
print(f"  baseline x* = rho_obs/(sigma^3 rhobar) = {xstar:.4f}")
print("  (5a) power-law prior p(rho_V) ~ rho_V^alpha near 0:")
for alpha in (-0.5, -0.2, -0.05, -0.01, 0.0, 0.01, 0.05, 0.2, 0.5, 1.0):
    Dal = MSWdist(1.0, alpha)
    print(f"     alpha={alpha:+5.2f}: P(<=rho_obs) = {Dal.P_le(xstar):.4f};  median rho_V/rho_obs = {Dal.quantile(0.5)/xstar:7.2f}")

print("\n  (5b) discrete Poisson ('BP discretuum') prior with mean density matching A1 (flat on average):")
print("       N_w = expected # vacua with 0 <= rho_V <= 10 sigma^3 rhobar.  P computed as")
print("       sum_{rho_i<=rho_obs} F(rho_i) / sum_i F(rho_i) for each realisation (vacua weighted by F).")
rng = np.random.default_rng(2026)
ygrid = D1.y
Fgrid = D1.F
Fint = lambda yy: np.interp(np.log(yy), np.log(ygrid), Fgrid)
for Nw in (3, 10, 30, 100, 1000, 100000):
    lam_density = Nw / 10.0           # per unit y
    ymax = 200.0
    vals = []
    nobelow = 0
    for r in range(4000 if Nw < 1000 else 400):
        n = rng.poisson(lam_density * ymax)
        if n == 0:
            continue
        yy = rng.uniform(0, ymax, n)
        yy = yy[yy > 1e-9]
        Fw = Fint(yy)
        tot = Fw.sum()
        if tot <= 0:
            continue
        below = Fw[yy <= xstar].sum()
        if below == 0:
            nobelow += 1
        vals.append(below / tot)
    vals = np.array(vals)
    print(f"     N_w={Nw:>6}: P(<=rho_obs) mean {vals.mean():.4f}, median {np.median(vals):.4f}, 16-84% [{np.quantile(vals,0.16):.4f},"
          f" {np.quantile(vals,0.84):.4f}];  fraction of landscapes with NO vacuum below rho_obs: {nobelow/len(vals):.3f}")

print("\n  (5c) gappy 'comb' prior (near-commensurate charges, A1 part d1): vacua only at y = (phi + j) * Delta, phi ~ U[0,1)")
for Delta in (0.03, 0.3, 1.0, 3.0):
    Ps = []
    for phi in np.linspace(0.0005, 0.9995, 400):
        yy = (phi + np.arange(0, int(200 / Delta) + 1)) * Delta
        yy = yy[(yy > 1e-9) & (yy < 200)]
        Fw = Fint(yy)
        Ps.append(Fw[yy <= xstar].sum() / Fw.sum())
    Ps = np.array(Ps)
    print(f"     Delta={Delta:5.2f} sigma^3 rhobar: P(<=rho_obs) ranges [{Ps.min():.4f}, {Ps.max():.4f}], phase-average {Ps.mean():.4f}")

json.dump(OUT, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "T3_A3_msw_anthropic_data.json"), "w"), default=float)
