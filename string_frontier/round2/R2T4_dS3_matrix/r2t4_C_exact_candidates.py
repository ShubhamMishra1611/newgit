"""R2-T4 script C: compare exact (all-order) proposals for Z(S^3) of pure dS3 gravity, in the invariant form
of ADLS 2009.12464 eq. (1.3):  log|Z| = S0 - a log S0 + b + sum_n c_n S0^{-2n}  (S0 chosen so that NO odd powers
of 1/S0 appear; odd powers are local counterterms / field redefinitions, ADLS around (1.3), (8.16)-(8.18)).

Candidates (all functions of the invariant S0):
 [CEM-h]  CEM 2501.01486 (85) [= SciPost (5.7)], |b|<1 branch, b^2 = i beta. Natural S0 = 2 pi/beta (shown below to
          have no odd powers). CEM's own dictionary S_GH = (pi/3) Im c = 2 pi (1/beta - beta).
 [CEM-p]  same construction but with the Virasoro S-kernel in the Liouville-momentum (p) basis instead of the
          h basis (CEM say they have no better justification for the h basis than an expectation; SciPost p.34).
 [ADLS]   ADLS (8.18): S0 + log| sqrt(4/K) sin(pi/K) |^2, K = 2 + i S0/2pi  (SO(4), bare l=0 => renormalized Re K = 2).
 [CCFZ0]  Castro-Coman-Fliss-Zukowski 2302.12281 (2.82)/(2.84) with renormalized delta-hat = 0 (Re K = 0), times 2
          for SO(4) vs SU(2)xSU(2):  Z = 2 * i (2/s) e^{2 pi s} sinh^2(pi/s),  S0 = 2 pi s.
Also: the 2009.12464 one-loop truncation S0 - 3 log S0 + 5 log 2pi.

Outputs: exact series coefficients (sympy) + independent high-precision numerical fits with error estimates,
the phases, and the exact ratio CEM-h / CCFZ0.
"""
import sympy as sp
import mpmath as mp

mp.mp.dps = 60
u = sp.symbols('u', positive=True)        # u = 1/S0
PI = sp.pi

# ---------------- perturbative parts of log|Z| (e^{-S0} terms dropped), as functions of u ----------------
beta = 2*PI*u                                   # CEM: S0 = 2 pi / beta
# |Z_CEM| = 16 beta sinh^2(pi beta) sinh^2(pi/beta)/(1+beta^2)^2 ; 2 log sinh(pi/beta) = 2pi/beta - 2log2 + O(e^-S0)
L_CEMh = sp.log(16*beta) + 2*sp.log(sp.sinh(PI*beta)) + (2*PI/beta - 2*sp.log(2)) - 2*sp.log(1+beta**2)
# p basis: S11^(p) = S_{1,0}-S_{1,1} = -8 sqrt2 sin(pi b^2) sin(pi b^-2)  ->  |Z| = 128 sinh^2(pi beta) sinh^2(pi/beta)
L_CEMp = sp.log(128) + 2*sp.log(sp.sinh(PI*beta)) + (2*PI/beta - 2*sp.log(2))
# CCFZ delta-hat=0, SO(4): |Z| = 2*(2/s) e^{2 pi s} sinh^2(pi/s), s = S0/2pi = 1/(2 pi u)
s_ = 1/(2*PI*u)
L_CCFZ0 = sp.log(4/s_) + 2*PI*s_ + 2*sp.log(sp.sinh(PI/s_))
# ADLS (8.18): S0 + log 4 - (1/2) log|K|^2 + log[ sin(pi/K) sin(pi/Kbar) ],  K = 2 + i x, x = S0/2pi
x_ = 1/(2*PI*u)
K = 2 + sp.I*x_; Kb = 2 - sp.I*x_
L_ADLS = 1/u + sp.log(4) - sp.Rational(1, 2)*sp.log(4 + x_**2) + sp.log(sp.sin(PI/K)*sp.sin(PI/Kb))

def series_report(name, L, order=7):
    R = sp.simplify(L - 1/u)                      # remove S0
    # isolate the log: R = a*log(u) + regular(u)   (log S0 = -log u)
    a = sp.limit(sp.diff(R, u)*u, u, 0)            # coefficient of log u  ( = +a_ADLS since -a log S0 = a log u )
    reg = sp.simplify(R - a*sp.log(u))
    ser = sp.series(reg, u, 0, order).removeO()
    ser = sp.expand(sp.simplify(ser))
    coeffs = [sp.simplify(ser.coeff(u, k)) if k > 0 else sp.simplify(ser.subs(u, 0)) for k in range(order)]
    print(f"\n[{name}]  log|Z| = S0 - ({a}) log S0 + b + sum_k d_k S0^-k")
    print(f"   a = {a};  b = {sp.simplify(coeffs[0])}  (5 log 2pi = {sp.N(5*sp.log(2*PI), 12)}; b = {sp.N(coeffs[0], 12)})")
    for k in range(1, order):
        print(f"   d_{k} = {sp.simplify(coeffs[k])}  = {sp.N(coeffs[k], 12)}")
    return a, coeffs

print("="*100)
print("EXACT SERIES (sympy) in the invariant variable S0 (no odd powers <=> local terms absent)")
res = {}
for name, L in [("CEM-h (85)", L_CEMh), ("CEM-p (p-basis variant)", L_CEMp),
                ("CCFZ0 (2.84), SO(4)", L_CCFZ0), ("ADLS (8.18)", L_ADLS)]:
    res[name] = series_report(name, L)

# ---------------- independent numerics: exact functions incl. non-perturbative terms, fits ----------------
def Z_CEMh(S0):
    beta = 2*mp.pi/S0; b2 = mp.mpc(0, beta)
    return -16*b2*mp.sin(mp.pi*b2)**2*mp.sin(mp.pi/b2)**2/(1-b2**2)**2
def Z_CEMp(S0):
    beta = 2*mp.pi/S0; b2 = mp.mpc(0, beta)
    return 128*mp.sin(mp.pi*b2)**2*mp.sin(mp.pi/b2)**2
def Z_ADLS(S0):
    Kc = 2 + 1j*S0/(2*mp.pi)
    return mp.exp(S0)*abs(mp.sqrt(4/Kc)*mp.sin(mp.pi/Kc))**2
def Z_CCFZ0(S0):
    s = S0/(2*mp.pi)
    return 2*1j*(2/s)*mp.exp(2*mp.pi*s)*mp.sinh(mp.pi/s)**2

def fit(Zf, S0s, nterms):
    """least squares for R(S0) = log|Z| - S0 = -a log S0 + b + sum_{k=1}^{nterms} d_k S0^-k"""
    rows, rhs = [], []
    for S0 in S0s:
        S0 = mp.mpf(S0)
        rows.append([-mp.log(S0), 1] + [S0**-k for k in range(1, nterms+1)])
        rhs.append(mp.log(abs(Zf(S0))) - S0)
    A = mp.matrix(rows); y = mp.matrix(rhs)
    sol = mp.lu_solve(A.T*A, A.T*y)
    return [sol[i] for i in range(len(rows[0]))]

print("\n" + "="*100)
print("NUMERICAL FITS (exact functions incl. e^{-S0} terms), R = log|Z| - S0 = -a log S0 + b + sum d_k S0^-k")
print("error estimate = spread over fit windows S0 in [40,4000],[80,8000],[160,16000] and nterms in {6,8,10}")
for name, Zf in [("CEM-h (85)", Z_CEMh), ("CEM-p", Z_CEMp), ("CCFZ0 (2.84)xSO4", Z_CCFZ0), ("ADLS (8.18)", Z_ADLS)]:
    ests = []
    for lo in [40, 80, 160]:
        S0s = [lo*mp.power(100, mp.mpf(j)/39) for j in range(40)]
        for nt in [6, 8, 10]:
            ests.append(fit(Zf, S0s, nt))
    keys = ['a', 'b', 'd1', 'd2', 'd3']
    out = []
    for i, k in enumerate(keys):
        vals = [e[i] for e in ests]
        m = sum(vals)/len(vals); err = max(abs(v-m) for v in vals)
        out.append(f"{k}={mp.nstr(m, 12)}(+-{mp.nstr(err, 2)})")
    print(f"  {name:18s}: " + "  ".join(out))

print("\nanalytic targets: 5 log 2pi =", mp.nstr(5*mp.log(2*mp.pi), 12),
      "| c1 values: CEM-h 4pi^2(pi^2/3-2) =", mp.nstr(4*mp.pi**2*(mp.pi**2/3-2), 12),
      "; CCFZ0 4pi^4/3 =", mp.nstr(4*mp.pi**4/3, 12),
      "; ADLS 4pi^2(pi^2/3-6) =", mp.nstr(4*mp.pi**2*(mp.pi**2/3-6), 12))

# ---------------- phases ----------------
print("\nPHASES arg Z/pi at S0 = 20, 200, 2000:")
for name, Zf in [("CEM-h (85), |b|<1", Z_CEMh), ("CEM-p", Z_CEMp), ("CCFZ0 (2.82)", Z_CCFZ0), ("ADLS (8.18) [canonical framing]", Z_ADLS)]:
    print(f"  {name:32s}: " + ", ".join(mp.nstr(mp.arg(Zf(mp.mpf(S0)))/mp.pi, 10) for S0 in [20, 200, 2000]))
# |b|>1 branch of (84): S11 = 4 i b^3 sin sin/(b^4-1), b^2 = i beta, beta large; S0 = 2 pi beta
def Z_CEMh_big(S0):
    beta = S0/(2*mp.pi); b2 = mp.mpc(0, beta)
    return -16*b2**3*mp.sin(mp.pi*b2)**2*mp.sin(mp.pi/b2)**2/(b2**2-1)**2
print("  CEM-h (85), |b|>1 branch        : " + ", ".join(mp.nstr(mp.arg(Z_CEMh_big(mp.mpf(S0)))/mp.pi, 10) for S0 in [20, 200, 2000]),
      "  | |Z| equal to |b|<1 branch:", all(abs(abs(Z_CEMh_big(mp.mpf(S0)))/abs(Z_CEMh(mp.mpf(S0)))-1) < 1e-40 for S0 in [20, 200]))

# ---------------- exact ratio CEM-h / CCFZ0 ----------------
print("\nEXACT RATIO |Z_CEM-h| / |Z_CCFZ0,SO4|  vs  (1-e^{-S0})^2/(1+4pi^2/S0^2)^2 :")
for S0 in [3, 10, 50, 300]:
    S0 = mp.mpf(S0)
    r = abs(Z_CEMh(S0))/abs(Z_CCFZ0(S0)); pred = (1-mp.exp(-S0))**2/(1+4*mp.pi**2/S0**2)**2
    print(f"  S0={mp.nstr(S0,4)}: ratio={mp.nstr(r, 20)}  pred={mp.nstr(pred, 20)}")

# ---------------- CEM in its own dictionary S_GH = (pi/3) Im c ----------------
print("\nCEM-h in CEM's own dictionary S_GH = (pi/3) Im c = 2pi(1/beta - beta):")
print("  log|Z| - [S_GH - 3 log S_GH + 5 log 2pi]  vs 4pi^2/S_GH (an odd = local 'two-loop' term):")
for S0 in [100, 1000, 10000]:
    S0 = mp.mpf(S0); beta = 2*mp.pi/S0; SGH = 2*mp.pi*(1/beta - beta)
    dev = mp.log(abs(Z_CEMh(S0))) - (SGH - 3*mp.log(SGH) + 5*mp.log(2*mp.pi))
    print(f"  S_GH={mp.nstr(SGH,8)}: dev={mp.nstr(dev,12)}  4pi^2/S_GH={mp.nstr(4*mp.pi**2/SGH,12)}  dev*S_GH/4pi^2={mp.nstr(dev*SGH/(4*mp.pi**2),10)}")

# ---------------- verify (84) from (83) by analytic continuation ----------------
print("\nCHECK (84) from (83): S_{1h} = 2 sqrt2 sinh(2 pi b p) sinh(2 pi p/b)/p, p^2 = h-(c-1)/24, S11 = S_{1,0} -+ S_{1,1}")
for beta in [mp.mpf('0.3'), mp.mpf('3')]:
    b = mp.sqrt(mp.mpc(0, beta)); b2 = b**2
    S1h = lambda p: 2*mp.sqrt(2)*mp.sinh(2*mp.pi*b*p)*mp.sinh(2*mp.pi*p/b)/p
    p0 = 0.5j*(b+1/b); p1 = 0.5j*(b-1/b)          # p^2 = -(b+-1/b)^2/4 for h=0,1
    small = 4j*b*mp.sin(mp.pi*b2)*mp.sin(mp.pi/b2)/(1-b2**2)
    big = 4j*b**3*mp.sin(mp.pi*b2)*mp.sin(mp.pi/b2)/(b2**2-1)
    m_ = S1h(p0) - S1h(p1); p_ = S1h(p0) + S1h(p1)
    print(f"  beta={mp.nstr(beta,3)}: (S10-S11)/(84,|b|>1) = {mp.nstr(m_/big, 12)};  (S10+S11)/(84,|b|<1) = {mp.nstr(p_/small, 12)}")
    S1p = lambda p: 4*mp.sqrt(2)*mp.sinh(2*mp.pi*b*p)*mp.sinh(2*mp.pi*p/b)
    print(f"             p-basis: S_{{1,0}}-S_{{1,1}} = {mp.nstr(S1p(p0)-S1p(p1), 10)};  -8 sqrt2 sin sin = {mp.nstr(-8*mp.sqrt(2)*mp.sin(mp.pi*b2)*mp.sin(mp.pi/b2), 10)}")
