"""
R1 independent derivation of the mapping radius of the equal-lambda hyperbolic three-string vertex
(Firat 2102.03936), NOT using his closed formula (5) for v.

Fuchsian equation psi'' + (1/2) T psi = 0, T = d/z^2 + d/(z-1)^2 - d/(z(z-1)), d = (1+lambda^2)/2
(exponents (1 +- i lambda)/2 at 0,1,inf: hyperbolic monodromy -e^{-+pi lambda}).
psi_pm = z^{(1 +- i l)/2}(1+...) (Frobenius at 0);  phi_pm(z) = psi_pm(1-z) (symmetry z<->1-z for equal lambdas).
Connection psi_+ = alpha phi_+ + beta phi_-.  Requiring real SL(2,R) monodromy around z=1 in the basis
e^{+-iv/2} psi_pm (Firat eq. (24)) gives  v = -arg(alpha beta)  (mod pi).
Scaled ratio rho_1 = e^{v/l}(psi_+/psi_-)^{1/(i l)};  on (0,1): Phi(z) = l log rho_1 = v + 2 arg psi_+(z) (continuous).
Metric e^phi = |d(l log rho)|^2 / sin^2(l log|rho|) is singular on l log|rho| in pi Z; the geodesic boundary of hole 1 is
the waist l log|rho_1| = pi(l~ + 1/2) of the annulus that contains the seam point z = 1/2 (no singular curve between
the waist and the pants).  Hence l~ = floor(Phi(1/2)/pi) and  l log r = pi(l~ + 1/2) - v ,  r = mapping radius
(z = r w + ..., |w|=1 on the geodesic).  Tachyon cubic coefficient (alpha'=2): (1/3) r^{-6}  (Firat eq. (7), first equality).
"""
import mpmath as mp
mp.mp.dps = 40

def frob(l, s, M=400):
    d = (1 + l*l)/2
    # T(z) z^2 = d + d z^2/(z-1)^2 - d z/(z-1) = sum_n t_n z^n
    t = [mp.mpc(0)]*M
    t[0] += d
    for n in range(M):              # z^2/(1-z)^2 = sum_{n>=2} (n-1) z^n
        if n >= 2: t[n] += d*(n - 1)
        if n >= 1: t[n] += d       # -d z/(z-1) = d z/(1-z) = d sum_{n>=1} z^n
    a = [mp.mpc(1)] + [mp.mpc(0)]*(M - 1)
    # (n+s)(n+s-1) a_n + (1/2) sum_k t_k a_{n-k} = 0
    for n in range(1, M):
        acc = sum(t[k]*a[n - k] for k in range(1, n + 1))/2
        a[n] = -acc/((n + s)*(n + s - 1) + t[0]/2)
    return a

def ev(a, s, z):
    S = mp.mpc(0); dS = mp.mpc(0)
    for n in range(len(a) - 1, -1, -1):
        S = S*z + a[n]
    for n in range(len(a) - 1, 0, -1):
        dS = dS*z + n*a[n]
    zs = mp.exp(s*mp.log(z))
    return zs*S, zs*(s*S/z + dS)

def mapping_radius(l, verbose=False):
    sp_, sm_ = (1 + 1j*l)/2, (1 - 1j*l)/2
    ap, am = frob(l, sp_), frob(l, sm_)
    z0 = mp.mpf('0.5')
    p, dp = ev(ap, sp_, z0)
    fp, dfp = ev(ap, sp_, 1 - z0); dfp = -dfp
    fm, dfm = ev(am, sm_, 1 - z0); dfm = -dfm
    # psi_+ = alpha phi_+ + beta phi_-
    det = fp*dfm - fm*dfp
    alpha = (p*dfm - fm*dp)/det
    beta = (fp*dp - p*dfp)/det
    v = -mp.arg(alpha*beta)
    # continuous arg of psi_+ on (0, 1/2]
    K = 60 + int(40*l)     # fine grid: arg of the regular factor must change < pi between samples
    zs = [z0*mp.exp(-mp.mpf(15)*k/K) for k in range(K, -1, -1)]
    zs = sorted(set([z for z in zs if z <= z0]))
    S0 = None; argacc = None
    for z in zs:
        S = sum(ap[n]*z**n for n in range(len(ap)))
        if argacc is None:
            argacc = mp.arg(S); last = mp.arg(S)
        else:
            a_ = mp.arg(S); dlt = a_ - last
            dlt -= 2*mp.pi*mp.nint(dlt/(2*mp.pi)); argacc += dlt; last = a_
    Phi_half = v + 2*((l/2)*mp.log(z0) + argacc)
    lt = mp.floor(Phi_half/mp.pi)
    loglr = mp.pi*(lt + mp.mpf(1)/2) - v
    r = mp.exp(loglr/l)
    if verbose:
        print(f"   v(mod pi) = {mp.nstr(v,12)}, Phi(1/2)/pi = {mp.nstr(Phi_half/mp.pi,8)}, l~ = {int(lt)}, r = {mp.nstr(r,10)}")
    return r, v, lt

lam = mp.asinh(1)/mp.pi
print("lambda* = arcsinh(1)/pi =", mp.nstr(lam, 12))
r, v, lt = mapping_radius(lam, verbose=True)
print("   cubic coefficient (1/3) r^-6 =", mp.nstr(r**-6/3, 10), "   [Firat eq.(8): 1.62187e8]")
# compare with Firat's printed formula (5) evaluated mod pi
def gam(x): return mp.gamma(1 - x)/mp.gamma(x)
i = 1j
X = mp.gamma(-i*lam)**2/mp.gamma(i*lam)**2*gam((1 + 3*i*lam)/2)/gam((1 - i*lam)/2)
v5 = mp.log(X)/(2*i)
Xp = mp.gamma(-i*lam)**2/mp.gamma(i*lam)**2*gam((1 - i*lam)/2)/gam((1 + 3*i*lam)/2)
v5p = mp.log(Xp)/(2*i)
print(f"   printed (5): v mod pi = {mp.nstr(mp.re(v5) % mp.pi,12)};  agent's variant (gamma ratio inverted): {mp.nstr(mp.re(v5p) % mp.pi,12)};  "
      f"monodromy-derived v mod pi = {mp.nstr(v % mp.pi,12)};  -v mod pi = {mp.nstr((-v) % mp.pi, 12)}")
print("\nminimal-area limit (6561/4096 = %s):" % mp.nstr(mp.mpf(6561)/4096, 10))
for l in [1, 2, 5, 10, 30, 100]:
    r, v, lt = mapping_radius(mp.mpf(l))
    print(f"   lambda={l:4d}: l~={int(lt):4d}  (1/3) r^-6 = {mp.nstr(r**-6/3, 10)}")
