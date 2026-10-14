"""RT3: late-time attractor check for claim 8 (all-positive sets) and claim 9 (R6(+) + O4(-) + F4(+)).
Canonical fields (rho_hat, tau_hat), M_P = 1, e-fold time N:  phi' = psi,  psi' = -(3 - eps)(psi + grad ln V),
eps = |psi|^2/2  (from H^2 = V/(3 - eps)).  Exponent vectors as derived in rt3_slopes.py."""
import numpy as np, math
from scipy.integrate import solve_ivp
s23, s2 = math.sqrt(2/3), math.sqrt(2)
def k(a, b):  # V ~ rho^a tau^b
    return np.array([-a * s23, -b / s2])
K = {'R6': k(-1, -2), 'H3': k(-3, -2), 'NS5': k(-2, -2), 'Q': k(1, -2),
     **{f"F{p}": k(3 - p, -4) for p in range(7)}}
# D-brane exponents: V ~ rho^{(p-6)/2} tau^-3
for p in range(3, 10):
    K[f'D{p}'] = k((p - 6) / 2, -3)
def run(terms, phi0, N1=400.0):
    names = [t[0] for t in terms]; A = np.array([t[1] for t in terms]); KK = np.array([K[n] for n in names])
    def rhs(N, y):
        phi, psi = y[:2], y[2:]
        lw = np.log(np.abs(A)) - KK @ phi
        w = np.sign(A) * np.exp(lw - lw.max())      # rescaled weights (log-sum-exp; V itself underflows)
        V = w.sum()
        g = -(w @ KK) / V
        eps = 0.5 * psi @ psi
        return np.r_[psi, -(3 - eps) * (psi + g)]
    sol = solve_ivp(rhs, (0, N1), np.r_[phi0, 0, 0], rtol=1e-10, atol=1e-12, method="LSODA")
    phi, psi = sol.y[:2, -1], sol.y[2:, -1]
    lw = np.log(np.abs(A)) - KK @ phi
    w = np.sign(A) * np.exp(lw - lw.max())
    return 0.5 * psi @ psi, w.sum(), phi, psi
rng = np.random.default_rng(3)
IIB = ['R6', 'H3', 'F1', 'F3', 'F5', 'D3', 'D5', 'D7', 'D9', 'NS5']
IIA = ['R6', 'H3', 'F0', 'F2', 'F4', 'F6', 'D4', 'D6', 'D8', 'NS5']
for lab, S, target in (("IIB all positive", IIB, 9/7), ("IIA all positive", IIA, 25/19), ("IIB+Q", IIB + ['Q'], 1.0),
                       ("IIB without D9", [n for n in IIB if n != 'D9'], 4/3)):
    eps = []
    for trial in range(4):
        terms = [(n, float(np.exp(rng.normal(0, 1)))) for n in S]
        e, V, phi, psi = run(terms, rng.normal(0, 1, 2))
        eps.append(e)
    print(f"  {lab:18s}: late-time eps = {np.round(eps, 4)}  (gamma_min^2/2 = {target:.4f})")
# mixed (+,-,+) collinear: R6(+) O4(-) F4(+): transverse minimum, flow along +rho_hat with eps = 1/3
# choose coefficients so that V_perp(tau_hat) has a local min with V>0: in x = sqrt2 tau_hat: e^-x - B e^-1.5x + C e^-2x, B=1.9, C=1
terms = [('R6', 1.0), ('D4', -1.9), ('F4', 1.0)]
x0 = 0.4463  # local minimum of V_perp found in rt3 numeric check
e, V, phi, psi = run(terms, np.array([0.0, x0 / math.sqrt(2)]), N1=60.0)
dphi10 = math.sqrt(1.5) * psi[0] - psi[1] / math.sqrt(2)
print(f"  R6(+)+O4(-)+F4(+) from transverse min: eps = {e:.4f} (1/3), V>0: {V>0}, d phi_10/dN = {dphi10:+.3f} (>0: coupling grows)")
