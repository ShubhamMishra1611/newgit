"""
R1 independent check of ABJM (k=1) Fermi-gas partition functions (Marino-Putrov 1110.4066):
 rho(x1,x2) = 1/(2 pi k) * [2cosh(x1/2)]^{-1/2} [2cosh((x1-x2)/(2k))]^{-1} [2cosh(x2/2)]^{-1/2}
 Z(N) = e_N(eigenvalues of rho)  (grand potential: Xi(z) = det(1 + z rho) = sum_N Z(N) z^N)
Airy: Z_Airy = e^A C^{-1/3} Ai[C^{-1/3}(N-B)], C = 2/(pi^2 k), B = 1/(3k) + k/24,
 A(1) = log(2)/4 - zeta(3)/(8 pi^2)  [constant map at k=1; e.g. Hanada et al. 1202.5300 / Hatsuda-Okuyama].
Discretisation: trapezoid on a uniform grid (exponentially convergent for this analytic, exponentially decaying kernel);
eigenvalues positive => e_N by stable DP (no cancellations). Two grids compared for error control.
"""
import numpy as np
import mpmath as mp

def eigs(k=1, h=0.05, L=90.0):
    x = np.arange(-L, L + h/2, h)
    a = 1/np.sqrt(2*np.cosh(x/2))
    K = 1/(2*np.cosh((x[:, None] - x[None, :])/(2*k)))
    M = (h/(2*np.pi*k))*a[:, None]*K*a[None, :]
    ev = np.linalg.eigvalsh(M)
    return np.sort(ev)[::-1]

def esym(ev, Nmax):
    e = np.zeros(Nmax + 1); e[0] = 1.0
    for lam in ev:
        if lam <= 0: continue
        e[1:] = e[1:] + lam*e[:-1]
    return e

mp.mp.dps = 30
k = 1
C = 2/(mp.pi**2*k); B = mp.mpf(1)/(3*k) + mp.mpf(k)/24
A = mp.log(2)/4 - mp.zeta(3)/(8*mp.pi**2)
def ZAiry(N): return mp.exp(A)*C**(-mp.mpf(1)/3)*mp.airyai(C**(-mp.mpf(1)/3)*(N - B))

Nmax = 8
e1 = esym(eigs(h=0.08, L=80.0), Nmax)
e2 = esym(eigs(h=0.05, L=100.0), Nmax)
print("Z(1) =", e2[1], " (1/4 =", 0.25, ")")
print("Z(2) =", e2[2], " (1/(16 pi) =", 1/(16*np.pi), ")")
print("Z(3) =", e2[3], " ((pi-3)/(64 pi) =", (np.pi - 3)/(64*np.pi), ")")
print("Z(4) =", e2[4], " ((10-pi^2)/(1024 pi^2) =", (10 - np.pi**2)/(1024*np.pi**2), ")")
print(" N   Z(N)                 grid-diff/Z    Z/Z_Airy - 1")
for N in range(1, Nmax + 1):
    print(f" {N}  {e2[N]:.15e}  {abs(e2[N]-e1[N])/e2[N]:.1e}   {e2[N]/float(ZAiry(N)) - 1:.4e}")
