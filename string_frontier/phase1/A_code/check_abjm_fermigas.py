"""ABJM U(N)_k x U(N)_{-k} S^3 partition function from the Fermi-gas density matrix (Marino-Putrov 1110.4066),
compared with the all-orders-in-1/N Airy function (Fuji-Hirano-Moriyama 1106.4631, Marino-Putrov 1110.4066).
rho(x,y) = 1/(2 pi k) * (2cosh(x/2))^{-1/2} * (2cosh((x-y)/(2k)))^{-1} * (2cosh(y/2))^{-1/2}
Z(N) from Fredholm expansion of det(1 + z rho). Airy: Z ~ e^{A} C^{-1/3} Ai[C^{-1/3}(N-B)], C = 2/(pi^2 k), B = k/24 + 1/(3k),
A(k) = 2 zeta(3)/(pi^2 k) (1 - k^3/16) + k^2/pi^2 int_0^inf dx x/(e^{kx}-1) log(1-e^{-2x})   [constants: eq. numbers not checked]
"""
import numpy as np, mpmath as mp
from numpy.polynomial.legendre import leggauss
k = 1
L, M = 80.0, 1600                       # x in [-L, L], Gauss-Legendre nodes
t, w = leggauss(M); x = L*t; w = L*w
def ch(u): return 2*np.cosh(u)
R = (1/(2*np.pi*k)) * np.outer(ch(x/2)**-0.5, ch(x/2)**-0.5) / ch((x[:, None]-x[None, :])/(2*k))
K = np.sqrt(w)[:, None]*R*np.sqrt(w)[None, :]          # symmetric Nystrom matrix
lam = np.linalg.eigvalsh(K)
Nmax = 8
# Z(N) = e_N(lam) elementary symmetric polynomials of the eigenvalues
e = np.zeros(Nmax+1); e[0] = 1.0
for l in lam:
    e[1:] = e[1:] + l*e[:-1]
C = 2/(np.pi**2*k); B = k/24 + 1/(3*k)
A = 2*mp.zeta(3)/(mp.pi**2*k)*(1 - k**3/16.) + k**2/mp.pi**2*mp.quad(lambda u: u/(mp.exp(k*u)-1)*mp.log(1-mp.exp(-2*u)), [0, mp.inf])
print("A(1) =", A, "  compare log(2)/4 - zeta(3)/(8 pi^2) =", mp.log(2)/4 - mp.zeta(3)/(8*mp.pi**2))
for N in range(1, Nmax+1):
    Zair = float(mp.exp(A)*C**(-1/3.)*mp.airyai(C**(-1/3.)*(N-B)))
    print(f"N={N}: Z_fermigas = {e[N]:.12e}   Z_Airy = {Zair:.12e}   ratio-1 = {e[N]/Zair-1:+.3e}")
print("exact small-N values (literature, not verified here): Z(1,k=1)=1/4 ->", 0.25, " Z(2,k=1)=1/(16 pi) ->", 1/(16*np.pi))
