"""
R2 independent check of E3 (Guerrieri-Penedones-Vieira 2102.02847 eqs. 4,5):
  alpha_IIB = E_{3/2}(tau)/2^6,  alpha_IIA = zeta(3)/(32 g^{3/2}) + pi^2 g^{1/2}/96 >= pi^{3/2} zeta(3)^{1/4}/(24 sqrt 3).
Normalisation fixed by matching eq.5 to the perturbative part of E_{3/2}:
  E_s(tau) = sum'_{(m,n)} tau2^s/|m+n tau|^{2s} = 2 zeta(2s) tau2^s + 2 sqrt(pi) Gamma(s-1/2) zeta(2s-1)/Gamma(s) tau2^{1-s} + ...
  (s=3/2: 2 zeta3 tau2^{3/2} + 4 zeta2 tau2^{-1/2}; with tau2=1/g_s and /64 this is exactly eq. 5).
At tau = rho = e^{i pi/3}: |m+n rho|^2 = m^2+mn+n^2 and sum' (m^2+mn+n^2)^{-s} = 6 zeta(s) L(s,chi_{-3}).
Three independent evaluations: (i) zeta*L closed form, (ii) Fourier (Bessel) expansion, (iii) direct lattice sum
with Ewald-free brute force + tail estimate.  Plus a 2D scan of the whole fundamental domain for the minimum.
"""
import mpmath as mp
mp.mp.dps = 30

def E_closed_rho(s):
    s = mp.mpf(s)
    L = mp.dirichlet(s, [0, 1, -1])
    return (mp.sqrt(3)/2)**s * 6*mp.zeta(s)*L

def E_fourier(s, tau, N=40):
    s = mp.mpf(s); t1, t2 = mp.re(tau), mp.im(tau)
    val = 2*mp.zeta(2*s)*t2**s + 2*mp.sqrt(mp.pi)*mp.gamma(s-mp.mpf(1)/2)*mp.zeta(2*s-1)/mp.gamma(s)*t2**(1-s)
    pref = 8*mp.pi**s/mp.gamma(s)*mp.sqrt(t2)
    for n in range(1, N+1):
        sig = mp.fsum(mp.mpf(d)**(1-2*s) for d in range(1, n+1) if n % d == 0)
        val += pref*mp.mpf(n)**(s-mp.mpf(1)/2)*sig*mp.besselk(s-mp.mpf(1)/2, 2*mp.pi*n*t2)*mp.cos(2*mp.pi*n*t1)
    return val

def E_brute(s, tau, R=400):
    # direct sum over |m|,|n| <= R  (float), plus continuum tail estimate  2 pi tau2^{s-1}... (rough)
    import numpy as np
    t1, t2 = float(mp.re(tau)), float(mp.im(tau))
    m = np.arange(-R, R+1)
    M, Nn = np.meshgrid(m, m, indexing='ij')
    q = ((M + Nn*t1)**2 + (Nn*t2)**2)/t2      # unit-determinant quadratic form
    q[R, R] = np.inf
    mask = (np.abs(M) <= R) & (np.abs(Nn) <= R)
    val = np.sum(q[mask]**(-float(s)))
    return val

if __name__ == "__main__":
    rho = mp.mpc(0.5, mp.sqrt(3)/2)
    for s in ['1.5', '2.5']:
        c = E_closed_rho(s); f = E_fourier(s, rho)
        print(f"E_{s}(rho): closed-form 6 zeta L = {mp.nstr(c,12)}, Fourier = {mp.nstr(f,12)}, brute(R=400, no tail) = {E_brute(s, rho):.6f}")
    print("alpha_IIB min = E_3/2(rho)/64 =", mp.nstr(E_closed_rho('1.5')/64, 10))
    aIIA = mp.pi**1.5*mp.zeta(3)**0.25/(24*mp.sqrt(3))
    g = mp.findroot(lambda g: mp.diff(lambda x: mp.zeta(3)/(32*x**1.5)+mp.pi**2*x**0.5/96, g), 1.2)
    print("alpha_IIA min closed form =", mp.nstr(aIIA,10), "; numeric min at g_s =", mp.nstr(g,8),
          "value", mp.nstr(mp.zeta(3)/(32*g**1.5)+mp.pi**2*g**0.5/96, 10), "; 3 sqrt(zeta3)/pi =", mp.nstr(3*mp.sqrt(mp.zeta(3))/mp.pi,8))
    # scan fundamental domain |tau|>=1, |tau1|<=1/2 for E_{3/2} and E_{5/2}
    for s in ['1.5', '2.5']:
        best = None
        for i in range(0, 26):
            t1 = mp.mpf(i)/50                  # 0..0.5 (E even in tau1)
            t2min = mp.sqrt(1 - t1**2)
            for j in range(0, 30):
                t2 = t2min + mp.mpf(j)**2/200
                v = E_fourier(s, mp.mpc(t1, t2), N=25)
                if best is None or v < best[0]:
                    best = (v, t1, t2)
        print(f"scan E_{s}: min {mp.nstr(best[0],10)} at tau=({mp.nstr(best[1],4)}, {mp.nstr(best[2],6)}); rho=(0.5, 0.866025)")
    # weak-coupling check of the D^4R^4 normalisation: tree-level VS coefficient of (s^2+t^2+u^2) is zeta(5);
    # E_{5/2} = 2 zeta(5) tau2^{5/2} + ...  so the matching constant c in  c*E_{5/2} -> zeta(5) g^{-5/2} is 1/2.
    print("E_5/2 leading coefficient 2 zeta(5) =", mp.nstr(2*mp.zeta(5),10), "; E_3/2 leading 2 zeta(3) =", mp.nstr(2*mp.zeta(3),10))
