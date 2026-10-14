"""Check GPV (arXiv:2102.02847) eqs.(4),(5),(30): alpha_min^IIB = E_{3/2}(e^{i pi/3})/64 ~ 0.1389, alpha_min^IIA ~ 0.1403.
Also brute-force lattice sum of E_s at tau = e^{i pi/3} and tau = i as a cross-check, and check where E_{3/2}, E_{5/2} are minimal
along the boundary of the fundamental domain (sample)."""
import mpmath as mp
mp.mp.dps = 20
E32_rho = mp.mpf(3)**0.25*mp.zeta(1.5)*(mp.zeta(1.5, mp.mpf(1)/3) - mp.zeta(1.5, mp.mpf(2)/3))/mp.sqrt(2)
print("eq.(30): E_3/2(rho) =", E32_rho, " /64 =", E32_rho/64)
z3 = mp.zeta(3)
print("eq.(5) IIA min:", mp.pi**1.5*z3**0.25/(24*mp.sqrt(3)), " at gs^2=9 z3/pi^2 ->",
      (lambda g: z3/(32*g**1.5) + g**0.5*mp.pi**2/96)(mp.sqrt(9*z3/mp.pi**2)))
# Fourier expansion of E_s (standard): E_s = 2 zeta(2s) y^s + 2 sqrt(pi) Gamma(s-1/2)/Gamma(s) zeta(2s-1) y^{1-s}
#   + (8 pi^s sqrt(y)/Gamma(s)) sum_{k>=1} k^{s-1/2} sigma_{1-2s}(k) K_{s-1/2}(2 pi k y) cos(2 pi k x)
def E(s, x, y, K=60):
    s = mp.mpf(s)
    val = 2*mp.zeta(2*s)*y**s + 2*mp.sqrt(mp.pi)*mp.gamma(s-0.5)/mp.gamma(s)*mp.zeta(2*s-1)*y**(1-s)
    tot = 0
    for k in range(1, K):
        sig = sum(mp.mpf(dd)**(1-2*s) for dd in range(1, k+1) if k % dd == 0)
        tot += k**(s-0.5)*sig*mp.besselk(s-0.5, 2*mp.pi*k*y)*mp.cos(2*mp.pi*k*x)
    return val + 8*mp.pi**s*mp.sqrt(y)/mp.gamma(s)*tot
print("Fourier E_3/2(rho) =", E(1.5, 0.5, mp.sqrt(3)/2), "  E_3/2(i) =", E(1.5, 0, 1))
# brute-force lattice sum check at tau=i (slowly convergent; just sanity)
def E_lat(s, x, y, N=400):
    tot = mp.mpf(0)
    for m in range(-N, N+1):
        for n in range(-N, N+1):
            if m == 0 and n == 0: continue
            tot += y**s/((m*x+n)**2 + (m*y)**2)**s
    return tot
print("lattice E_5/2(i), N=150:", E_lat(2.5, 0, 1, 150), " Fourier:", E(2.5, 0, 1))
# boundary scan |tau|=1, x in [0,1/2] and line x=1/2, y >= sqrt3/2
for s in [1.5, 2.5]:
    pts = [(mp.cos(th), mp.sin(th)) for th in mp.linspace(mp.pi/3, mp.pi/2, 7)] + [(mp.mpf('0.5'), y) for y in [0.9, 1.0, 1.2, 1.5]]
    vals = [(mp.nstr(x, 4), mp.nstr(y, 4), mp.nstr(E(s, x, y), 10)) for x, y in pts]
    print("E_%s on boundary:" % s, vals)
