"""R2-T4 script A: the CEM (arXiv:2501.01486) microstate count -- what is input, what is output.

Equation numbers: arXiv-HTML numbering (73)-(78) == published SciPost Phys. 18, 131 numbering (4.32)-(4.37);
(44)-(46) == (4.3)-(4.5); (135) is the eigenvalue density in appendix C.

Checks:
 1. sympy: (T-hat_{1,1})^2 / C_{S^2} is b-independent (= 2/pi^4)  [eq. (75) vs (45), condition (77)].
 2. closed form (74): int_2^{E0} rho_0(E) dE = T-hat/(2 pi i), with rho_0 from (135), E0 = 2cos(pi b^-2);
    numerical quadrature at several b^2 = i*beta; E0 is the first zero of rho_0.
 3. 'input test': feed in ANY sphere partition function Z_in(beta) (correct one-loop, CEM (85), ADLS (8.18),
    and deliberately WRONG ones with a different log coefficient and an arbitrary phase), fix e^{S0} by
    (78) g_s^{-1} = T_{1,1} = 2 pi i |Z_in|^{1/2}, then COUNT eigenvalues numerically and compare
    log N_eff^2 with log|Z_in|.  If the count reproduces every input, (i)-(iv) are inputs, not outputs.
"""
import sympy as sp
import mpmath as mp

mp.mp.dps = 40

# ---------- 1. symbolic identity (77) ----------
b = sp.symbols('b', nonzero=True)
That = 8*b**2*sp.sin(sp.pi*b**2)*sp.sin(sp.pi/b**2)/(1-b**4)            # (75)
CS2 = 32*sp.pi**4*(sp.sin(sp.pi*b**2)*sp.sin(sp.pi/b**2)/(b**2-b**-2))**2  # (45)
ratio = sp.simplify(That**2/CS2)
print("1. (T-hat_11)^2 / C_S2 =", ratio, " (b-independent:", sp.diff(ratio, b) == 0, ")")

# ---------- 2. N_eff closed form (74) vs quadrature of (135) ----------
def rho0(E, beta):
    """(135): rho_0(E) = (2/pi) sinh(-i pi b^2) sin(-i b^2 arccosh(E/2)), with b^2 = i beta."""
    b2 = 1j*beta
    return (2/mp.pi)*mp.sinh(-1j*mp.pi*b2)*mp.sin(-1j*b2*mp.acosh(E/2))

def That_num(beta):
    b2 = mp.mpc(0, beta)
    return 8*b2*mp.sin(mp.pi*b2)*mp.sin(mp.pi/b2)/(1-b2**2)

print("\n2. N_eff/e^{S0}: quadrature of (135) on [2,E0] vs closed form (74) = T-hat/(2 pi i)")
for beta in [0.05, 0.1, 0.25, 0.5, 0.9, 2.0, 5.0]:
    b2 = mp.mpc(0, beta)
    E0 = 2*mp.cos(mp.pi/b2)                       # = 2 cosh(pi/beta), real
    assert abs(mp.im(E0)) < mp.mpf(10)**-30
    E0 = mp.re(E0)
    # quadrature in theta (E = 2 cosh theta) to avoid the sqrt edge
    f = lambda th: rho0(2*mp.cosh(th), beta)*2*mp.sinh(th)
    Iq = mp.quad(f, [0, mp.pi/beta])
    cf = That_num(beta)/(2j*mp.pi)
    # first zero check: rho0 > 0 just below E0 and < 0 just above
    th0 = mp.pi/beta
    s_below = mp.re(rho0(2*mp.cosh(th0*(1-1e-6)), beta)); s_above = mp.re(rho0(2*mp.cosh(th0*(1+1e-6)), beta))
    print(f"  beta={beta:5}: quad={mp.nstr(Iq,15)}  closed={mp.nstr(cf,15)}  rel.diff={mp.nstr(abs(Iq/cf-1),3)}"
          f"  sign(rho) below/above E0: {mp.sign(s_below)},{mp.sign(s_above)}")

# ---------- 3. input test ----------
# Invariant tree-level entropy for the |b|<1 branch: S0 = 2 pi / beta (see script C; S_GH = 2pi(1/beta - beta)).
def logZ_in(name, beta):
    S0 = 2*mp.pi/beta
    if name == 'one-loop (ADLS 1.23, phase -i)':
        return S0 - 3*mp.log(S0) + 5*mp.log(2*mp.pi) - 1j*mp.pi/2
    if name == 'CEM (85)':
        b2 = mp.mpc(0, beta)
        return mp.log(-16*b2*mp.sin(mp.pi*b2)**2*mp.sin(mp.pi/b2)**2/(1-b2**2)**2)
    if name == 'ADLS (8.18)':
        x = S0/(2*mp.pi); K = 2+1j*x
        return S0 + 2*mp.log(abs(mp.sqrt(4/K)*mp.sin(mp.pi/K)))
    if name == 'WRONG: a=-7, b=0, phase 0.3':
        return S0 + 7*mp.log(S0) + 0.3j
    if name == 'WRONG: S0/2 + 11 log S0':
        return S0/2 + 11*mp.log(S0)
    raise ValueError

print("\n3. Count eigenvalues after fixing e^{S0} by (78); compare log N_eff^2 with log|Z_in|")
for name in ['one-loop (ADLS 1.23, phase -i)', 'CEM (85)', 'ADLS (8.18)',
             'WRONG: a=-7, b=0, phase 0.3', 'WRONG: S0/2 + 11 log S0']:
    worst = 0
    for beta in [0.02, 0.05, 0.1, 0.2]:
        lZ = logZ_in(name, beta)
        absZ = mp.exp(mp.re(lZ))
        # (78): T_{1,1} = e^{S0} T-hat = 2 pi i |Z|^{1/2}  ->  e^{S0}
        eS0 = 2j*mp.pi*mp.sqrt(absZ)/That_num(beta)
        assert abs(mp.im(eS0)) < 1e-25*abs(eS0) and mp.re(eS0) > 0
        eS0 = mp.re(eS0)
        f = lambda th: eS0*rho0(2*mp.cosh(th), beta)*2*mp.sinh(th)
        Neff = mp.re(mp.quad(f, [0, mp.pi/beta]))
        dev = abs(2*mp.log(Neff) - mp.re(lZ))
        worst = max(worst, dev)
    print(f"  {name:35s}: max_beta |log N_eff^2 - log|Z_in|| = {mp.nstr(worst, 3)}")

# e^{S0} implied by CEM (85): (87) says e^{S0} = pi/|b|
print("\n   CEM (87) check: e^{S0} implied by (85)+(78) vs pi/|b|:")
for beta in [0.05, 0.3]:
    absZ = mp.exp(mp.re(logZ_in('CEM (85)', beta)))
    eS0 = mp.re(2j*mp.pi*mp.sqrt(absZ)/That_num(beta))
    print(f"   beta={beta}: e^S0={mp.nstr(eS0,15)}, pi/|b|={mp.nstr(mp.pi/mp.sqrt(beta),15)}")
print("\nConclusion: the count returns whatever Z is fed in (including wrong log coefficients and phases):\n"
      "S_GH, -3 log S_GH, 5 log 2pi and the phase are INPUTS of the counting argument; the only b-dependent\n"
      "content is the identity (77), (T-hat)^2/C_S2 = 2/pi^4, internal to the complex Liouville string/matrix model.")
