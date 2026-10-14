"""Test the Elvang-Herderschee-Morales nonlinear constraints (arXiv:2601.11705 eq.(1.6)):
   g^2 a20 = (2/5) a00^2,  g^2 a21 = (1/10) a00^2,  g^2 a31 = 2 g^2 a30 - a10 a00
on the planar CHR lambda-family A_lambda (arXiv:2408.03362v1 eqs.(16)-(18),(21)), identified as F(s,u) with g^2 = rho.
rho-independent test: a21/a20 = 1/4.  Second test (rho fixed by first relation): a31 - 2 a30 + a10 a00/rho = 0.
Coefficients from the dual-resonant sum (17): A = -rho/(s t) + sum_{n>=1} Rhat(n,t)/(mu(n)-s),
Rhat(n,t) = (a+lam t)_{n-1}/(b)_{n-1}, a=(1+lam)/2, b=(1+3lam)/2, mu(n)=(n+lam-1)/lam.
a_{k,q} = coefficient of s^{k-q} t^q = sum_n [t^q]Rhat(n,t) / mu(n)^{k-q+1}."""
import mpmath as mp
mp.mp.dps = 25

def coeffs(lam):
    lam = mp.mpf(lam); a = (1+lam)/2; b = (1+3*lam)/2
    R0 = lambda n: mp.gammaprod([a+n-1, b], [a, b+n-1])
    psi1 = lambda n: mp.digamma(a+n-1) - mp.digamma(a)          # sum_{k=0}^{n-2} 1/(a+k)
    psi2 = lambda n: -(mp.polygamma(1, a+n-1) - mp.polygamma(1, a))  # sum 1/(a+k)^2
    mu = lambda n: (n+lam-1)/lam
    # [t^1] Rhat = R0 * lam*psi1 ; [t^2] Rhat = R0*( (lam psi1)^2 - lam^2 psi2 )/2
    T1 = lambda n: R0(n)*lam*psi1(n)
    T2 = lambda n: R0(n)*((lam*psi1(n))**2 - lam**2*psi2(n))/2
    S = lambda f: mp.nsum(f, [1, mp.inf])
    a00 = mp.hyper([1, lam, a], [1+lam, b], 1)                 # eq.(21) closed form
    a00_sum = S(lambda n: R0(n)/mu(n))
    a10 = S(lambda n: R0(n)/mu(n)**2); a11 = S(lambda n: T1(n)/mu(n))
    a20 = S(lambda n: R0(n)/mu(n)**3); a21 = S(lambda n: T1(n)/mu(n)**2); a22 = S(lambda n: T2(n)/mu(n))
    a30 = S(lambda n: R0(n)/mu(n)**4); a31 = S(lambda n: T1(n)/mu(n)**3); a32 = S(lambda n: T2(n)/mu(n)**2)
    return dict(a00=a00, a00_sum=a00_sum, a10=a10, a11=a11, a20=a20, a21=a21, a22=a22, a30=a30, a31=a31, a32=a32)

print("Veneziano check (lambda=1): expect a00=z2, a10=z3, a20=z4, a21=z4/4, a30=z5, a31=2z5-z2z3 (EHM eq.1.8)")
z = mp.zeta
for lam in ['0.25', '0.5', '0.9', '1', '1.1', '2', '4']:
    c = coeffs(lam)
    rho_fix = mp.mpf(2)/5*c['a00']**2/c['a20']       # g^2 from first relation
    test2 = c['a31'] - 2*c['a30'] + c['a10']*c['a00']/rho_fix
    print("lam=%-4s a00=%s (sum %s) a10=%s a11=%s | a20=%s a22=%s a21=%s | a21/a20=%s (need 0.25) | rho*=%s  rel3 residual=%s" % (
        lam, mp.nstr(c['a00'], 9), mp.nstr(c['a00_sum'], 9), mp.nstr(c['a10'], 9), mp.nstr(c['a11'], 9),
        mp.nstr(c['a20'], 9), mp.nstr(c['a22'], 9), mp.nstr(c['a21'], 9), mp.nstr(c['a21']/c['a20'], 9),
        mp.nstr(rho_fix, 9), mp.nstr(test2, 6)))
print("zeta: z2=%s z3=%s z4=%s z4/4=%s z5=%s 2z5-z2z3=%s" % tuple(mp.nstr(x, 9) for x in (z(2), z(3), z(4), z(4)/4, z(5), 2*z(5)-z(2)*z(3))))
