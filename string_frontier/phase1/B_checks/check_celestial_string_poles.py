# Stieberger-Taylor 1806.05688 eq.(3.30)/(3.31): non-delta part of celestial type-I 4-gluon energy integral
#   I(r,beta) - 2pi delta = i*pi/(1-exp(-2 pi i beta)) * sum_n (n/r)^(-beta) Res_{s=n} B(-s,1+s/r),
#   Res_{s=n} B(-s,1+s/r) = Gamma(n - n/r)/(n! Gamma(-n/r)).
# Question: at positive integer beta=k, 1/(1-e^{-2 pi i beta}) has a pole; does the residue sum vanish there?
# (APRS 2012.04208 argue no poles at Re beta>0 for exponentially soft UV.)  Quick numerical look only.
import mpmath as mp
mp.mp.dps=30
def S(r,beta,N=600):
    return mp.fsum((mp.mpf(n)/r)**(-beta)*mp.gamma(n-n/r)*mp.rgamma(-n/r)/mp.factorial(n) for n in range(1,N))
for r in [mp.mpf(2),mp.mpf('3.7')]:
    print('convergence check N=300 vs 600, beta=1:', mp.nstr(S(r,1,300),15), mp.nstr(S(r,1,600),15))
    for k in [1,2,3]:
        print(f"r={r}, beta={k}:  residue-sum S = {mp.nstr(S(r,k),12)}")
    for b in [-2,-4]:
        print(f"r={r}, beta={b}:  S = {mp.nstr(S(r,b),12)}")
