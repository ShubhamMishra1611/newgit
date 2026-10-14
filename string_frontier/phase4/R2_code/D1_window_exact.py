"""R2: exact edges of the KKLT dS window (W0=-1e-4, A=1, a=0.1), uplift D/sigma^n.
Lower edge: local min with V=0 and V'=0  ->  solve {V'(s)=0, V(s)=0} for (s, D).
Upper edge: min and max merge             ->  solve {V'(s)=0, V''(s)=0} for (s, D)."""
import mpmath as mp
mp.mp.dps = 40
W0, A, a = mp.mpf('-1e-4'), mp.mpf(1), mp.mpf('0.1')
def VF(s):
    e = mp.exp(-a*s); W = W0 + A*e
    return a**2*A**2*e**2/(6*s) + a*A*e*W/(2*s**2)
dVF  = lambda s: mp.diff(VF, s)
d2VF = lambda s: mp.diff(VF, s, 2)
for n in [3, 2]:
    # V=0 & V'=0:  D = -VF s^n ;  dVF - n D s^{-n-1} = 0  ->  dVF + n VF/s = 0
    sl = mp.findroot(lambda s: dVF(s) + n*VF(s)/s, 116)
    Dl = -VF(sl)*sl**n
    # V'=0 & V''=0:  D = dVF s^{n+1}/n ;  d2VF + n(n+1) D s^{-n-2} = 0 -> d2VF + (n+1) dVF/s = 0
    su = mp.findroot(lambda s: d2VF(s) + (n+1)*dVF(s)/s, 125)
    Du = dVF(su)*su**(n+1)/n
    print(f"n={n}: lower edge D={mp.nstr(Dl,6)} at sigma={mp.nstr(sl,6)};  upper edge D={mp.nstr(Du,6)} at sigma={mp.nstr(su,6)}")
