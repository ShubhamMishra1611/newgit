# Direct numerical test of Stieberger-Taylor 1806.05688 eq.(3.31) (Hankel/residue representation) for Re(beta)<0,
# where the energy integral (3.15) I = 1/2 int_0^inf w^{-beta-1} F_I(r w, -w) dw converges at w->0 (no delta term).
# s-channel poles at s=rw=n are handled with +i eps by rotating the w-contour slightly into the upper half plane.
import mpmath as mp
mp.mp.dps=25
def F(r,w):
    s=r*w; u=-w
    return mp.exp(mp.loggamma(1-s)+mp.loggamma(1-u)-mp.loggamma(1-s-u))
def I_direct(r,beta,phi=mp.mpf('0.15')):
    e=mp.expjpi(phi/mp.pi) if False else mp.exp(1j*phi)
    f=lambda t: (t*e)**(-beta-1)*F(r,t*e)*e
    return mp.mpf(1)/2*mp.quad(f,[0,0.5,1,2,4,8,16,32,64,mp.inf])
def I_res(r,beta,N=400):
    S=mp.fsum((mp.mpf(n)/r)**(-beta)*mp.gamma(n-n/r)*mp.rgamma(-n/r)/mp.factorial(n) for n in range(1,N))
    return 1j*mp.pi/(1-mp.exp(-2j*mp.pi*beta))*S
for r in [mp.mpf(2),mp.mpf('3.7')]:
    for beta in [mp.mpf('-0.5')+0.2j, mp.mpf('-1.3'), mp.mpf('-0.25')-0.4j]:
        d=I_direct(r,beta); d2=I_direct(r,beta,mp.mpf('0.3')); rr=I_res(r,beta)
        print(f"r={r} beta={mp.nstr(beta,4)}: direct={mp.nstr(d,10)}  direct(phi=.3)={mp.nstr(d2,10)}  residue-formula={mp.nstr(rr,10)}")
