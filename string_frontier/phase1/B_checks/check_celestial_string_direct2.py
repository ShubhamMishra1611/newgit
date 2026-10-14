# Variant tests for Stieberger-Taylor (3.30)/(3.31): contour above/below the s-channel poles, and both branches of (-r)^beta.
import mpmath as mp
mp.mp.dps=25
def F(r,w):
    s=r*w; u=-w
    return mp.exp(mp.loggamma(1-s)+mp.loggamma(1-u)-mp.loggamma(1-s-u))
def I_direct(r,beta,phi):
    e=mp.exp(1j*phi)
    return mp.mpf(1)/2*mp.quad(lambda t:(t*e)**(-beta-1)*F(r,t*e)*e,[0,0.5,1,2,4,8,16,32,64,mp.inf])
def Ssum(r,beta,N=400):
    return mp.fsum((mp.mpf(n)/r)**(-beta)*mp.gamma(n-n/r)*mp.rgamma(-n/r)/mp.factorial(n) for n in range(1,N))
for r in [mp.mpf(2),mp.mpf('3.7')]:
    for beta in [mp.mpf('-1.3'), mp.mpf('-0.5')+0.2j]:
        up=I_direct(r,beta,mp.mpf('0.2')); dn=I_direct(r,beta,mp.mpf('-0.2')); S=Ssum(r,beta)
        A=mp.pi/2*mp.exp(1j*mp.pi*beta)/mp.sin(mp.pi*beta)*S
        B=mp.pi/2*mp.exp(-1j*mp.pi*beta)/mp.sin(mp.pi*beta)*S
        print(f"r={r} beta={mp.nstr(beta,4)}\n  above={mp.nstr(up,10)} below={mp.nstr(dn,10)} avg={mp.nstr((up+dn)/2,10)} (above-below)={mp.nstr(up-dn,10)}")
        print(f"  (3.31) branch e^+ipb: {mp.nstr(A,10)}   branch e^-ipb: {mp.nstr(B,10)}   -i*pi*S*(n/r)^-b: {mp.nstr(-1j*mp.pi*S,10)}")
