# DSSYK chord transfer matrix (Berkooz et al 1811.02584; Verlinde 2402.00635 eq.(2)):
# H|n> = |n+1> + [n]_q |n-1>, symmetrized: off-diag sqrt([n]_q). Check spectrum edges 2/sqrt(1-q)
# and density vs q-Hermite weight  rho(theta) ∝ (q;q)_inf |(e^{2i theta};q)_inf|^2 , E = 2cos(theta)/sqrt(1-q).
import numpy as np
def qpoch(a,q,K=4000):
    out=np.ones_like(a,dtype=complex)
    for k in range(K): out*= (1-a*q**k)
    return out
for lam in [0.3,1.0]:
    q=np.exp(-lam); Nmax=1500
    n=np.arange(1,Nmax)
    off=np.sqrt((1-q**n)/(1-q))
    T=np.diag(off,1)+np.diag(off,-1)
    w,v=np.linalg.eigh(T)
    Emax=2/np.sqrt(1-q)
    # spectral measure of |0> : weights v[0,:]^2  -> moments; compare moments with q-Hermite weight
    wt=v[0,:]**2
    th=np.linspace(1e-6,np.pi-1e-6,20001)
    rho=np.real(qpoch(np.array([q+0j]),q)[0]*qpoch(np.exp(2j*th),q)*qpoch(np.exp(-2j*th),q))/(2*np.pi)
    E=2*np.cos(th)/np.sqrt(1-q)
    norm=np.trapezoid(rho,th)
    mom_exact=[np.trapezoid(rho*E**k,th)/norm for k in (2,4,6)]
    mom_mat=[np.sum(wt*w**k) for k in (2,4,6)]
    # chord-diagram moments: <H^2>=1, <H^4>=2+q, <H^6>=5+6q+3q^2+q^3 (Touchard-Riordan)
    mom_chord=[1,2+q,5+6*q+3*q**2+q**3]
    print(f"lam={lam}: max eig={w.max():.6f}, 2/sqrt(1-q)={Emax:.6f}, weight-normalization={norm:.6f}")
    print("  moments matrix  :",np.round(mom_mat,8))
    print("  moments q-Herm  :",np.round(mom_exact,8))
    print("  moments chords  :",np.round(mom_chord,8))
