"""RT2 independent re-derivation of T2 claims 1-5 (data, cs backreaction, BDGL, masses, LR).
Reads ONLY the decoded JSON (never the pickle). Formulas typed from primary sources:
 - 2406.13751 eqs (3.11),(3.32),(B.2),(B.9),(B.16),(B.18) ; Table 2 ; eqs 5.9,5.26,5.39,5.53,5.67 (SUSY pts)
 - 1809.06861 (BDGL) eqs (2.17),(2.32),(2.33),(3.2),(3.3)
 - 2201.03572 (Junghans) eqs (2.1),(2.5),(2.8),(2.9),(2.13),(2.16)
 - authors' code kahler_stabilisation.py (fetched via Exa): V += c_up/Vs^(4/3) + gs^3*CSF/Vs^2, Vs = string-frame vol.
"""
import json, numpy as np
from scipy.integrate import quad
J = json.load(open('/home/user/newgit/string_frontier/phase3/targets/T2/data_extracted_data_dS.json'))
cabs = lambda v: abs(complex(v['re'], v['im'])) if isinstance(v, dict) else abs(v)
N = 5
M = np.array([J['M vector'][i][0] for i in range(N)], float)   # conifold class = first basis element
Kp = np.array([26/5, 93/19, 40/11, 17/5, 29/10])                 # Table 2 of 2406.13751
gs, gsS = np.array(J['gs']), np.array(J['SUSY gs'])
z, zS = np.array([cabs(v) for v in J['zcf']]), np.array([cabs(v) for v in J['SUSY zcf']])
W0, W0S = np.array([cabs(v) for v in J['W0']]), np.array([cabs(v) for v in J['SUSY W0']])
Vt, VtS = np.array(J['volTilde']), np.array(J['SUSY volTilde'])
VE, VEA = np.array(J['corrected CY volume dS']), np.array(J['corrected CY volume AdS'])
VdS, VAdS = np.array(J['V dS']), np.array(J['V AdS'])
cup, CSF = np.array(J['c uplift']), np.array(J['CSFterm'])
H = np.sqrt(VdS/3)
pr = lambda name, arr, f='%11.4g': print('%-46s' % name + ''.join(f % a for a in arr))

print('=== A. data vs paper text (2406.13751 Table 2 and Sec.5) ===')
paper_dS = dict(gs=[0.0657,0.0571,0.0442,0.0404,0.0466], W0=[0.0115,0.00490,0.0222,0.0539,0.0304],
                gsM=[1.051,0.913,0.796,0.808,0.746], z=[2.822e-8,7.934e-9,8.730e-8,1.965e-6,8.703e-7],
                V0=[1.937e-19,1.692e-20,4.983e-19,2.341e-15,2.113e-15])
paper_SUSY = dict(gs=[0.0732,0.0595,0.0450,0.0410,0.0486], z=[1.390e-7,1.591e-8,1.121e-7,2.369e-6,1.468e-6],
                  W0=[0.0103,0.00465,0.0217,0.0525,0.0291], gsM=[1.171,0.953,0.811,0.821,0.778])
pr('M (repo M-vector[0])', M, '%11.0f')
pr('repo gs/paper - 1 (dS)', gs/paper_dS['gs']-1); pr('repo W0/paper - 1 (dS)', W0/paper_dS['W0']-1)
pr('repo |z|/paper - 1 (dS)', z/paper_dS['z']-1); pr('repo V_dS/paper - 1', VdS/paper_dS['V0']-1)
pr('repo gs*M / paper g_sM - 1', gs*M/paper_dS['gsM']-1)
pr('repo SUSY gs/paper -1', gsS/paper_SUSY['gs']-1); pr('repo SUSY |z|/paper -1', zS/paper_SUSY['z']-1)
pr('repo SUSY W0/paper -1', W0S/paper_SUSY['W0']-1); pr('repo SUSY gsM/paper -1', gsS*M/paper_SUSY['gsM']-1)
pr('g_s M^2 (dS)', gs*M**2, '%11.3f'); pr('g_s M^2 (SUSY)', gsS*M**2, '%11.3f')

print('\n=== B. claim 1: F-flat relation  ln(2pi|z|) = -2pi K\'/(g_s M n_cf), n_cf=2 (3.11 with Q_throat~MK\') ===')
res = lambda zz, g: np.log(2*np.pi*zz) + np.pi*Kp/(g*M)
pr('residual SUSY', res(zS, gsS), '%11.4f'); pr('residual dS', res(z, gs), '%11.4f')
pr('change of residual dS-SUSY', res(z, gs)-res(zS, gsS), '%11.4f')
pr('dln g_s', np.log(gs/gsS), '%11.4f'); pr('dln|z|', np.log(z/zS), '%11.4f')
pr('pred. dln|z| = (piK\'/(M))*(1/gsS-1/gs)', np.pi*Kp/M*(1/gsS-1/gs), '%11.4f')
# with paper's (5.39) for Ex3 instead of repo
zS3p, gsS3p = 1.121e-7, 0.0450
print('Ex3 with paper (5.39) SUSY values: residual SUSY = %.4f' % (np.log(2*np.pi*zS3p)+np.pi*Kp[2]/(gsS3p*M[2])))

print('\n=== C. claim 2: uplift bookkeeping (authors code: V_up = c_code/Vs^{4/3}, V_cs = gs^3 CSF/Vs^2, Vs=gs^1.5 VE) ===')
eta = 3*114.037/128  # (3.32) vs (B.10)/(B.11): eta = 3 zeta/128
Sig = 2**(2/3)*quad(lambda x: (x/np.tanh(x)-1)/np.sinh(x)**2*(np.sinh(2*x)-2*x)**(1/3), 0, 50)[0]
zeta_c = 128/Sig*(np.pi/3)**(1/3); eta = 3*zeta_c/128
print('Sigma = %.5f  zeta = %.4f  eta = %.5f' % (Sig, zeta_c, eta))
c_pap = lambda zz, g, vt: eta*zz**(4/3)/(g*M**2*vt**(2/3))
pr('c_code/(g_s^2 c(3.32)) at dS', cup/(gs**2*c_pap(z, gs, Vt)), '%11.5f')
Vs = gs**1.5*VE
Vup = cup/Vs**(4/3); Vcs = gs**3*CSF/Vs**2
pr('V_antiD3 (dS)', Vup); pr('V_cs = CSF/VE^2', Vcs); pr('V_cs/V_antiD3', Vcs/Vup, '%11.3f')
pr('remainder V_F,Kahler = VdS - Vup - Vcs', VdS-Vup-Vcs)
R_cs = c_pap(zS, gsS, VtS)/c_pap(z, gs, Vt)
pr('R_cs = c(SUSY cs)/c(dS cs)', R_cs, '%11.3f')
R3p = (zS3p/zS[2])**(4/3)*(gsS[2]/gsS3p)
print('Ex3 R_cs if paper (5.39) SUSY z,g_s used (Vt from repo): %.3f' % (R_cs[2]*R3p))
pr('Vt*g_s^3 SUSY', VtS*gsS**3, '%11.4f'); pr('Vt*g_s^3 dS', Vt*gs**3, '%11.4f')
pr('|F_cs|/|W| := sqrt(3 Vcs/|V_AdS|)', np.sqrt(3*Vcs/abs(VAdS)), '%11.3f')
pr('W0 VE^{2/3}', W0*VE**(2/3), '%11.2f')

print('\n=== D. claim 3: BDGL/Junghans critical-point condition ===')
I0 = quad(lambda x: (x/np.tanh(x)-1)/np.sinh(x)**2*(np.sinh(2*x)-2*x)**(1/3), 0, 50)[0]
cpp = 2**(1/3)/I0; cp = 1.18
print('I(0)=%.5f c\'\'=2^{1/3}/I(0)=%.5f ; bound 64 pi c\' c\'\'/9 = %.3f ; sqrt-form (8/3)sqrt(pi c\'c\'\')=%.4f' % (I0, cpp, 64*np.pi*cp*cpp/9, 8/3*np.sqrt(np.pi*cp*cpp)))
beta = 4*np.pi*cp*cpp/(gs*M**2)
pr('beta = 4pi c\'c\'\'/(g_sM^2)  (need <= 9/16)', beta, '%11.4f')
pr('bound/g_sM^2', 64*np.pi*cp*cpp/9/(gs*M**2), '%11.3f')
pr('c\' max for a critical point', 9*gs*M**2/(64*np.pi*cpp), '%11.4f')
# direct numerical scan of BDGL (2.17)+(2.32) in strongly-warped regime (3.1): V ~ |S|^{4/3}[ (M L/2pi)^2/c' + pi^{-1} c''/(g_s M^2) * ... ]
# Write V_KS = A |S|^{4/3}/(c' g_s M^2) * g_s * |M L/(2pi)|^2 ,  V_D3 = A pi^{-1} c'' |S|^{4/3}/(g_s M^2), L = ln(S/s_KS)
for i in range(N):
    g, m = gs[i], M[i]
    L = np.linspace(-40, 5, 200001)
    V = np.exp(4*L/3)*(np.pi*g*(m*L/(2*np.pi))**2/(cp*g*m**2) + cpp/(g*m**2))
    dV = np.diff(V); crit = np.sum(np.diff(np.sign(dV)) != 0)
    print('  Ex%d: numerical V(L) on L in[-40,5]: #critical points = %d ; monotonic increasing in L: %s' % (i+1, crit, bool(np.all(dV > 0))))
# check that the normalisation reproduces BDGL (3.2) at large g_sM^2 (sanity)
g, m = 0.1, 40.0  # g_s M^2 = 160
L = np.linspace(-5, 1, 600001); V = np.exp(4*L/3)*(np.pi*g*(m*L/(2*np.pi))**2/(cp*g*m**2) + cpp/(g*m**2))
iL = np.argmin(V[:400000]+0); 
from scipy.signal import argrelextrema
mins = L[argrelextrema(V, np.less)[0]]; maxs = L[argrelextrema(V, np.greater)[0]]
b = 4*np.pi*cp*cpp/(g*m**2)
print('  sanity g_sM^2=160: numeric min L=%s max L=%s ; BDGL(3.2): %.5f, %.5f' % (mins, maxs, -0.75+np.sqrt(9/16-b), -0.75-np.sqrt(9/16-b)))
# strong-warping condition (3.1) via 2406.13751 (B.9): e^{-4A_IR} vs t = VE^{2/3}
gam = 8*np.pi*np.sqrt(3)/Sig**0.75
zhat = z/np.sqrt(8*Vt)
e3A = gam*zhat/(np.sqrt(gs)*M)**1.5
pr('e^{-4A_IR}/t (B.9,B.3)', e3A**(-4/3)/VE**(2/3), '%11.3g')

print('\n=== E. claim 4: masses (Junghans A at c\'\'->0), m_wKK (B.18), Kahler spectrum (repo) ===')
mz2 = np.sqrt(2*np.pi*Vup/(cpp*gs*M**2))/(8*np.pi**2*cp**2)
mz = np.sqrt(mz2)
mwkk = np.sqrt(Sig**-0.5*(3/(8*np.pi))**(1/3)*z**(2/3)/((np.sqrt(gs)*M)**3*VE**(2/3)*Vt**(1/3)))
pr('m_zeta/H', mz/H, '%11.0f'); pr('m_wKK(B.18)/H', mwkk/H, '%11.0f'); pr('m_zeta/m_wKK', mz/mwkk, '%11.3f')
# independent S^3 (l=1) KK estimate from metric (B.3): m^2 = l(l+2) e^{2A_IR}/(4 pi t R^2), R from (B.6),(B.7)
R = (Sig**0.75/(2*np.pi*np.sqrt(6))*(np.sqrt(gs)*M)**1.5/(2*np.pi**2))**(1/3)
mS3 = np.sqrt(3*e3A**(2/3)/(4*np.pi*VE**(2/3)*R**2))
pr('tip S^3 radius R (l_s, Einstein)', R, '%11.3f'); pr('m_KK(S3,l=1)/H', mS3/H, '%11.0f'); pr('m_KK(S3)/m_wKK(B.18)', mS3/mwkk, '%11.2f')
fr, rmin = [], []
for i in range(N):
    km = np.array(J['km mass spectrum in Hubble units'][i])
    fr.append(np.mean(km > mz[i]/H[i])); rmin.append(mz[i]/H[i]/km.min())
pr('fraction Kahler masses > m_zeta', fr, '%11.3f'); pr('m_zeta / lightest Kahler', rmin, '%11.1f')
fr2 = [np.mean(np.array(J['km mass spectrum in Hubble units'][i]) > mwkk[i]/H[i]) for i in range(N)]
pr('fraction Kahler masses > m_wKK(B.18)', fr2, '%11.3f')
# check of Junghans A formula against direct second derivative
for i in [0]:
    zeta_ = 1e-6; Vv = 1e4
    Vz = lambda zz: zz**(4/3)/(2*np.pi**2*cp*Vv**(4/3))*0.25*np.log(zz/zeta_)**2
    h = 1e-9*zeta_
    d2 = (Vz(zeta_+h)-2*Vz(zeta_)+Vz(zeta_-h))/h**2
    Q = cp*1.0*zeta_**(-4/3)/Vv**(2/3)   # g_sM^2 = 1
    print('  check m^2=V\'\'/(2Q) vs zeta^{2/3}/(8pi^2 c\'^2 gsM^2 V^{2/3}): %.6g vs %.6g' % (d2/(2*Q), zeta_**(2/3)/(8*np.pi**2*cp**2*Vv**(2/3))))

print('\n=== F. claim 5: Lust-Randall S/S0 = 0.94 at g_sM^2 ~ 46 (2206.04708 eq 7.5), linear-response extrapolations ===')
x = 46.25/(gs*M**2)
pr('ln-linear: (0.94)^(46.25/gsM^2)', 0.94**x, '%11.3f'); pr('S-linear: 1-0.06*46.25/gsM^2', 1-0.06*x, '%11.3f')
