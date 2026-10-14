"""T2 task (2): heavy-conifold self-consistency for the five candidate dS vacua of 2406.13751.
Inputs: data published in 2406.13751 (Table 2, Sec. 5, App. B) and the authors' repository
(AndreasSchachner/kklt_de_sitter_vacua, data/dS_examples/data.p, read code-free by safe_read_pickle.py
-> ../data_extracted_data_dS.json).
Conventions: 2406.13751 units l_s = 2 pi sqrt(alpha') = 1, M_pl = 1; V_E = Einstein-frame CY volume,
Vt = mirror volume ('volTilde'), z = conifold modulus, n_cf = 2 conifolds.
Uplift as implemented in the authors' code (kahler_stabilisation.py): V_up = c_code / Vs^(4/3) with
Vs = g_s^(3/2) V_E the string-frame volume;  c_code = g_s^2 * c_paper with c_paper from (3.32).
CS-sector F-term energy (code): g_s^3 * CSFterm / Vs^2 = CSFterm / V_E^2.
Warped-throat EFT for the conifold modulus: Junghans 2201.03572 eqs.(2.1),(2.2),(2.5),(2.13),(2.16) and
BDGL 1809.06861 eqs.(2.17),(2.24),(2.32),(3.2),(3.3); c' = 1.18 (BDGL 2.18), c'' = 2^(1/3)/I(0).
"""
import json, os
import numpy as np
from scipy.integrate import quad

here = os.path.dirname(os.path.abspath(__file__))
D = json.load(open(os.path.join(here, '..', 'data_extracted_data_dS.json')))
cx = lambda v: complex(v['re'], v['im']) if isinstance(v, dict) else v

I0 = quad(lambda x: (x/np.tanh(x)-1)/np.sinh(x)**2*(np.sinh(2*x)-2*x)**(1/3), 1e-8, 60, limit=400)[0]
Sigma = 2**(2/3)*I0
cpp = 2**(1/3)/I0           # c''
cp = 1.18                   # c'
zeta_B = 128/Sigma*(np.pi/3)**(1/3)   # (B.11) 114.037
eta = 3*zeta_B/128                     # (3.32) 2.6727
gamma = 8*np.pi*np.sqrt(3)/Sigma**0.75  # (B.9) 39.4612
print(f"Sigma={Sigma:.5f} zeta={zeta_B:.3f} eta={eta:.4f} gamma={gamma:.4f} c''={cpp:.4f} c'={cp}")
Kp = [26/5, 93/19, 40/11, 17/5, 29/10]   # K' from Table 2 / eqs. (5.6),(5.23),(5.36),(5.50),(5.64)
QO = [160, 160, 160, 100, 100]          # Q_O = 2 + h11 + h21 (4.5)
ncf = 2

rows = []
for i in range(5):
    M = int(D['M vector'][i][0])  # conifold curve = (1,0,..) basis => M = M^1 (3.7)
    gS, gD = D['SUSY gs'][i], D['gs'][i]
    W0S, W0D = abs(cx(D['SUSY W0'][i])), abs(cx(D['W0'][i]))
    zS, zD = abs(cx(D['SUSY zcf'][i])), abs(cx(D['zcf'][i]))
    VtS, VtD = D['SUSY volTilde'][i], D['volTilde'][i]
    VEA, VED = D['corrected CY volume AdS'][i], D['corrected CY volume dS'][i]
    VAdS, VdS = D['V AdS'][i], D['V dS'][i]
    ccode, CSF = D['c uplift'][i], D['CSFterm'][i]
    km = np.sort(np.abs(np.array(D['km mass spectrum in Hubble units'][i], float)))
    csm = np.sort(np.abs(np.array(D['csm mass spectrum in Hubble units'][i], float)))
    H = np.sqrt(VdS/3)
    r = dict(ID=i+1, M=M, Kp=Kp[i], gS=gS, gD=gD, W0S=W0S, W0D=W0D, zS=zS, zD=zD, VtS=VtS, VtD=VtD,
             VEA=VEA, VED=VED, VAdS=VAdS, VdS=VdS, H=H)
    # --- (a) z_cf tracks its F-flat relation (3.11)/(4.15) with n_cf=2 at both g_s values
    r['resS'] = np.log(2*np.pi*zS) + 2*np.pi*Kp[i]/(gS*M*ncf)
    r['resD'] = np.log(2*np.pi*zD) + 2*np.pi*Kp[i]/(gD*M*ncf)
    # --- (b) uplift energies
    cpaperD = eta*zD**(4/3)/(gD*M**2*VtD**(2/3))
    VsD = gD**1.5*VED
    r['ccheck'] = ccode/(gD**2*cpaperD)                 # should be 1
    r['Vup'] = ccode/VsD**(4/3)
    r['Vcs'] = CSF/VED**2
    r['VkmF'] = VdS - r['Vup'] - r['Vcs']
    # precursor alignment Xi: uplift at SUSY cs values and AdS Kahler point, over |V_AdS|
    cpaperS = eta*zS**(4/3)/(gS*M**2*VtS**(2/3))
    r['Xi_exact'] = cpaperS/VEA**(4/3)/abs(VAdS)
    r['Xi_B16'] = zeta_B*zS**(4/3)*VEA**(2/3)*VtS**(1/3)/(gS*M)**2/W0S**2
    # |V_F| approximation (B.2) with ||Omega||^2 = 8 Vt (B.15)
    r['VF_B2'] = 3/16*gS/VEA**2*W0S**2/(8*VtS)
    # cs relaxation: uplift at SUSY cs point vs dS cs point (same V_E)
    r['R_cs'] = (cpaperS/cpaperD)
    r['dlng'] = np.log(gD/gS)
    r['dlnz'] = np.log(zD/zS)
    r['Vt_g3'] = (VtS*gS**3, VtD*gD**3)
    r['dlnVup_dlng'] = 4*np.pi/3*Kp[i]/(gD*M) + 1   # V_up ∝ z^{4/3} g_s with z ∝ exp(-pi K'/(g_s M)), Vt ∝ g_s^-3
    # --- (c) BDGL / Junghans criterion
    gM2 = gD*M**2
    r['gM2'] = gM2; r['gM2S'] = gS*M**2
    r['beta'] = 4*np.pi*cp*cpp/gM2
    r['disc'] = 9 - 64*np.pi*cp*cpp/gM2      # Junghans (2.8): need >= 0
    r['ratio_bound'] = (64*np.pi*cp*cpp/9)/gM2
    r['cp_crit'] = 9*gM2/(64*np.pi*cpp)      # c' below which a minimum would exist
    # --- (d) warped-EFT conifold mass without anti-D3 (Junghans A at c''=0):
    # m^2 = zeta^{2/3}/(8 pi^2 c'^2 g_s M^2 V^{2/3}); normalisation of zeta fixed by matching the
    # physical anti-D3 energy V_up = c'' zeta^{4/3}/(2 pi g_s M^2 V^{4/3}) (Junghans 2.5) to the repo value.
    zeta_over_V_23 = np.sqrt(2*np.pi*gM2*r['Vup']/cpp)
    m2z = zeta_over_V_23/(8*np.pi**2*cp**2*gM2)
    r['mz_H'] = np.sqrt(m2z)/H
    r['frac_km_heavier'] = np.mean(km > r['mz_H'])
    r['km_min'], r['km_med'], r['km_max'] = km[0], np.median(km), km[-1]
    r['csm_min'], r['csm_max'] = csm[0], csm[-1]
    # Junghans mapping zeta = (eta/(c''/2pi))^{3/4} |z| / sqrt(Vt)  (consistency only)
    r['zeta_J'] = (eta/(cpp/(2*np.pi)))**0.75*zD/np.sqrt(VtD)
    # --- (e) warped KK scale: 2406.13751 (B.18) and an S^3-radius estimate
    m2wKK = Sigma**-0.5*(3/(8*np.pi))**(1/3)*zD**(2/3)/((np.sqrt(gD)*M)**3*VED**(2/3)*VtD**(1/3))
    r['mwKK_H'] = np.sqrt(m2wKK)/H
    # S^3 estimate: lowest S^3 harmonic l=1: k^2 = 3/R^2; R from Vol(A)=2 pi^2 R^3 = (B.7);
    # 4d mass^2 = e^{2A_IR} k^2 l_s^2/(4 pi t) with metric (B.3), t = V_E^{2/3}, e^{3A_IR} from (B.9)
    zhat = zD/np.sqrt(8*VtD)
    e3A = gamma*zhat/(np.sqrt(gD)*M)**1.5
    R = (Sigma**0.75/(2*np.pi*np.sqrt(6))*(np.sqrt(gD)*M)**1.5/(2*np.pi**2))**(1/3)
    m2S3 = e3A**(2/3)*3/R**2/(4*np.pi*VED**(2/3))
    r['mS3_H'] = np.sqrt(m2S3)/H
    r['R_S3_ls'] = R
    # strong-warping ratio e^{-4A_IR}/t  (B.3: e^{-4A} ~ t in the bulk)
    r['warp_ratio'] = e3A**(-4/3)/VED**(2/3)
    # hierarchy per throat (B.17) at Xi=1
    r['rIR_rUV3'] = (3*np.pi)**0.75/2**(15/4)*W0D**1.5*(gD/VtD)**0.75
    r['Qthroat'] = M*Kp[i]; r['QO'] = QO[i]
    r['Kp_eff_per_throat'] = Kp[i]/ncf
    r['warp_e4A'] = np.exp(-8*np.pi*(Kp[i]/ncf)/(3*gD*M))
    # Lust-Randall 2206.04708 eq.(7.5): S/S0 ~ 0.94 at g_s M^2 ~ 46; crude linear-response extrapolation
    r['LR_S_over_S0'] = np.exp(np.log(0.94)*46.0/gM2)
    rows.append(r)

def p(fmt, key, label):
    print(f"{label:46s}" + "".join(fmt.format(r[key]) for r in rows))

print("\n================ DATA (repo + paper) ================")
print(f"{'':46s}" + "".join(f"{'Ex'+str(r['ID']):>12s}" for r in rows))
p("{:12d}", 'M', 'M'); p("{:12.4f}", 'Kp', "K'")
p("{:12.5f}", 'gS', 'g_s (SUSY cs point)'); p("{:12.5f}", 'gD', 'g_s (dS)')
p("{:12.5f}", 'W0S', 'W0 (SUSY)'); p("{:12.5f}", 'W0D', 'W0 (dS)')
p("{:12.3e}", 'zS', '|z_cf| (SUSY)'); p("{:12.3e}", 'zD', '|z_cf| (dS)')
p("{:12.2f}", 'VtS', 'mirror vol Vt (SUSY)'); p("{:12.2f}", 'VtD', 'mirror vol Vt (dS)')
p("{:12.0f}", 'VEA', 'V_E (AdS precursor)'); p("{:12.0f}", 'VED', 'V_E (dS)')
p("{:12.3e}", 'VAdS', 'V_AdS'); p("{:12.3e}", 'VdS', 'V_dS'); p("{:12.3e}", 'H', 'H_dS')
print("\n================ (a) z_cf slaved to its F-flat relation ================")
p("{:12.2e}", 'resS', 'ln(2pi|z|)+pi K\'/(g_s M) at SUSY'); p("{:12.2e}", 'resD', 'same at dS')
print("\n================ (b) uplift bookkeeping ================")
p("{:12.4f}", 'ccheck', 'c_code/(g_s^2 c_paper(3.32))  [=1?]')
p("{:12.3e}", 'Vup', 'V_antiD3 at dS point'); p("{:12.3e}", 'Vcs', 'V_cs F-term energy at dS')
p("{:12.3e}", 'VkmF', 'V_F(Kahler sector) at dS (remainder)')
for r in rows: r['Vcs_over_Vup'] = r['Vcs']/r['Vup']; r['Vup_over_VdS'] = r['Vup']/r['VdS']
p("{:12.3f}", 'Vcs_over_Vup', 'V_cs / V_antiD3'); p("{:12.3f}", 'Vup_over_VdS', 'V_antiD3 / V_dS')
p("{:12.3f}", 'Xi_exact', 'Xi precursor = V_up(SUSY cs,V_E^AdS)/|V_AdS|'); p("{:12.3f}", 'Xi_B16', 'Xi precursor via (B.16)')
for r in rows: r['VFratio'] = r['VF_B2']/abs(r['VAdS'])
p("{:12.3f}", 'VFratio', '(B.2)|V_F| / exact |V_AdS|')
p("{:12.2f}", 'R_cs', 'V_up(SUSY cs)/V_up(dS cs) at fixed V_E')
for r in rows:
    r['Xi_relaxed'] = r['Xi_exact']/r['R_cs']
    r['pos_over_AdS'] = (r['Vup'] + r['Vcs'])/abs(r['VAdS'])
    r['dres'] = r['resD'] - r['resS']
p("{:12.3f}", 'Xi_relaxed', 'Xi_prec / R_cs (uplift after cs relaxation)')
p("{:12.3f}", 'pos_over_AdS', '(V_antiD3 + V_cs)/|V_AdS|')
p("{:12.4f}", 'dres', 'change of F-flat residual SUSY->dS')
p("{:12.4f}", 'dlng', 'ln(g_s^dS/g_s^SUSY)'); p("{:12.3f}", 'dlnz', 'ln(z^dS/z^SUSY)')
p("{:12.2f}", 'dlnVup_dlng', 'dlnV_up/dlng_s (|z| slaved, Vt~g_s^-3)')
print('Vt*g_s^3 (SUSY, dS):', [tuple(round(x, 4) for x in r['Vt_g3']) for r in rows])
print("\n================ (c) BDGL/Junghans conifold criterion at the dS point ================")
p("{:12.2f}", 'gM2', 'g_s M^2 (dS)'); p("{:12.2f}", 'gM2S', 'g_s M^2 (SUSY)')
p("{:12.3f}", 'beta', "beta=4 pi c'c''/(g_sM^2) [min needs<=9/16]")
p("{:12.2f}", 'disc', "9-64pi c'c''/(g_sM^2) (Junghans 2.8)")
p("{:12.2f}", 'ratio_bound', '46.25/(g_s M^2)'); p("{:12.3f}", 'cp_crit', "c' needed for a minimum (c''=1.7546)")
p("{:12.3f}", 'LR_S_over_S0', 'Lust-Randall lin.-response S/S0 (crude)')
print("\n================ (d,e) mass scales in Hubble units (warped EFT, no anti-D3) ================")
p("{:12.0f}", 'mz_H', 'm_conifold/H (Junghans A, c\'\'->0)')
p("{:12.0f}", 'mwKK_H', 'm_wKK/H (2406.13751 B.18)'); p("{:12.0f}", 'mS3_H', 'm_KK,S3/H (S^3 l=1 estimate)')
for r in rows: r['mz_over_mwKK'] = r['mz_H']/r['mwKK_H']; r['mz_over_mS3'] = r['mz_H']/r['mS3_H']
p("{:12.3f}", 'mz_over_mwKK', 'm_conifold/m_wKK(B.18)'); p("{:12.4f}", 'mz_over_mS3', 'm_conifold/m_KK,S3')
p("{:12.2f}", 'km_min', 'Kahler: lightest m/H'); p("{:12.0f}", 'km_med', 'Kahler: median m/H'); p("{:12.0f}", 'km_max', 'Kahler: heaviest m/H')
p("{:12.3f}", 'frac_km_heavier', 'fraction of Kahler modes heavier than m_conifold')
for r in rows: r['mz_over_kmmin'] = r['mz_H']/r['km_min']
p("{:12.0f}", 'mz_over_kmmin', 'm_conifold / lightest Kahler mass')
p("{:12.2f}", 'csm_min', 'bulk cs+dilaton: lightest m/H'); p("{:12.0f}", 'csm_max', 'bulk cs+dilaton: heaviest m/H')
p("{:12.3e}", 'warp_ratio', 'e^{-4A_IR}/t (strong-warping ratio)')
p("{:12.3f}", 'R_S3_ls', 'tip S^3 radius in l_s (Einstein frame)')
p("{:12.3e}", 'rIR_rUV3', '(r_IR/r_UV)^3 at Xi=1 (B.17)')
p("{:12.1f}", 'Qthroat', "Q_throat ~ M K'"); p("{:12d}", 'QO', 'Q_O')
p("{:12.3f}", 'Kp_eff_per_throat', "K'/n_cf (eff. K per throat)"); p("{:12.3e}", 'warp_e4A', "e^{-8pi (K'/2)/(3 g_s M)}")
json.dump([{k: (v if not isinstance(v, tuple) else list(v)) for k, v in r.items()} for r in rows],
          open(os.path.join(here, '..', 't2_heavy_conifold_results.json'), 'w'), indent=1, default=float)

print("\n================ (f) partial F-term uplift: links to 2512.17995 (HLSS) ================")
for r in rows:
    r['eps_over_W0'] = np.sqrt(3*r['Vcs']/abs(r['VAdS']))   # |F_cs|/|W| defined via V_cs/|V_AdS| = |F|^2/(3|W|^2) (OUR definition)
    r['W0V23'] = r['W0D']*r['VED']**(2/3)                   # HLSS eq.(5.67) wants W0 V^{2/3} >> 1 (V: Einstein-frame CY volume here)
    r['gN_over_W0V23'] = r['gD']*r['QO']/r['W0V23']         # HLSS eq.(5.61) parametric ratio g_s N/(eps V^{2/3}) with eps~W0, N~Q_O
p("{:12.2f}", 'eps_over_W0', 'eps/W0 := sqrt(3 V_cs/|V_AdS|)')
p("{:12.1f}", 'W0V23', 'W0 V_E^{2/3}  (HLSS 5.67: >>1 needed)')
p("{:12.2f}", 'gN_over_W0V23', 'g_s Q_O/(W0 V_E^{2/3}) (HLSS 5.61 estimate)')
json.dump([{k: (v if not isinstance(v, tuple) else list(v)) for k, v in r.items()} for r in rows],
          open(os.path.join(here, '..', 't2_heavy_conifold_results.json'), 'w'), indent=1, default=float)
