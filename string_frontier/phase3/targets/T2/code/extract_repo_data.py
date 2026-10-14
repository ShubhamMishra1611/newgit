"""Extract scalar per-vacuum data (and Kahler mass spectra) from the repo pickle, safely.
Run: python3 -I extract_repo_data.py"""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import safe_read_pickle as S
here = os.path.dirname(os.path.abspath(__file__))
src = os.path.join(here, '..', 'sources', 'repo_data_untrusted', sys.argv[1] if len(sys.argv) > 1 else 'data_dS.p')
out = os.path.join(here, '..', 'data_extracted_' + os.path.basename(src).replace('.p', '') + '.json')
keys = ['paper name','paper section','h11','h21','name','SUSY W0','SUSY gs','SUSY gsM','SUSY volTilde','SUSY zcf','gs','W0','gsM','zcf','volTilde','c uplift','CSFterm','V AdS','V dS','corrected CY volume AdS','corrected CY volume dS','lambda','smallest potent curve volume','N_C_min','M vector','K vector','P vector','tau','SUSY tau','km mass spectrum in Hubble units','csm mass spectrum in Hubble units','km mass spectrum in units of 1/lAdS','csm mass spectrum in units of 1/lAdS']
cols = S.main_extract(src, out, keys)
scal = ['paper name','paper section','h11','h21','name','SUSY W0','SUSY gs','SUSY gsM','SUSY volTilde','SUSY zcf',
        'gs','W0','gsM','zcf','volTilde','c uplift','CSFterm','V AdS','V dS','corrected CY volume AdS','corrected CY volume dS',
        'lambda','smallest potent curve volume','N_C_min','M vector','K vector','P vector','tau','SUSY tau']
for k in scal:
    if k in cols:
        v = cols[k]
        print(f'{k:28s}:', [ (complex(x) if isinstance(x,complex) else x) if not isinstance(x, np.ndarray) else ('arr', x.shape) for x in v])
for k in ['km mass spectrum in Hubble units','csm mass spectrum in Hubble units','km mass spectrum in units of 1/lAdS','csm mass spectrum in units of 1/lAdS']:
    if k in cols:
        for r, v in enumerate(cols[k]):
            a = np.sort(np.abs(np.asarray(v, dtype=float)))
            print(f'{k} row {r}: n={a.size} min={a[0]:.4g} 2nd={a[1]:.4g} median={np.median(a):.4g} max={a[-1]:.4g}')
pf = cols.get('Pfaffians')
if pf is not None:
    for r, v in enumerate(pf):
        a = np.asarray(v, dtype=float); print(f'Pfaffians row {r}: n={a.size}, unique nonzero values={np.unique(a[a!=0])}, zeros={int((a==0).sum())}')
dc = cols.get('dual coxeter numbers')
if dc is not None:
    for r, v in enumerate(dc):
        a = np.asarray(v); u, c = np.unique(a, return_counts=True); print(f'dual coxeter row {r}:', dict(zip(u.tolist(), c.tolist())))
