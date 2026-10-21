#!/usr/bin/env python3
"""Compute hCom2* consumption-matrix spectral and support descriptors.

Pass an exported consumer-by-resource matrix from the pinned public source.
The community prediction and resource fitting use that source model separately.
"""
exec(open(__file__.replace('run_hcom2_analysis.py','_bootstrap.py')).read())
import argparse,pandas as pd,numpy as np
from microdim.metrics import participation_dimension,entropy_dimension,energy_rank,leading_energy
from microdim.support_strength import ordinary_nmf,support_latent
ap=argparse.ArgumentParser(); ap.add_argument('--consumption-csv',required=True); a=ap.parse_args()
if not a.consumption_csv:
 raise SystemExit('A 63-by-98 consumption CSV is required. Obtain the matrix from the pinned model source; see docs/HCOM2.md')
C=pd.read_csv(a.consumption_csv,index_col=0).to_numpy(float)
if C.ndim!=2 or not np.isfinite(C).all() or C.min()<0:
    raise SystemExit('Expected a finite, nonnegative consumption matrix.')
rows=[]
for name,X in [('original',C),('binary_support',(C>0).astype(float)),('row_normalized',C/np.maximum(np.linalg.norm(C,axis=1,keepdims=True),1e-15))]:
 rows.append([name,participation_dimension(X),entropy_dimension(X),energy_rank(X,.9),energy_rank(X,.95),leading_energy(X),np.linalg.matrix_rank(X)])
(ROOT/'results').mkdir(exist_ok=True)
pd.DataFrame(rows,columns=['matrix','participation_dimension','entropy_dimension','d90','d95','leading_energy','numerical_rank']).to_csv(ROOT/'results/hcom2_geometry.csv',index=False)
print(pd.DataFrame(rows,columns=['matrix','participation_dimension','entropy_dimension','d90','d95','leading_energy','numerical_rank']).to_string(index=False))
