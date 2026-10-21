#!/usr/bin/env python3
"""Fit ordinary NMF and support-preserving latent strength approximations.

This script produces reconstruction-selected models and diagnostics. Biological
functional sufficiency additionally requires per-candidate resource re-fitting
and forward prediction with the source-study ecological kernel.
"""
from pathlib import Path
import sys,argparse
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from microdim.support_strength import ordinary_nmf,support_latent
from microdim.metrics import participation_dimension

def parse_ints(s):
    try: x=[int(v.strip()) for v in s.split(',')]
    except ValueError as e: raise argparse.ArgumentTypeError('Expected comma-separated integers') from e
    if not x or min(x)<0: raise argparse.ArgumentTypeError('Use nonnegative integer dimensions')
    return sorted(set(x))

p=argparse.ArgumentParser()
p.add_argument('--consumption-csv',required=True,help='Numeric consumption CSV with strain names as first column')
p.add_argument('--ordinary-ranks',type=parse_ints,default=[1,2,5,10,20,40,50,51,52,53,54,55,56,57,58,59,60])
p.add_argument('--support-dimensions',type=parse_ints,default=[0,1,2,3,4,5,6,8,10,12,15,20])
p.add_argument('--seeds',type=parse_ints,default=list(range(20)))
p.add_argument('--output',type=Path,default=ROOT/'results'/'compression')
a=p.parse_args()
frame=pd.read_csv(a.consumption_csv,index_col=0)
C=frame.to_numpy(dtype=float)
if C.ndim!=2 or C.size==0 or not np.isfinite(C).all() or np.any(C<0):
    raise SystemExit('Consumption input must be a nonnegative finite matrix')
if np.any(C>1):
    raise SystemExit('Support-conditioned sigmoid strengths require input consumption values in [0,1]')
a.output.mkdir(parents=True,exist_ok=True)
rows=[]
for family,dimensions in [('ordinary',a.ordinary_ranks),('support',a.support_dimensions)]:
    for k in dimensions:
        if family=='ordinary' and (k<1 or k>min(C.shape)):
            raise SystemExit(f'Ordinary NMF rank {k} outside [1,min(shape)]')
        candidates=[]
        for seed in a.seeds:
            if family=='ordinary':
                fitted,error=ordinary_nmf(C,k,seed)
                loss=float(error**2)
                success=True
            else:
                fitted,error,success=support_latent(C,k,seed)
                loss=float(error)
            candidates.append((loss,seed,fitted,success))
            rows.append({'family':family,'dimension':k,'seed':seed,'reconstruction_loss':loss,'converged':success,'participation_dimension':participation_dimension(fitted)})
        loss,seed,fitted,ok=min(candidates,key=lambda item:(item[0],item[1]))
        pd.DataFrame(fitted,index=frame.index,columns=frame.columns).to_csv(a.output/f'{family}_k{k}_selected.csv')
        print(f'{family}: dimension={k:2d} seed={seed:2d} reconstruction_loss={loss:.6g} converged={ok}')
pd.DataFrame(rows).to_csv(a.output/'reconstruction_grid.csv',index=False)
print('Biological prediction is not evaluated by the reconstruction objective. See docs/HCOM2.md.')
