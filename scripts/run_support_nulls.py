#!/usr/bin/env python3
"""Generate reproducible null binary support matrices for external model evaluation."""
from pathlib import Path
import argparse
import sys
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from microdim.support_controls import permute_species,degree_preserving_rewire

parser=argparse.ArgumentParser()
parser.add_argument('--support-csv',required=True,help='Binary support CSV; first column contains strain labels')
parser.add_argument('--output',type=Path,default=ROOT/'results'/'support_nulls')
parser.add_argument('--replicates',type=int,default=100)
parser.add_argument('--seed',type=int,default=20261005)
parser.add_argument('--rewire-attempts',type=int,default=1000)
a=parser.parse_args()
if a.replicates<1: raise SystemExit('--replicates must be positive')
df=pd.read_csv(a.support_csv,index_col=0)
S=df.to_numpy()
a.output.mkdir(parents=True,exist_ok=True)
rows=[]
for i in range(a.replicates):
    for family in ('species_shuffle','degree_preserving'):
        if family=='species_shuffle': M=permute_species(S,seed=a.seed+i); accepted=0
        else: M,accepted=degree_preserving_rewire(S,a.rewire_attempts,seed=a.seed+i)
        assert M.shape==S.shape
        fname=f'{family}_{i:03d}.csv'
        pd.DataFrame(M,index=df.index,columns=df.columns).to_csv(a.output/fname)
        rows.append({'family':family,'replicate':i,'seed':a.seed+i,'edge_count':int(M.sum()),'swaps_accepted':accepted,'file':fname})
pd.DataFrame(rows).to_csv(a.output/'manifest.csv',index=False)
print(f'Generated {len(rows)} support nulls: {a.output}')
