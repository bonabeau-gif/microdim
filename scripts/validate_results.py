#!/usr/bin/env python3
"""Report discrepancies between rerun result tables and manuscript reference tables.

This is a comparison, not a proof of exact computational equivalence: optimizer
randomness and preprocessing details can produce expected variation.
"""
from pathlib import Path
import pandas as pd
import numpy as np
root=Path(__file__).resolve().parents[1]
for p in sorted((root/'results').glob('*.csv')) if (root/'results').exists() else []:
    ref=root/'expected'/p.name
    if not ref.exists(): continue
    x,y=pd.read_csv(ref),pd.read_csv(p)
    nc=[c for c in x.columns.intersection(y.columns) if pd.api.types.is_numeric_dtype(x[c]) and pd.api.types.is_numeric_dtype(y[c])]
    if x.shape[0]!=y.shape[0]:
        print(f'{p.name}: row count {len(y)} vs reference {len(x)}; report separately')
        continue
    if not nc: continue
    diff=np.abs(x[nc].to_numpy(dtype=float)-y[nc].to_numpy(dtype=float))
    print(f'{p.name}: maximum numeric absolute difference: {np.nanmax(diff):.6g}')
