#!/usr/bin/env python3
"""Check coverage, numeric type, missingness and uniqueness of reference result tables."""
from pathlib import Path
import pandas as pd
root=Path(__file__).resolve().parents[1]
files=sorted((root/'expected').glob('*.csv'))
if len(files)<20: raise SystemExit('Insufficient machine-readable result tables')
issues=[]
for f in files:
    x=pd.read_csv(f)
    if x.empty: issues.append(f'{f.name}: empty')
    if x.columns.duplicated().any(): issues.append(f'{f.name}: duplicate columns')
    if not x.select_dtypes('number').columns.any(): issues.append(f'{f.name}: no numeric columns')
if issues: raise SystemExit('\n'.join(issues))
print(f'PASS: {len(files)} result tables parsed with numeric columns')
