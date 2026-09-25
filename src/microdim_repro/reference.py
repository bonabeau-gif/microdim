from __future__ import annotations
from pathlib import Path
import pandas as pd
import numpy as np


def compare_csv(actual,reference,key=None,rtol=1e-6,atol=1e-8):
    a=pd.read_csv(actual);r=pd.read_csv(reference)
    if key:
        a=a.sort_values(key).reset_index(drop=True);r=r.sort_values(key).reset_index(drop=True)
    common=[c for c in r.columns if c in a.columns]
    problems=[]
    for c in common:
        if pd.api.types.is_numeric_dtype(r[c]) and pd.api.types.is_numeric_dtype(a[c]):
            x=a[c].to_numpy(float);y=r[c].to_numpy(float)
            if len(x)!=len(y) or not np.allclose(x,y,rtol=rtol,atol=atol,equal_nan=True):problems.append(c)
    return problems
