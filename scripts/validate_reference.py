"""Lightweight audit of bundled manuscript reference results."""
from pathlib import Path
import pandas as pd
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
R=ROOT/'reference_results'

checks=[]
s4=pd.read_csv(R/'table_S4_cohort_benchmark.csv')
checks.append(('Table S4 has eight cohorts',len(s4)==8))
checks.append(('RRR dimensions are 1--5',s4.RRR_dimension.between(1,5).all()))
s5=pd.read_csv(R/'table_S5_niche_dimensions.csv')
checks.append(('Table S5 has eight cohorts',len(s5)==8))
checks.append(('Reported eta range',np.isclose(s5.eta_D_mean.min(),1.43) and np.isclose(s5.eta_D_mean.max(),4.05)))
s9=pd.read_csv(R/'table_S9_micom_dimensions.csv')
checks.append(('MICOM has four communities',len(s9)==4))
checks.append(('MICOM dPR range',s9.flux_participation_rank.min()>1.5 and s9.flux_participation_rank.max()<2.7))
cal=pd.read_csv(R/'sample_matched_primary.csv')
checks.append(('Primary calibration has 112 conditions',len(cal)==112))
for name,ok in checks:
    print(('PASS' if ok else 'FAIL')+': '+name)
if not all(ok for _,ok in checks): raise SystemExit(1)
