#!/usr/bin/env python3
from pathlib import Path
import pandas as pd, numpy as np, sys
ROOT=Path(__file__).resolve().parents[1]
checks=[]
def ck(path,cond,msg): checks.append((path,bool(cond),msg))
h=pd.read_csv(ROOT/'expected/hcom2_geometry.csv'); ck('hcom2_geometry',abs(h.loc[h['matrix']=='original','participation_dimension'].iloc[0]-3.270)<1e-12,'D_PR=3.270')
n=pd.read_csv(ROOT/'expected/hcom2_nmf_threshold.csv'); ck('hcom2_nmf',int(n.loc[(n.retained_improvement>=.95)&((.773-n.p3_pearson)<=.03),'rank'].min())==52,'selected-seed threshold 52')
s=pd.read_csv(ROOT/'expected/hcom2_support_strength.csv'); ck('hcom2_support',int(s.loc[(s.retained_improvement>=.95)&((.773-s.p3_pearson)<=.03),'latent_dimension'].min())==6,'support latent threshold 6')
ho=pd.read_csv(ROOT/'expected/ho_support_controls.csv'); m0=float(ho.loc[ho.support_model=='Published support','primary_error'].iloc[0]); sh=float(ho.loc[ho.support_model=='Species-label shuffle','primary_error'].iloc[0]); ck('ho_shuffle',(sh-m0)/sh>.276-1e-3,'>27.6% lower error'); ck('ho_wins',abs(float(ho.loc[ho.support_model=='Species-label shuffle','assembly_win_fraction'].iloc[0])-56/64)<1e-12,'56/64 wins')
cal=pd.read_csv(ROOT/'expected/rank_thresholds.csv'); ck('rank_threshold',cal.iloc[0].per_direction_share==.05 and cal.iloc[-1].per_direction_share==.01,'5% to 1%')
for name,ok,msg in checks: print(('PASS' if ok else 'FAIL'),name,msg)
if not all(x[1] for x in checks): sys.exit(1)
