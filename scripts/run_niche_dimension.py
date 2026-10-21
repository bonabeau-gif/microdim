#!/usr/bin/env python3
exec(open(__file__.replace('run_niche_dimension.py','_bootstrap.py')).read())
import json,pandas as pd,numpy as np
from microdim.niche import eta_dimension, metabolite_decay_dimension
cfg=json.loads((ROOT/'config/cohorts.json').read_text()); rows=[]; rng=np.random.default_rng(20260925)
for key,c in cfg.items():
 f=ROOT/'data'/'third_party'/'borenstein'/key; g=pd.read_csv(f/'genera.tsv',sep='\t',index_col=0); m=pd.read_csv(f/'mtb.tsv',sep='\t',index_col=0,na_values=['NA','']); common=g.index.intersection(m.index); g=g.loc[common]; m=m.loc[common]
 keep=g.mean(0)>=.001; g=g.loc[:,keep]; nsub=50 if len(g)>=50 else int(np.floor(.70*len(g))); vals=[]
 for b in range(100):
  ix=rng.choice(len(g),nsub,replace=False); P=g.iloc[ix].to_numpy(float); P=P/np.maximum(P.sum(1,keepdims=True),1e-15); eta,_=eta_dimension(P,5,seed=20260925+b); vals.append(eta)
 X=m.fillna(m.median()).clip(lower=0).to_numpy(float); rows.append({'cohort':c['label'].split()[0],'samples':len(g),'genera':g.shape[1],'subsample':nsub,'niche_mean':np.mean(vals),'niche_sd':np.std(vals),'metabolite_space_dimension':metabolite_decay_dimension(X)})
out=pd.DataFrame(rows); (ROOT/'results').mkdir(exist_ok=True)
out.to_csv(ROOT/'results/niche_dimensions.csv',index=False); print(out.to_string(index=False))
