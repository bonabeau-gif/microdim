#!/usr/bin/env python3
exec(open(__file__.replace('run_cohort_benchmark.py','_bootstrap.py')).read())
import json,pandas as pd,numpy as np
from microdim.cohort import benchmark
cfg=json.loads((ROOT/'config/cohorts.json').read_text()); rows=[]
for key,c in cfg.items():
    f=ROOT/'data'/'third_party'/'borenstein'/key
    if not f.exists(): raise SystemExit('Run scripts/fetch_sources.py --cohorts first')
    d=benchmark(f,c['group']); rows.append({'cohort':c['label'],'rrr_mean':d.rrr.mean(),'rrr_sd':d.rrr.std(ddof=0),'pcr_mean':d.pcr.mean(),'pcr_sd':d.pcr.std(ddof=0),'ridge_mean':d.ridge.mean(),'ridge_sd':d.ridge.std(ddof=0),'rrr_dimension':float(d.rrr_rank.median())})
out=pd.DataFrame(rows); (ROOT/'results').mkdir(exist_ok=True); out.to_csv(ROOT/'results/cohort_benchmark.csv',index=False); print(out.to_string(index=False))
