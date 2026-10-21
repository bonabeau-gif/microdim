#!/usr/bin/env python3
exec(open(__file__.replace('run_rank_calibration.py','_bootstrap.py')).read())
import argparse,pandas as pd
from microdim.calibration import grid
ap=argparse.ArgumentParser(); ap.add_argument('--smoke',action='store_true'); a=ap.parse_args(); ns=[60,176] if a.smoke else [35,60,72,77,176,222,230,278]; r2=[.1,.2] if a.smoke else [.10,.15,.20,.30,.50]
rows=grid(ns,r2,reps=5 if a.smoke else 100); out=pd.DataFrame(rows,columns=['n','population_r2','median_rank','mean_rank']); (ROOT/'results').mkdir(exist_ok=True)
out.to_csv(ROOT/'results/rank_calibration.csv',index=False); print(out.to_string(index=False))
