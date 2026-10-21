#!/usr/bin/env python3
exec(open(__file__.replace('run_consumer_resource_verification.py','_bootstrap.py')).read())
import numpy as np,pandas as pd
from scipy.integrate import solve_ivp
from microdim.cr_theory import interaction,qss,rhs_full,rhs_red
rng=np.random.default_rng(20260925); rows=[]
for eps in [16,32,64,128,256]:
 errs=[]
 for seed in range(30):
  r=np.random.default_rng(20260925+seed); F=r.uniform(.2,1,(6,3)); G=r.uniform(.2,1,(5,2)); H=r.uniform(.1,.6,(3,2)); C=F@H@G.T; y=r.uniform(.5,1.5,5); a=r.uniform(.8,1.2,5); N0=r.uniform(.01,.03,6); rho=C.T@N0+a*r.uniform(.5,1.,5); m=.5*(C@(y*qss(rho,a,C,N0))); R0=qss(rho,a,C,N0); z0=np.r_[N0,R0]; tt=np.linspace(0,1,101); f=solve_ivp(rhs_full,[0,1],z0,t_eval=tt,args=(C,y,m,rho,a,eps),rtol=1e-10,atol=1e-12); g=solve_ivp(rhs_red,[0,1],N0,t_eval=tt,args=(C,y,m,rho,a),rtol=1e-10,atol=1e-12); errs.append(np.max(np.abs(np.log(f.y[:6])-np.log(g.y))))
 rows.append([eps,np.median(errs),np.quantile(errs,.025),np.quantile(errs,.975)])
out=pd.DataFrame(rows,columns=['timescale_factor','median_error','q025','q975']); (ROOT/'results').mkdir(exist_ok=True)
out.to_csv(ROOT/'results/cr_timescale.csv',index=False); print(out.to_string(index=False))
slope=float(np.polyfit(np.log(out.timescale_factor),np.log(out.median_error),1)[0])
print(f'Log-log slope of trajectory tracking error: {slope:.4f}; expected asymptotic slope: -1')
