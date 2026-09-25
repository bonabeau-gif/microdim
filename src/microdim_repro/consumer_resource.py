from __future__ import annotations
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.integrate import solve_ivp

ETAS=np.array([16.,32.,64.,128.,256.])

def generate_system(seed, uptake_scale=0.365):
    """Generate one deterministic feature-factorized consumer-resource system.

    The manuscript retained the dimensions and qualitative random construction,
    but not the exact original random seed/range file. These ranges are therefore
    an explicit reconstruction, not a claim about lost historical settings.
    """
    rng=np.random.default_rng(seed);S,P,K,L=6,5,3,2
    Phi_n=rng.uniform(.2,1.0,(S,K));Phi_R=rng.uniform(.2,1.0,(P,L))
    C=rng.uniform(.05,.2,(K,L))*uptake_scale;c=Phi_n@C@Phi_R.T
    v=rng.uniform(.7,1.3,P);s=rng.uniform(.8,1.2,P)
    mu=rng.uniform(.03,.07,K);m=Phi_n@mu
    rho=rng.uniform(.8,1.2,L);r=Phi_R@rho+0.8
    n0=rng.uniform(.03,.08,S)
    return dict(S=S,P=P,K=K,L=L,Phi_n=Phi_n,Phi_R=Phi_R,C=C,c=c,v=v,s=s,mu=mu,m=m,rho=rho,r=r,n0=n0)

def reduced_parameters(sys):
    c,v,s,r,m=sys['c'],sys['v'],sys['s'],sys['r'],sys['m']
    growth=c@(v*r/s)-m
    A=-c@np.diag(v/s)@c.T
    M=-sys['C']@sys['Phi_R'].T@np.diag(v/s)@sys['Phi_R']@sys['C'].T
    Af=sys['Phi_n']@M@sys['Phi_n'].T
    return growth,A,M,Af

def simulate_one(seed,etas=ETAS,T=3.0,uptake_scale=0.365):
    sys=generate_system(seed,uptake_scale);S=sys['S'];c=sys['c'];v=sys['v'];s=sys['s'];r=sys['r'];m=sys['m'];n0=sys['n0']
    Req=lambda n:(r-c.T@n)/s
    R0=Req(n0);growth,A,M,Af=reduced_parameters(sys)
    t=np.linspace(0,T,101)
    def red(_,n):return n*(growth+A@n)
    sr=solve_ivp(red,(0,T),n0,t_eval=t,rtol=1e-10,atol=1e-12)
    errors=[];minR=[]
    for eta in etas:
        def full(_,y):
            n=y[:S];R=y[S:]
            return np.r_[n*(c@(v*R)-m),eta*R*(r-s*R-c.T@n)]
        sf=solve_ivp(full,(0,T),np.r_[n0,R0],t_eval=t,rtol=1e-10,atol=1e-12)
        errors.append(np.max(np.abs(np.log(sf.y[:S])-np.log(sr.y))))
        minR.append(sf.y[S:].min())
    relfact=np.linalg.norm(A-Af,'fro')/max(np.linalg.norm(A,'fro'),1e-300)
    return dict(seed=seed,errors=np.array(errors),min_resource=float(min(minR)),factorization_relative_error=float(relfact),
                numerical_rank=int(np.linalg.matrix_rank(A,tol=1e-10)),largest_eigenvalue=float(np.linalg.eigvalsh(A).max()))

def run(outdir,nseeds=30):
    outdir=Path(outdir);outdir.mkdir(parents=True,exist_ok=True)
    sims=[simulate_one(i) for i in range(nseeds)]
    rows=[]
    for s in sims:
        for eta,e in zip(ETAS,s['errors']): rows.append(dict(seed=s['seed'],eta=eta,max_log_abundance_error=e,min_resource=s['min_resource'],
            factorization_relative_error=s['factorization_relative_error'],numerical_rank=s['numerical_rank'],largest_eigenvalue=s['largest_eigenvalue']))
    df=pd.DataFrame(rows);df.to_csv(outdir/'fast_resource_reconstructed_seed_results.csv',index=False)
    curves=df.groupby('eta').max_log_abundance_error.agg(median='median',q025=lambda x:np.quantile(x,.025),q975=lambda x:np.quantile(x,.975)).reset_index()
    curves.to_csv(outdir/'fast_resource_reconstructed_curve.csv',index=False)
    slopes=[]
    for seed,g in df.groupby('seed'):
        slopes.append(np.polyfit(np.log(g.eta),np.log(g.max_log_abundance_error),1)[0])
    summary=pd.DataFrame([dict(nseeds=nseeds,max_factorization_relative_error=df.factorization_relative_error.max(),
        ranks=','.join(map(str,sorted(df.numerical_rank.unique()))),largest_eigenvalue_max=df.largest_eigenvalue.max(),
        minimum_resource=df.min_resource.min(),median_error_eta16=float(curves.loc[curves.eta==16,'median'].iloc[0]),
        median_error_eta256=float(curves.loc[curves.eta==256,'median'].iloc[0]),median_loglog_slope=float(np.median(slopes)),
        slope_q025=float(np.quantile(slopes,.025)),slope_q975=float(np.quantile(slopes,.975)))])
    summary.to_csv(outdir/'fast_resource_reconstructed_summary.csv',index=False)
    return df,curves,summary
