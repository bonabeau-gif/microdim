from __future__ import annotations
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.special import logsumexp
from scipy.stats import spearmanr

from .io import load_cohort
from .config import cohorts


def _normalize_compositions(X):
    X=np.asarray(X,float); X=np.nan_to_num(X,nan=0.0); X[X<0]=0
    rs=X.sum(1,keepdims=True); keep=rs[:,0]>0
    return X[keep]/rs[keep], keep


def fit_lowrank_multinomial(q,K,seed=0,max_iter=20000,tol=1e-2,lr=0.05):
    """Fit the Srinivasan-Plata-Dixit low-rank multinomial model.

    The paper specifies the likelihood, gradients and stopping rule but not a
    unique initialization/learning-rate schedule. This implementation uses a
    deterministic small-normal initialization and normalized gradient ascent
    with backtracking. It is therefore an executable reconstruction of the
    published estimator, while the manuscript's frozen Table S5 remains the
    regression target supplied in reference_results/.
    """
    q=np.asarray(q,float); S,O=q.shape
    rng=np.random.default_rng(seed)
    z=rng.normal(scale=0.05,size=(S,K)); theta=rng.normal(scale=0.05,size=(K,O))

    def objective(z,theta):
        logits=z@theta
        return float(np.sum(q*(logits-logsumexp(logits,axis=1,keepdims=True))))
    old=objective(z,theta)
    for _ in range(max_iter):
        logits=z@theta
        x=np.exp(logits-logsumexp(logits,axis=1,keepdims=True))
        resid=q-x
        gz=resid@theta.T
        gt=z.T@resid
        rel=(np.sum(np.abs(gz))/(np.sum(np.abs(z))+1e-12)+
             np.sum(np.abs(gt))/(np.sum(np.abs(theta))+1e-12))
        if rel<=tol: break
        # Normalize by sample/taxon count so lr is scale-stable.
        gz/=max(1,O); gt/=max(1,S)
        step=lr
        accepted=False
        for _bt in range(12):
            zn=z+step*gz; tn=theta+step*gt
            new=objective(zn,tn)
            if new>=old-1e-12:
                z,theta,old=zn,tn,new; accepted=True; break
            step*=0.5
        if not accepted: break
        # Gauge stabilization: center z columns; preserves model after intercept-free refit approximately.
        # Do not impose arbitrary rotations; only clip extreme values for numerical safety.
        z=np.clip(z,-50,50); theta=np.clip(theta,-50,50)
    logits=z@theta
    x=np.exp(logits-logsumexp(logits,axis=1,keepdims=True))
    return z,theta,x


def mean_kl(q,x):
    mask=q>0
    vals=np.zeros_like(q)
    vals[mask]=q[mask]*(np.log(q[mask])-np.log(np.maximum(x[mask],1e-300)))
    return float(vals.sum(1).mean())


def eta_from_kl(Ks,kls):
    Ks=np.asarray(Ks,float); kls=np.asarray(kls,float)
    slope,intercept=np.polyfit(Ks,np.log(np.maximum(kls,1e-300)),1)
    return float(-1/slope),float(slope),float(intercept)


def estimate_eta(q,seed=0,Kmax=5,**fit_kwargs):
    kls=[]
    for K in range(1,Kmax+1):
        _,_,x=fit_lowrank_multinomial(q,K,seed=seed+1009*K,**fit_kwargs)
        kls.append(mean_kl(q,x))
    eta,slope,intercept=eta_from_kl(np.arange(1,Kmax+1),kls)
    return eta,np.array(kls),slope,intercept


def prepare_abundance_for_niche(taxa,abundance_filter=0.001):
    X=taxa.to_numpy(float)
    X=np.nan_to_num(X,nan=0.0); X[X<0]=0
    rs=X.sum(1,keepdims=True); rs[rs==0]=1
    rel=X/rs
    mean=rel.mean(0)
    keep=mean>=abundance_filter
    rel=rel[:,keep]; rel/=rel.sum(1,keepdims=True)
    return rel,keep


def estimate_cohort_niche(taxa,n_subsamples=100,seed=20260924,abundance_filter=0.001):
    rel,keep=prepare_abundance_for_niche(taxa,abundance_filter)
    n=len(rel); size=50 if n>=50 else max(2,int(round(0.70*n)))
    rng=np.random.default_rng(seed)
    rows=[]
    for b in range(n_subsamples):
        idx=rng.choice(n,size=size,replace=False)
        eta,kls,_,_=estimate_eta(rel[idx],seed=seed+100000*b)
        row={"replicate":b,"eta_D":eta}
        row.update({f"KL_K{k+1}":v for k,v in enumerate(kls)})
        rows.append(row)
    return pd.DataFrame(rows), int(keep.sum()), size


def metabolite_space_dimension(mtb):
    Y=mtb.to_numpy(float)
    vals=Y[np.isfinite(Y)&(Y>0)]
    offset=float(np.quantile(vals,0.10)) if vals.size else 1e-12
    Y=np.where(np.isfinite(Y),Y,np.nan)
    med=np.nanmedian(Y,axis=0); inds=np.where(~np.isfinite(Y)); Y[inds]=med[inds[1]]
    Y=np.maximum(Y+offset,1e-300); Y=np.log(Y)
    sd=Y.std(0); sd[sd==0]=1; Y=(Y-Y.mean(0))/sd
    s=np.linalg.svd(Y,compute_uv=False,full_matrices=False)
    k=np.arange(1,len(s)+1)
    # Fit leading spectrum as in Appendix B; use all nonzero singular values.
    nz=s>0; slope,_=np.polyfit(k[nz],np.log(s[nz]),1)
    return float(-1/slope)


def run_all_niche(data_dir,outdir,config_path=None,n_subsamples=100,seed=20260924,abundance_filter=0.001):
    outdir=Path(outdir); outdir.mkdir(parents=True,exist_ok=True)
    summary=[]; reps=[]
    for ci,c in enumerate(cohorts(config_path)):
        taxa,mtb,_=load_cohort(data_dir,c.id)
        rr,ng,ns=estimate_cohort_niche(taxa,n_subsamples,seed+10000*ci,abundance_filter)
        rr.insert(0,"cohort",c.label); reps.append(rr)
        summary.append(dict(cohort=c.label,samples=len(taxa),genera_retained=ng,samples_used=ns,
                            eta_D_mean=rr.eta_D.mean(),eta_D_sd=rr.eta_D.std(ddof=1),
                            metabolite_space_dimension=metabolite_space_dimension(mtb)))
    reps=pd.concat(reps,ignore_index=True); sm=pd.DataFrame(summary)
    reps.to_csv(outdir/"niche_subsamples.csv",index=False); sm.to_csv(outdir/"table_S5_recomputed.csv",index=False)
    return reps,sm


def run_sensitivity(data_dir,outdir,config_path=None,n_subsamples=100,seed=20260924):
    """Re-run the manuscript-style sensitivity checks.

    Historical variant seed integers were not retained. The variants here use
    three explicitly documented seeds and abundance filters 0.5% and 0%, so the
    workflow can be independently stress-tested. The exact manuscript-reported
    sensitivity correlations are retained in reference_results/.
    """
    outdir=Path(outdir);outdir.mkdir(parents=True,exist_ok=True)
    rrr={c.label:r for c,r in zip(cohorts(config_path),[2,5,2,1,1,1,2,2])}
    variants=[('seed_A',seed,0.001),('seed_B',seed+1,0.001),('seed_C',seed+2,0.001),('filter_0.005',seed,0.005),('filter_none',seed,0.0)]
    rows=[]
    for name,s,filt in variants:
        _,sm=run_all_niche(data_dir,outdir/name,config_path,n_subsamples,s,filt)
        x=sm.set_index('cohort').eta_D_mean
        y=pd.Series(rrr)
        common=x.index.intersection(y.index)
        rho,p=spearmanr(x.loc[common],y.loc[common])
        rows.append(dict(variant=name,seed=s,abundance_filter=filt,n=len(common),spearman_eta_vs_rrr=rho,p_value=p,
                         median_eta=float(x.median())))
    res=pd.DataFrame(rows);res.to_csv(outdir/'niche_sensitivity_recomputed.csv',index=False)
    return res
