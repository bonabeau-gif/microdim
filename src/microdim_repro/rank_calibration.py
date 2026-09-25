from __future__ import annotations
import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt

COHORTS=[
('Erawijantari gastric cancer',96,0.064,2),('Franzosa IBD',220,0.215,5),('He infants',277,0.091,2),
('Jacobs IBD families',90,0.112,1),('Kang autism',44,-0.015,1),('Wandro preterms',75,0.181,1),
('Wang ESRD',287,0.096,2),('Yachida CRC',347,0.098,2)]
R2_GRID=np.array([0.02,0.05,0.10,0.15,0.20,0.30,0.50])
RANKS=[2,10]; NSEEDS=100; RIDGE=10.0; MAXRANK=12

def outer_train_size(N): return int(round(0.8*N))
def inner_train_approx(n): return int(round(0.8*n))

def make_B(rng,p,q,r,R2,spectrum='equal',ratio=0.7):
    U,_=np.linalg.qr(rng.normal(size=(p,r))); V,_=np.linalg.qr(rng.normal(size=(q,r)))
    total_signal=R2/(1-R2)
    if spectrum=='equal': s2=np.repeat(q*total_signal/r,r)
    else:
        w=ratio**np.arange(r); s2=q*total_signal*w/w.sum()
    return U[:,:r]@np.diag(np.sqrt(s2))@V[:,:r].T

def folds5(n,seed):
    rng=np.random.default_rng(seed+987654321); perm=rng.permutation(n)
    sizes=np.full(5,n//5,dtype=int); sizes[:n%5]+=1
    out=[]; start=0; allidx=np.arange(n)
    for sz in sizes:
        va=perm[start:start+sz]; mask=np.ones(n,dtype=bool); mask[va]=False
        out.append((allidx[mask],va)); start+=sz
    return out

def cv_select(X,Y,seed,return_projection=False):
    maxrank=min(MAXRANK,X.shape[1],Y.shape[1]); sses=np.zeros(maxrank+1); sse0=0.0
    for tr,va in folds5(len(X),seed):
        Xtr=X[tr];Ytr=Y[tr];Xv=X[va];Yv=Y[va]
        xm=Xtr.mean(0);ym=Ytr.mean(0);Xc=Xtr-xm;Yc=Ytr-ym;Xvc=Xv-xm
        resid=Yv-ym;e0=np.sum(resid*resid);sses[0]+=e0;sse0+=e0
        B=np.linalg.solve(Xc.T@Xc+RIDGE*np.eye(X.shape[1]),Xc.T@Yc)
        _,_,Vt=np.linalg.svd(Xc@B,full_matrices=False);coeff=(Xvc@B)@Vt.T;R=resid.copy()
        for k in range(maxrank):
            R-=np.outer(coeff[:,k],Vt[k]);sses[k+1]+=np.sum(R*R)
    sel=int(np.argmin(sses)); r2=1-sses[sel]/sse0
    if not return_projection: return sel,r2
    xm=X.mean(0);ym=Y.mean(0);Xc=X-xm;Yc=Y-ym
    Bhat=np.linalg.solve(Xc.T@Xc+RIDGE*np.eye(X.shape[1]),Xc.T@Yc)
    _,_,Vt=np.linalg.svd(Xc@Bhat,full_matrices=False)
    P=Vt[:sel].T@Vt[:sel] if sel>0 else np.zeros((Y.shape[1],Y.shape[1]))
    return sel,r2,P

def sim_condition(n,p,q,r,R2,base_seed,spectrum='equal',nseeds=NSEEDS,subspace=False):
    sels=np.empty(nseeds,dtype=int);cvs=np.empty(nseeds);capt=np.full(nseeds,np.nan)
    for j in range(nseeds):
        seed=base_seed+j;rng=np.random.default_rng(seed);X=rng.normal(size=(n,p));B=make_B(rng,p,q,r,R2,spectrum);Y=X@B+rng.normal(size=(n,q))
        if subspace:
            sels[j],cvs[j],P=cv_select(X,Y,seed,True);capt[j]=np.sum((B@P)**2)/np.sum(B**2)
        else: sels[j],cvs[j]=cv_select(X,Y,seed)
    return sels,cvs,capt

def summarize(sels,cvs,r,capt=None):
    d=dict(median_selected_rank=float(np.median(sels)),mean_selected_rank=float(np.mean(sels)),exact_rank_rate=float(np.mean(sels==r)),
           rank_zero_rate=float(np.mean(sels==0)),fraction_true_directions_retained=float(np.mean(np.minimum(sels,r)/r)),
           median_cv_R2=float(np.median(cvs)),mean_cv_R2=float(np.mean(cvs)))
    if capt is not None and np.isfinite(capt).any(): d["signal_variance_captured"]=float(np.nanmean(capt))
    return d

def run(outdir):
    outdir=Path(outdir);outdir.mkdir(parents=True,exist_ok=True)
    rows=[];cond=0
    for cohort,N,obs_r2,obs_rank in COHORTS:
        n=outer_train_size(N)
        for r in RANKS:
            for R2 in R2_GRID:
                cond+=1;sels,cvs,_=sim_condition(n,24,24,r,R2,10_000_000+cond*1000)
                d=dict(cohort=cohort,N_total=N,n_rank_selection=n,n_inner_train_approx=inner_train_approx(n),p=24,q=24,true_rank=r,
                       population_R2=R2,per_direction_variance_share=R2/r,spectrum='equal',seeds=NSEEDS,observed_cohort_R2=obs_r2,observed_cohort_rank=obs_rank)
                d.update(summarize(sels,cvs,r));rows.append(d)
    primary=pd.DataFrame(rows);primary.to_csv(outdir/'sample_matched_primary.csv',index=False)

    rows=[];cond=0
    for cohort,N,obs_r2,obs_rank in COHORTS:
        n=outer_train_size(N);dim=max(2,min(24,n//4))
        for r in RANKS:
            if r>dim:
                rows.append(dict(cohort=cohort,N_total=N,n_rank_selection=n,n_inner_train_approx=inner_train_approx(n),p=dim,q=dim,true_rank=r,
                    population_R2=np.nan,per_direction_variance_share=np.nan,spectrum='equal',seeds=0,observed_cohort_R2=obs_r2,observed_cohort_rank=obs_rank,note='true rank exceeds preprocessing feature cap'))
                continue
            for R2 in [0.10,0.20,0.50]:
                cond+=1;sels,cvs,_=sim_condition(n,dim,dim,r,R2,20_000_000+cond*1000)
                d=dict(cohort=cohort,N_total=N,n_rank_selection=n,n_inner_train_approx=inner_train_approx(n),p=dim,q=dim,true_rank=r,population_R2=R2,
                       per_direction_variance_share=R2/r,spectrum='equal',seeds=NSEEDS,observed_cohort_R2=obs_r2,observed_cohort_rank=obs_rank,note='')
                d.update(summarize(sels,cvs,r));rows.append(d)
    cap=pd.DataFrame(rows);cap.to_csv(outdir/'sample_and_cap_matched_sensitivity.csv',index=False)

    rows=[];n=outer_train_size(220);cond=0
    for spec in ['equal','geometric']:
        for R2 in R2_GRID:
            cond+=1;sels,cvs,_=sim_condition(n,24,24,10,R2,30_000_000+cond*1000,spec)
            d=dict(size_label='Franzosa-sized',n_rank_selection=n,n_inner_train_approx=inner_train_approx(n),p=24,q=24,true_rank=10,
                   population_R2=R2,per_direction_variance_share_equal=R2/10,spectrum=spec,seeds=NSEEDS);d.update(summarize(sels,cvs,10));rows.append(d)
    specdf=pd.DataFrame(rows);specdf.to_csv(outdir/'spectrum_sensitivity.csv',index=False)

    cohort_table=pd.DataFrame([dict(cohort=c,N_total=N,n_outer_train_approx=outer_train_size(N),n_inner_train_approx=inner_train_approx(outer_train_size(N)),
      feature_cap_quarter_outer_train=outer_train_size(N)//4,observed_R2=r2,observed_RRR_rank=rr) for c,N,r2,rr in COHORTS])
    cohort_table.to_csv(outdir/'cohort_sampling_design.csv',index=False)

    threshold=[]
    for cohort,N,obs_r2,obs_rank in COHORTS:
        for r in RANKS:
            sub=primary[(primary.cohort==cohort)&(primary.true_rank==r)].sort_values('population_R2')
            for k in ([2] if r==2 else [1,5,8,10]):
                ok=sub[sub.median_selected_rank>=k];thr=float(ok.population_R2.iloc[0]) if len(ok) else np.nan
                threshold.append(dict(cohort=cohort,N_total=N,n_rank_selection=outer_train_size(N),n_inner_train_approx=inner_train_approx(outer_train_size(N)),
                    true_rank=r,target_median_selected_rank=k,min_grid_population_R2=thr,corresponding_per_direction_share=(thr/r if np.isfinite(thr) else np.nan)))
    pd.DataFrame(threshold).to_csv(outdir/'detection_threshold_summary.csv',index=False)

    # Table S8: use the primary-condition seed bases used above.
    selected=[(60,.10),(176,.10),(278,.10),(176,.15),(176,.20)]
    s8=[]
    # Determine primary condition seed base by reusing the condition numbering.
    cond_map={};cond=0
    for cohort,N,_,_ in COHORTS:
        n0=outer_train_size(N)
        for r in RANKS:
            for R2 in R2_GRID:
                cond+=1;cond_map[(n0,r,float(R2))]=10_000_000+cond*1000
    for n,R2 in selected:
        sels,cvs,capt=sim_condition(n,24,24,10,R2,cond_map[(n,10,float(R2))],subspace=True)
        s8.append(dict(n_rank_selection=n,population_R2=R2,per_direction_share=R2/10,
                       dimensions_retained=float(np.mean(np.minimum(sels,10)/10)),signal_variance_captured=float(np.mean(capt))))
    pd.DataFrame(s8).to_csv(outdir/'table_S8_subspace_recovery.csv',index=False)
    return primary,cap,specdf
