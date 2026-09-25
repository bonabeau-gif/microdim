import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt

OUT=Path('/mnt/data/calibration_outputs')
OUT.mkdir(exist_ok=True)

COHORTS=[
('Erawijantari gastric cancer',96,0.064,2),('Franzosa IBD',220,0.215,5),('He infants',277,0.091,2),
('Jacobs IBD families',90,0.112,1),('Kang autism',44,-0.015,1),('Wandro preterms',75,0.181,1),
('Wang ESRD',287,0.096,2),('Yachida CRC',347,0.098,2)]
R2_GRID=np.array([0.02,0.05,0.10,0.15,0.20,0.30,0.50])
RANKS=[2,10]; NSEEDS=100; RIDGE=10.0; MAXRANK=12

def outer_train_size(N): return int(round(0.8*N))
def inner_train_approx(n): return int(round(0.8*n))

def make_B(rng,p,q,r,R2,spectrum='equal',ratio=0.7):
    U,_=np.linalg.qr(rng.normal(size=(p,r)))
    V,_=np.linalg.qr(rng.normal(size=(q,r)))
    total_signal=R2/(1-R2)  # mean signal variance per outcome, noise var=1
    if spectrum=='equal':
        s2=np.repeat(q*total_signal/r,r)
    else:
        w=ratio**np.arange(r); s2=q*total_signal*w/w.sum()
    return U[:,:r]@np.diag(np.sqrt(s2))@V[:,:r].T

def folds5(n,seed):
    rng=np.random.default_rng(seed+987654321)
    perm=rng.permutation(n)
    sizes=np.full(5,n//5,dtype=int); sizes[:n%5]+=1
    out=[]; start=0
    allidx=np.arange(n)
    for sz in sizes:
        va=perm[start:start+sz]; mask=np.ones(n,dtype=bool); mask[va]=False; tr=allidx[mask]
        out.append((tr,va)); start+=sz
    return out

def cv_select(X,Y,seed):
    maxrank=min(MAXRANK,X.shape[1],Y.shape[1]); sses=np.zeros(maxrank+1); sse0=0.0
    for tr,va in folds5(len(X),seed):
        Xtr=X[tr]; Ytr=Y[tr]; Xv=X[va]; Yv=Y[va]
        xm=Xtr.mean(0); ym=Ytr.mean(0); Xc=Xtr-xm; Yc=Ytr-ym; Xvc=Xv-xm
        resid=Yv-ym; e0=np.sum(resid*resid); sses[0]+=e0; sse0+=e0
        B=np.linalg.solve(Xc.T@Xc + RIDGE*np.eye(X.shape[1]), Xc.T@Yc)
        _,_,Vt=np.linalg.svd(Xc@B,full_matrices=False)
        coeff=(Xvc@B)@Vt.T
        R=resid.copy()
        for k in range(maxrank):
            R-=np.outer(coeff[:,k],Vt[k])
            sses[k+1]+=np.sum(R*R)
    sel=int(np.argmin(sses)); return sel,1-sses[sel]/sse0

def sim_condition(n,p,q,r,R2,base_seed,spectrum='equal'):
    sels=np.empty(NSEEDS,dtype=int); cvs=np.empty(NSEEDS)
    for j in range(NSEEDS):
        seed=base_seed+j; rng=np.random.default_rng(seed)
        X=rng.normal(size=(n,p)); B=make_B(rng,p,q,r,R2,spectrum); Y=X@B+rng.normal(size=(n,q))
        sels[j],cvs[j]=cv_select(X,Y,seed)
    return sels,cvs

def summarize(sels,cvs,r):
    return dict(median_selected_rank=float(np.median(sels)),mean_selected_rank=float(np.mean(sels)),
                exact_rank_rate=float(np.mean(sels==r)),rank_zero_rate=float(np.mean(sels==0)),
                fraction_true_directions_retained=float(np.mean(np.minimum(sels,r)/r)),
                median_cv_R2=float(np.median(cvs)),mean_cv_R2=float(np.mean(cvs)))

rows=[]
cond=0
for cohort,N,obs_r2,obs_rank in COHORTS:
    n=outer_train_size(N)
    for r in RANKS:
        for R2 in R2_GRID:
            cond+=1; sels,cvs=sim_condition(n,24,24,r,R2,10_000_000+cond*1000)
            d=dict(cohort=cohort,N_total=N,n_rank_selection=n,n_inner_train_approx=inner_train_approx(n),p=24,q=24,
                   true_rank=r,population_R2=R2,per_direction_variance_share=R2/r,spectrum='equal',seeds=NSEEDS,
                   observed_cohort_R2=obs_r2,observed_cohort_rank=obs_rank); d.update(summarize(sels,cvs,r)); rows.append(d)
primary=pd.DataFrame(rows); primary.to_csv(OUT/'sample_matched_primary.csv',index=False)

# Cap sensitivity at three key signal levels only.
rows=[]; cond=0
for cohort,N,obs_r2,obs_rank in COHORTS:
    n=outer_train_size(N); dim=max(2,min(24,n//4))
    for r in RANKS:
        if r>dim:
            rows.append(dict(cohort=cohort,N_total=N,n_rank_selection=n,n_inner_train_approx=inner_train_approx(n),p=dim,q=dim,true_rank=r,
                population_R2=np.nan,per_direction_variance_share=np.nan,spectrum='equal',seeds=0,observed_cohort_R2=obs_r2,
                observed_cohort_rank=obs_rank,note='true rank exceeds preprocessing feature cap'))
            continue
        for R2 in [0.10,0.20,0.50]:
            cond+=1; sels,cvs=sim_condition(n,dim,dim,r,R2,20_000_000+cond*1000)
            d=dict(cohort=cohort,N_total=N,n_rank_selection=n,n_inner_train_approx=inner_train_approx(n),p=dim,q=dim,true_rank=r,
                   population_R2=R2,per_direction_variance_share=R2/r,spectrum='equal',seeds=NSEEDS,observed_cohort_R2=obs_r2,
                   observed_cohort_rank=obs_rank,note=''); d.update(summarize(sels,cvs,r)); rows.append(d)
cap=pd.DataFrame(rows); cap.to_csv(OUT/'sample_and_cap_matched_sensitivity.csv',index=False)

# Spectrum sensitivity at Franzosa-sized outer-training n.
rows=[]; n=outer_train_size(220); cond=0
for spec in ['equal','geometric']:
    for R2 in R2_GRID:
        cond+=1; sels,cvs=sim_condition(n,24,24,10,R2,30_000_000+cond*1000,spec)
        d=dict(size_label='Franzosa-sized',n_rank_selection=n,n_inner_train_approx=inner_train_approx(n),p=24,q=24,true_rank=10,
               population_R2=R2,per_direction_variance_share_equal=R2/10,spectrum=spec,seeds=NSEEDS); d.update(summarize(sels,cvs,10)); rows.append(d)
specdf=pd.DataFrame(rows); specdf.to_csv(OUT/'spectrum_sensitivity.csv',index=False)

# Cohort design and descriptive nearest-grid summary.
cohort_table=pd.DataFrame([dict(cohort=c,N_total=N,n_outer_train_approx=outer_train_size(N),n_inner_train_approx=inner_train_approx(outer_train_size(N)),
  feature_cap_quarter_outer_train=outer_train_size(N)//4,observed_R2=r2,observed_RRR_rank=rr) for c,N,r2,rr in COHORTS])
cohort_table.to_csv(OUT/'cohort_sampling_design.csv',index=False)
near=[]
for cohort,N,obs_r2,obs_rank in COHORTS:
    if obs_r2<=0: continue
    target=float(R2_GRID[np.argmin(np.abs(R2_GRID-obs_r2))])
    for r in RANKS:
        row=primary[(primary.cohort==cohort)&(primary.true_rank==r)&(primary.population_R2==target)].iloc[0]
        near.append(dict(cohort=cohort,N_total=N,n_rank_selection=int(row.n_rank_selection),n_inner_train_approx=int(row.n_inner_train_approx),
             observed_heldout_R2=obs_r2,nearest_population_R2=target,true_rank=r,median_selected_rank=row.median_selected_rank,
             mean_selected_rank=row.mean_selected_rank,exact_rank_rate=row.exact_rank_rate,rank_zero_rate=row.rank_zero_rate,
             fraction_true_directions_retained=row.fraction_true_directions_retained,
             caveat='Observed held-out R2 is not population R2; nearest-grid comparison is descriptive only.'))
near=pd.DataFrame(near); near.to_csv(OUT/'cohort_nearest_signal_comparison.csv',index=False)

# Grid thresholds.
threshold=[]
for cohort,N,obs_r2,obs_rank in COHORTS:
    for r in RANKS:
        sub=primary[(primary.cohort==cohort)&(primary.true_rank==r)].sort_values('population_R2')
        for k in ([2] if r==2 else [1,5,8,10]):
            ok=sub[sub.median_selected_rank>=k]; thr=float(ok.population_R2.iloc[0]) if len(ok) else np.nan
            threshold.append(dict(cohort=cohort,N_total=N,n_rank_selection=outer_train_size(N),n_inner_train_approx=inner_train_approx(outer_train_size(N)),
                true_rank=r,target_median_selected_rank=k,min_grid_population_R2=thr,
                corresponding_per_direction_share=(thr/r if np.isfinite(thr) else np.nan)))
th=pd.DataFrame(threshold); th.to_csv(OUT/'detection_threshold_summary.csv',index=False)

# Figures
fig,ax=plt.subplots(figsize=(8.3,5.2))
for cohort,N,_,_ in COHORTS:
    sub=primary[(primary.cohort==cohort)&(primary.true_rank==10)].sort_values('population_R2')
    ax.plot(sub.population_R2,sub.median_selected_rank,marker='o',label=f'{cohort.split()[0]} (n={outer_train_size(N)})')
ax.axhline(10,linestyle='--',linewidth=1); ax.set_xlabel('Population variance share carried by rank-10 signal'); ax.set_ylabel('Median selected rank')
ax.set_title('Sample size strongly controls recovery of a true rank-10 signal'); ax.set_ylim(-0.3,10.6); ax.legend(fontsize=7,ncol=2,frameon=False)
fig.tight_layout(); fig.savefig(OUT/'figure_rank10_sample_matched.png',dpi=220); plt.close(fig)

fig,ax=plt.subplots(figsize=(7.2,5.0))
for cohort in ['Wandro preterms','Franzosa IBD','Yachida CRC']:
    sub=primary[(primary.cohort==cohort)&(primary.true_rank==10)].sort_values('per_direction_variance_share')
    ax.plot(100*sub.per_direction_variance_share,sub.fraction_true_directions_retained,marker='o',label=f'{cohort} (rank-selection n≈{int(sub.n_rank_selection.iloc[0])})')
ax.set_xlabel('Variance share per true direction (%)'); ax.set_ylabel('Mean fraction of true directions retained'); ax.set_ylim(-0.02,1.02)
ax.set_title('Direction recovery depends on signal per direction and samples per fold'); ax.legend(fontsize=8,frameon=False)
fig.tight_layout(); fig.savefig(OUT/'figure_per_direction_recovery.png',dpi=220); plt.close(fig)

fig,ax=plt.subplots(figsize=(6.6,4.8))
for spec in ['equal','geometric']:
    sub=specdf[specdf.spectrum==spec].sort_values('population_R2')
    ax.plot(sub.population_R2,sub.median_selected_rank,marker='o',label=spec)
ax.axhline(10,linestyle='--',linewidth=1); ax.set_xlabel('Population variance share carried by rank-10 signal'); ax.set_ylabel('Median selected rank')
ax.set_title('Unequal signal across directions increases under-selection'); ax.legend(frameon=False)
fig.tight_layout(); fig.savefig(OUT/'figure_spectrum_sensitivity.png',dpi=220); plt.close(fig)

print('DONE')
print(cohort_table.to_string(index=False))
print('\nRank-10 median>=5 threshold:')
print(th[(th.true_rank==10)&(th.target_median_selected_rank==5)].to_string(index=False))
print('\nNearest signal summary:')
print(near.to_string(index=False))
