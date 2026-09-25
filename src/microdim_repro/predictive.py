from __future__ import annotations
import numpy as np
import pandas as pd
from pathlib import Path
from sklearn.decomposition import PCA

from .preprocess import FoldPreprocessor
from .io import load_cohort
from .config import cohorts


def balanced_group_folds(groups, n_splits=5, seed=0):
    """Randomized, approximately balanced group-wise folds."""
    groups = np.asarray(groups).astype(str)
    uniq, counts = np.unique(groups, return_counts=True)
    if len(uniq) < n_splits:
        raise ValueError(f"Need at least {n_splits} independent groups; got {len(uniq)}")
    rng = np.random.default_rng(seed)
    order = rng.permutation(len(uniq))
    # Allocate largest groups greedily, with random ordering breaking ties.
    order = order[np.argsort(counts[order])[::-1]]
    bins = [[] for _ in range(n_splits)]
    sizes = np.zeros(n_splits, int)
    for j in order:
        k = int(np.argmin(sizes))
        bins[k].append(uniq[j]); sizes[k] += counts[j]
    out=[]
    allidx=np.arange(len(groups))
    for b in bins:
        val=np.flatnonzero(np.isin(groups,b)); mask=np.ones(len(groups),bool); mask[val]=False
        out.append((allidx[mask],val))
    return out


def r2_multivariate(Y, pred, training_mean):
    num=np.sum((Y-pred)**2)
    den=np.sum((Y-training_mean)**2)
    return float(1-num/den) if den>0 else np.nan


def fit_ridge(X,Y,alpha):
    xm=X.mean(0); ym=Y.mean(0); Xc=X-xm; Yc=Y-ym
    B=np.linalg.solve(Xc.T@Xc + alpha*np.eye(X.shape[1]), Xc.T@Yc)
    return xm,ym,B


def predict_ridge(model,X):
    xm,ym,B=model
    return ym+(X-xm)@B


def fit_rrr(X,Y,rank,alpha=10.0):
    xm,ym,B=fit_ridge(X,Y,alpha)
    if rank==0:
        return xm,ym,np.zeros_like(B)
    Xc=X-xm
    _,_,Vt=np.linalg.svd(Xc@B,full_matrices=False)
    k=min(rank,Vt.shape[0])
    P=Vt[:k].T@Vt[:k]
    return xm,ym,B@P


def predict_rrr(model,X): return predict_ridge(model,X)


def fit_pcr(X,Y,k):
    xm=X.mean(0); ym=Y.mean(0); Xc=X-xm; Yc=Y-ym
    if k==0:
        return xm,ym,None,np.zeros((0,Y.shape[1]))
    pca=PCA(n_components=min(k,X.shape[1],max(1,X.shape[0]-1)),svd_solver="full")
    Z=pca.fit_transform(Xc)
    coef=np.linalg.lstsq(Z,Yc,rcond=None)[0]
    return xm,ym,pca,coef


def predict_pcr(model,X):
    xm,ym,pca,coef=model
    if pca is None: return np.tile(ym,(len(X),1))
    return ym+pca.transform(X-xm)@coef


def _inner_select(Xraw,Yraw,groups,model_kind,seed,alpha_grid=None,max_rank=12):
    folds=balanced_group_folds(groups,5,seed)
    if model_kind=="ridge":
        candidates=list(alpha_grid or [1e-3,1e-2,1e-1,1,10,100,1000])
    else:
        # Dimension is bounded fold-by-fold by the fitted transformed matrices.
        candidates=list(range(0,max_rank+1))
    sse=np.zeros(len(candidates)); denom=0.0
    for tr,va in folds:
        pp=FoldPreprocessor().fit(Xraw[tr],Yraw[tr])
        Xtr,Ytr=pp.transform_X(Xraw[tr]),pp.transform_Y(Yraw[tr])
        Xv,Yv=pp.transform_X(Xraw[va]),pp.transform_Y(Yraw[va])
        base=Ytr.mean(0)
        denom += np.sum((Yv-base)**2)
        for ci,c in enumerate(candidates):
            if model_kind=="ridge":
                model=fit_ridge(Xtr,Ytr,c); pred=predict_ridge(model,Xv)
            elif model_kind=="pcr":
                if c>min(Xtr.shape[1],max(0,Xtr.shape[0]-1)):
                    sse[ci]+=np.inf; continue
                model=fit_pcr(Xtr,Ytr,c); pred=predict_pcr(model,Xv)
            else:
                if c>min(Xtr.shape[1],Ytr.shape[1]):
                    sse[ci]+=np.inf; continue
                model=fit_rrr(Xtr,Ytr,c,10.0); pred=predict_rrr(model,Xv)
            sse[ci]+=np.sum((Yv-pred)**2)
    return candidates[int(np.argmin(sse))]


def benchmark_cohort(taxa, metabolites, metadata, group_column="Subject", seed=20260924):
    common=np.asarray(taxa.index.astype(str))
    Xraw=taxa.loc[common].to_numpy(float); Yraw=metabolites.loc[common].to_numpy(float)
    md=metadata.set_index("Sample").loc[common]
    groups=(md[group_column] if group_column in md.columns else md.index.to_series()).astype(str).to_numpy()
    outer=balanced_group_folds(groups,5,seed)
    rows=[]
    for fold,(tr,te) in enumerate(outer):
        pp=FoldPreprocessor().fit(Xraw[tr],Yraw[tr])
        Xtr,Ytr=pp.transform_X(Xraw[tr]),pp.transform_Y(Yraw[tr])
        Xte,Yte=pp.transform_X(Xraw[te]),pp.transform_Y(Yraw[te])
        gtr=groups[tr]
        alpha=_inner_select(Xraw[tr],Yraw[tr],gtr,"ridge",seed+100+fold)
        kpcr=_inner_select(Xraw[tr],Yraw[tr],gtr,"pcr",seed+200+fold)
        krrr=_inner_select(Xraw[tr],Yraw[tr],gtr,"rrr",seed+300+fold)
        specs=[("Ridge",alpha),("PCR",kpcr),("RRR",krrr)]
        for kind,param in specs:
            if kind=="Ridge": pred=predict_ridge(fit_ridge(Xtr,Ytr,param),Xte)
            elif kind=="PCR": pred=predict_pcr(fit_pcr(Xtr,Ytr,param),Xte)
            else: pred=predict_rrr(fit_rrr(Xtr,Ytr,param,10.0),Xte)
            rows.append(dict(fold=fold,model=kind,selected=param,
                             r2=r2_multivariate(Yte,pred,Ytr.mean(0)),
                             n_train=len(tr),n_test=len(te),p=Xtr.shape[1],q=Ytr.shape[1]))
    return pd.DataFrame(rows)


def run_all_cohorts(data_dir, outdir, config_path=None, seed=20260924):
    outdir=Path(outdir); outdir.mkdir(parents=True,exist_ok=True)
    fold_rows=[]; summaries=[]
    for ci,c in enumerate(cohorts(config_path)):
        taxa,mtb,meta=load_cohort(data_dir,c.id)
        df=benchmark_cohort(taxa,mtb,meta,c.group_column,seed+1000*ci)
        df.insert(0,"cohort",c.label); fold_rows.append(df)
        s=df.groupby("model").r2.agg(["mean","std"]).reset_index()
        rrr=df[df.model=="RRR"].selected.median()
        rec={"cohort":c.label,"RRR_dimension":float(rrr)}
        for _,row in s.iterrows():
            rec[f"{row['model']}_mean_R2"]=row["mean"]; rec[f"{row['model']}_sd_R2"]=row["std"]
        summaries.append(rec)
    folds=pd.concat(fold_rows,ignore_index=True); summary=pd.DataFrame(summaries)
    folds.to_csv(outdir/"cohort_fold_results.csv",index=False)
    summary.to_csv(outdir/"table_S4_recomputed.csv",index=False)
    return folds,summary
