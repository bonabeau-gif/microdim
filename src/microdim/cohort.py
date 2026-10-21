from pathlib import Path
import numpy as np, pandas as pd
from sklearn.model_selection import GroupKFold, KFold
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import Ridge
from .metrics import clr, multivariate_r2
from .rrr import rrr_predict, pcr_predict

def load_tables(folder):
    f=Path(folder); g=pd.read_csv(f/'genera.tsv',sep='\t',index_col=0); m=pd.read_csv(f/'mtb.tsv',sep='\t',index_col=0,na_values=['NA','NaN','']); meta=pd.read_csv(f/'metadata.tsv',sep='\t',index_col='Sample')
    common=g.index.intersection(m.index).intersection(meta.index); return g.loc[common],m.loc[common],meta.loc[common]
def prep_fold(g,m,tr,te):
    gt,mt=g.iloc[tr],m.iloc[tr]
    taxa=(gt>0).mean(0)>=.25; mets=mt.notna().mean(0)>=.75
    taxa_names=gt.columns[taxa]; met_names=mt.columns[mets]
    cap=max(1,len(tr)//4)
    if len(taxa_names)>cap: taxa_names=gt[taxa_names].mean(0).sort_values(ascending=False).index[:cap]
    if len(met_names)>cap: met_names=mt[met_names].var(0,skipna=True).sort_values(ascending=False).index[:cap]
    minpos=gt[taxa_names].where(gt[taxa_names]>0).min().min(); pc=float(minpos/2) if np.isfinite(minpos) and minpos>0 else 1e-9
    Xtr=clr(gt[taxa_names].to_numpy(),pc); Xte=clr(g.iloc[te][taxa_names].to_numpy(),pc)
    med=mt[met_names].median(0); Ytr=np.log1p(mt[met_names].fillna(med).clip(lower=0).to_numpy()); Yte=np.log1p(m.iloc[te][met_names].fillna(med).clip(lower=0).to_numpy())
    sx=StandardScaler().fit(Xtr); sy=StandardScaler().fit(Ytr); return sx.transform(Xtr),sx.transform(Xte),sy.transform(Ytr),sy.transform(Yte)
def inner_score(model,X,Y,groups=None,seed=1):
    if groups is not None and len(np.unique(groups))>=5: cv=GroupKFold(5).split(X,Y,groups)
    else: cv=KFold(5,shuffle=True,random_state=seed).split(X,Y)
    vals=[]
    for tr,va in cv:
        pred=model(X[tr],Y[tr],X[va]); vals.append(multivariate_r2(Y[va],pred,Y[tr].mean(0)))
    return np.mean(vals)
def benchmark(folder,group_col='Subject',seed=20260925):
    g,m,meta=load_tables(folder); groups=meta[group_col].astype(str).to_numpy() if group_col in meta else np.arange(len(meta)).astype(str)
    outer=list(GroupKFold(5).split(g,np.zeros(len(g)),groups)) if len(np.unique(groups))>=5 else list(KFold(5,shuffle=True,random_state=seed).split(g))
    rows=[]
    for fold,(tr,te) in enumerate(outer):
        Xtr,Xte,Ytr,Yte=prep_fold(g,m,tr,te); gr=groups[tr]
        # ridge
        alphas=[.01,.1,1,10,100]; a=max(alphas,key=lambda z:inner_score(lambda a,b,c:Ridge(alpha=z).fit(a,b).predict(c),Xtr,Ytr,gr,seed))
        pr=Ridge(alpha=a).fit(Xtr,Ytr).predict(Xte)
        # PCR
        comps=list(range(1,min(12,Xtr.shape[1],len(Xtr)-1)+1)); nc=max(comps,key=lambda z:inner_score(lambda a,b,c:pcr_predict(a,b,c,z),Xtr,Ytr,gr,seed))
        pp=pcr_predict(Xtr,Ytr,Xte,nc)
        # RRR fixed ridge 10
        ranks=list(range(0,min(12,Xtr.shape[1],Ytr.shape[1])+1)); rr=max(ranks,key=lambda z:inner_score(lambda a,b,c:rrr_predict(a,b,c,z,10),Xtr,Ytr,gr,seed))
        pq=rrr_predict(Xtr,Ytr,Xte,rr,10)
        mu=Ytr.mean(0); rows.append({'fold':fold,'rrr':multivariate_r2(Yte,pq,mu),'pcr':multivariate_r2(Yte,pp,mu),'ridge':multivariate_r2(Yte,pr,mu),'rrr_rank':rr})
    return pd.DataFrame(rows)
