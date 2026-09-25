from __future__ import annotations
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import spearmanr

EPS=np.finfo(float).eps

def rank_metrics(mat,rel_tol=1e-8):
    a=np.asarray(mat,float)
    if a.size==0:return dict(numerical_rank=0,participation_rank=0.,entropy_rank=0.,d90=0,d95=0,leading_fraction=np.nan,stable_rank=0.),np.array([])
    s=np.linalg.svd(a,full_matrices=False,compute_uv=False)
    if not len(s) or s[0]<=0:return dict(numerical_rank=0,participation_rank=0.,entropy_rank=0.,d90=0,d95=0,leading_fraction=np.nan,stable_rank=0.),s
    tol=max(rel_tol*s[0],max(a.shape)*EPS*s[0]);nr=int((s>tol).sum());e=s*s;p=e/e.sum();nz=p>0
    return dict(numerical_rank=nr,participation_rank=float(1/np.sum(p*p)),entropy_rank=float(np.exp(-(p[nz]*np.log(p[nz])).sum())),
                d90=int(np.searchsorted(np.cumsum(p),.90)+1),d95=int(np.searchsorted(np.cumsum(p),.95)+1),leading_fraction=float(p[0]),stable_rank=float(e.sum()/s[0]**2)),s

def shannon(x):
    a=np.asarray(x,float);a=a[np.isfinite(a)&(a>0)]
    if not len(a):return np.nan
    p=a/a.sum();return float(-(p*np.log(p)).sum())

def prepare(df):
    req={'sample_id','taxon','flux','direction'};miss=req-set(df.columns)
    if miss:raise ValueError(f'Missing columns: {sorted(miss)}')
    if 'metabolite' not in df.columns:df=df.assign(metabolite=df['reaction'].astype(str))
    df=df[(df.direction.astype(str).str.lower()=='import')&(df.taxon.astype(str).str.lower()!='medium')].copy()
    df['flux_abs']=pd.to_numeric(df.flux,errors='coerce').abs();df['abundance']=pd.to_numeric(df.get('abundance'),errors='coerce')
    return df[np.isfinite(df.flux_abs)]

def matrix_for(sub,cutoff=1e-6,mode='flux'):
    x=sub[sub.flux_abs>cutoff].copy()
    if mode=='binary':x['value']=1.0
    elif mode=='weighted':x['value']=x.flux_abs*x.abundance.fillna(0)
    else:x['value']=x.flux_abs
    mat=x.pivot_table(index='taxon',columns='metabolite',values='value',aggfunc=('max' if mode=='binary' else 'sum'),fill_value=0.0)
    mat=mat.loc[(mat.abs().sum(1)>0),(mat.abs().sum(0)>0)]
    return mat

def run(exchanges_path,diversity_path,outdir):
    outdir=Path(outdir);outdir.mkdir(parents=True,exist_ok=True)
    df=prepare(pd.read_csv(exchanges_path));div=pd.read_csv(diversity_path).rename(columns={'supplier':'sample_id','shannon_diversity':'source_shannon'})
    rows=[];sing=[]
    for sid,sub in df.groupby('sample_id'):
        base={}
        for mode in ['binary','flux','weighted']:
            mat=matrix_for(sub,1e-6,mode);met,s=rank_metrics(mat.values)
            for k,v in enumerate(s,1):sing.append(dict(sample_id=sid,mode=mode,component=k,singular_value=v))
            base[mode]=(mat,met)
        ab=sub[['taxon','abundance']].drop_duplicates().abundance.dropna().to_numpy()
        row=dict(sample_id=sid,n_taxa=base['flux'][0].shape[0],n_resources=base['flux'][0].shape[1],shannon_from_abundance=shannon(ab))
        for mode in ['binary','flux','weighted']:
            for k,v in base[mode][1].items():row[f'{mode}_{k}']=v
        rows.append(row)
    dims=pd.DataFrame(rows).merge(div,on='sample_id',how='left').sort_values('sample_id')
    dims.to_csv(outdir/'table_S9_micom_dimensions.csv',index=False);pd.DataFrame(sing).to_csv(outdir/'micom_singular_values.csv',index=False)
    sens=[]
    for cutoff in [1e-6,1e-5,1e-4,1e-3]:
        for sid,sub in df.groupby('sample_id'):
            mat=matrix_for(sub,cutoff,'flux');met,_=rank_metrics(mat.values)
            sens.append(dict(sample_id=sid,flux_cutoff=cutoff,n_taxa=mat.shape[0],n_resources=mat.shape[1],participation_rank=met['participation_rank']))
    sens=pd.DataFrame(sens);sens.to_csv(outdir/'micom_cutoff_sensitivity.csv',index=False)
    z=dims.dropna(subset=['source_shannon']);rho,p=spearmanr(z.source_shannon,z.flux_participation_rank) if len(z)>=3 else (np.nan,np.nan)
    pd.DataFrame([dict(n=len(z),spearman_rho=rho,p_value=p,note='Descriptive only; n=4 in the manuscript analysis.')]).to_csv(outdir/'micom_diversity_correlation.csv',index=False)
    return dims,sens
