import numpy as np
from scipy.stats import spearmanr

def singular_values(x): return np.linalg.svd(np.asarray(x,float),full_matrices=False,compute_uv=False)
def participation_dimension(x):
    s=singular_values(x); e=s*s; den=np.sum(e*e); return float(np.sum(e)**2/den) if den>0 else np.nan
def entropy_dimension(x):
    s=singular_values(x); e=s*s; z=e.sum()
    if z<=0: return np.nan
    p=e[e>0]/z; return float(np.exp(-np.sum(p*np.log(p))))
def energy_rank(x, frac=.90):
    s=singular_values(x); e=s*s; z=e.sum()
    return int(np.searchsorted(np.cumsum(e)/z,frac)+1) if z>0 else 0
def leading_energy(x):
    s=singular_values(x); e=s*s; return float(e[0]/e.sum()) if e.sum()>0 else np.nan
def multivariate_r2(y_true,y_pred,y_train_mean):
    y=np.asarray(y_true,float); p=np.asarray(y_pred,float); mu=np.asarray(y_train_mean,float)
    num=np.sum((y-p)**2); den=np.sum((y-mu)**2); return float(1-num/den)
def log10_rmse(obs,pred,floor=1e-12):
    a=np.log10(np.maximum(np.asarray(obs,float),floor)); b=np.log10(np.maximum(np.asarray(pred,float),floor)); return float(np.sqrt(np.mean((a-b)**2)))
def retained_improvement(rmse_k, rmse_full, rmse_null): return float((rmse_null-rmse_k)/(rmse_null-rmse_full))
def clr(a, pseudocount):
    x=np.asarray(a,float).copy(); x[x<=0]=pseudocount; lx=np.log(x); return lx-lx.mean(axis=1,keepdims=True)
def support(x,tol=0): return (np.asarray(x)>tol).astype(int)
