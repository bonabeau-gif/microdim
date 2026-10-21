import numpy as np
from sklearn.linear_model import Ridge
from sklearn.decomposition import PCA

def ridge_coef(X,Y,alpha=10.0):
    m=Ridge(alpha=alpha,fit_intercept=False).fit(X,Y); return m.coef_.T

def rrr_predict(Xtr,Ytr,Xte,rank,alpha=10.0):
    mu=Ytr.mean(0)
    if rank==0: return np.tile(mu,(len(Xte),1))
    B=ridge_coef(Xtr,Ytr,alpha); F=Xtr@B
    _,_,Vt=np.linalg.svd(F,full_matrices=False); V=Vt[:rank].T
    Br=B@V@V.T; return Xte@Br

def pcr_predict(Xtr,Ytr,Xte,ncomp):
    p=PCA(n_components=ncomp,svd_solver='full').fit(Xtr); Z=p.transform(Xtr); zt=p.transform(Xte)
    B=np.linalg.lstsq(np.c_[np.ones(len(Z)),Z],Ytr,rcond=None)[0]
    return np.c_[np.ones(len(zt)),zt]@B
