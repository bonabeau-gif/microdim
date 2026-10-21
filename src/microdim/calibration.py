import numpy as np
from sklearn.model_selection import KFold
from .rrr import rrr_predict
from .metrics import multivariate_r2

def orth(rng,n,r): return np.linalg.qr(rng.normal(size=(n,r)))[0][:,:r]
def simulate(n,p=24,q=24,r=10,r2=.2,seed=0):
    rng=np.random.default_rng(seed); U=orth(rng,p,r); V=orth(rng,q,r); s=np.sqrt(q*r2/(r*(1-r2))); B=U@(np.eye(r)*s)@V.T; X=rng.normal(size=(n,p)); Y=X@B+rng.normal(size=(n,q)); return X,Y,B,V
def select_rank(X,Y,alpha=10,rmax=12,seed=0):
    vals=[]; cv=KFold(5,shuffle=True,random_state=seed)
    for rank in range(rmax+1):
        sc=[]
        for tr,va in cv.split(X): sc.append(multivariate_r2(Y[va],rrr_predict(X[tr],Y[tr],X[va],rank,alpha),Y[tr].mean(0)))
        vals.append(np.mean(sc))
    return int(np.argmax(vals))
def grid(ns,r2s,reps=100,seed=20260925):
    out=[]
    for n in ns:
      for r2 in r2s:
        rr=[]
        for rep in range(reps):
          X,Y,_,_=simulate(n,r2=r2,seed=seed+100000*n+1000*int(r2*100)+rep); rr.append(select_rank(X,Y,seed=seed+rep))
        out.append((n,r2,float(np.median(rr)),float(np.mean(rr))))
    return out
