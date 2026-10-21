import numpy as np
from scipy.optimize import minimize

def softmax(a):
    z=a-a.max(1,keepdims=True); e=np.exp(z); return e/e.sum(1,keepdims=True)
def fit_latent_multinomial(P,d,seed=0,maxiter=500):
    P=np.asarray(P,float); P=P/P.sum(1,keepdims=True); n,s=P.shape; rng=np.random.default_rng(seed)
    # intercept + sample scores + taxon loadings; center scores implicitly by penalty
    b=np.log(np.maximum(P.mean(0),1e-12)); b-=b.mean(); Z=rng.normal(scale=.05,size=(n,d)); W=rng.normal(scale=.05,size=(s,d)); x=np.r_[b,Z.ravel(),W.ravel()]
    def fg(x):
        b=x[:s]; Z=x[s:s+n*d].reshape(n,d); W=x[s+n*d:].reshape(s,d); L=b+Z@W.T; Q=softmax(L)
        f=-np.sum(P*np.log(np.maximum(Q,1e-15)))+1e-5*(np.sum(Z*Z)+np.sum(W*W))
        D=Q-P; gb=D.sum(0); gZ=D@W+2e-5*Z; gW=D.T@Z+2e-5*W
        return f,np.r_[gb,gZ.ravel(),gW.ravel()]
    o=minimize(lambda q:fg(q),x,jac=True,method='L-BFGS-B',options={'maxiter':maxiter,'ftol':1e-10}); x=o.x; b=x[:s]; Z=x[s:s+n*d].reshape(n,d); W=x[s+n*d:].reshape(s,d); return softmax(b+Z@W.T)
def eta_dimension(P,max_d=5,seed=0):
    errs=[]
    for d in range(1,max_d+1):
        Q=fit_latent_multinomial(P,d,seed+d); kl=np.sum(P*np.log(np.maximum(P,1e-15)/np.maximum(Q,1e-15)),axis=1).mean(); errs.append(kl)
    slope=np.polyfit(np.arange(1,max_d+1),np.log(np.maximum(errs,1e-15)),1)[0]; return float(-1/slope),np.array(errs)
def metabolite_decay_dimension(X):
    X=np.asarray(X,float); nz=X[X>0]; off=np.percentile(nz,10) if len(nz) else 1.; Z=np.log(X+off); Z=(Z-Z.mean(0))/np.where(Z.std(0)>0,Z.std(0),1); s=np.linalg.svd(Z,compute_uv=False); d=np.arange(1,min(20,len(s))+1); slope=np.polyfit(d,np.log(np.maximum(s[:len(d)],1e-15)),1)[0]; return float(-1/slope)
