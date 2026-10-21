import numpy as np
from scipy.optimize import minimize
from sklearn.decomposition import NMF
from scipy.special import expit

def ordinary_nmf(C,k,seed,max_iter=5000):
    m=NMF(n_components=k,init='random',random_state=seed,max_iter=max_iter,tol=1e-4,solver='cd',beta_loss='frobenius'); W=m.fit_transform(C); H=m.components_; return W@H, float(m.reconstruction_err_)
def support_latent(C,k,seed=0,l2=1e-6,maxiter=2000):
    C=np.asarray(C,float); S=C>0; ii,jj=np.where(S); y=C[ii,jj]; n,m=C.shape; rng=np.random.default_rng(seed)
    a=np.zeros(n); b=np.zeros(m); U=rng.normal(scale=.01,size=(n,k)); V=rng.normal(scale=.01,size=(m,k)); x=np.r_[a,b,U.ravel(),V.ravel()]
    def unpack(x):
        a=x[:n]; b=x[n:n+m]; U=x[n+m:n+m+n*k].reshape(n,k); V=x[n+m+n*k:].reshape(m,k); return a,b,U,V
    def fg(x):
        a,b,U,V=unpack(x); z=a[ii]+b[jj]+(U[ii]*V[jj]).sum(1) if k else a[ii]+b[jj]; p=expit(z); d=p-y; dz=2*d*p*(1-p)/len(y); f=np.mean(d*d)+l2*np.sum(x*x)
        ga=np.bincount(ii,dz,minlength=n); gb=np.bincount(jj,dz,minlength=m); gU=np.zeros_like(U); gV=np.zeros_like(V)
        if k:
            np.add.at(gU,ii,dz[:,None]*V[jj]); np.add.at(gV,jj,dz[:,None]*U[ii])
        g=np.r_[ga,gb,gU.ravel(),gV.ravel()]+2*l2*x; return f,g
    o=minimize(lambda z:fg(z),x,jac=True,method='L-BFGS-B',options={'maxiter':maxiter,'ftol':1e-12}); a,b,U,V=unpack(o.x); Z=a[:,None]+b[None,:]+(U@V.T if k else 0); Chat=S*expit(Z); return Chat,float(np.mean((Chat[S]-C[S])**2)),bool(o.success)
