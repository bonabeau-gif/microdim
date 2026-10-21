import numpy as np
from scipy.integrate import solve_ivp

def interaction(C,y,a): return C@np.diag(y/a)@C.T
def factorized_interaction(F,G,H,y,a):
    # C = F H G^T
    C=F@H@G.T; return interaction(C,y,a)
def qss(rho,a,C,N): return (rho-C.T@N)/a
def rhs_full(t,z,C,y,m,rho,a,eps):
    S=C.shape[0]; N=z[:S]; R=z[S:]; dN=N*(C@(y*R)-m); dR=eps*R*(rho-a*R-C.T@N); return np.r_[dN,dR]
def rhs_red(t,N,C,y,m,rho,a):
    R=qss(rho,a,C,N); return N*(C@(y*R)-m)
