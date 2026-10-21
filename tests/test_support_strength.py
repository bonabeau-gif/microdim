import numpy as np
from microdim.support_strength import support_latent

def test_support_exact():
 C=np.array([[.2,0,.8],[0,.4,.3],[.1,.2,0.]])
 H,err,ok=support_latent(C,2,0,maxiter=500); assert np.array_equal(H>0,C>0); assert H.min()>=0 and H.max()<=1
