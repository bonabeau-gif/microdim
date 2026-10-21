import numpy as np
from microdim.rrr import rrr_predict

def test_rrr_shapes():
 rng=np.random.default_rng(1); X=rng.normal(size=(30,5)); B=rng.normal(size=(5,4)); Y=X@B
 for r in [0,1,3]: assert rrr_predict(X,Y,X[:3],r).shape==(3,4)
