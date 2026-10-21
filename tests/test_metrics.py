import numpy as np
from microdim.metrics import participation_dimension,entropy_dimension,energy_rank,multivariate_r2

def test_participation_rank1(): assert abs(participation_dimension(np.array([[1.,2.],[2.,4.]]))-1)<1e-12
def test_identity_dimension(): assert abs(participation_dimension(np.eye(4))-4)<1e-12
def test_energy_rank(): assert energy_rank(np.diag([3.,1.]),.9)==1
def test_r2_perfect():
 y=np.array([[1.,2.],[2.,3.]]); assert abs(multivariate_r2(y,y,np.zeros(2))-1)<1e-12
