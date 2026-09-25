import numpy as np
from microdim_repro.rank_calibration import make_B, cv_select
from microdim_repro.micom_proxy import rank_metrics
from microdim_repro.consumer_resource import generate_system, reduced_parameters


def test_make_B_rank_and_population_scale():
    rng=np.random.default_rng(1);B=make_B(rng,24,24,10,.2)
    assert np.linalg.matrix_rank(B,tol=1e-10)==10
    # total signal variance per outcome = R2/(1-R2)
    assert np.isclose(np.sum(B*B)/24,.2/.8,rtol=1e-10)


def test_participation_rank_equal_spectrum():
    A=np.diag([2.,2.,2.])
    m,_=rank_metrics(A)
    assert np.isclose(m['participation_rank'],3.0)
    assert m['d90']==3


def test_consumer_resource_factorization_and_rank():
    s=generate_system(3);growth,A,M,Af=reduced_parameters(s)
    assert np.linalg.norm(A-Af)/np.linalg.norm(A)<1e-12
    assert np.linalg.matrix_rank(A,tol=1e-10)<=2
    assert np.linalg.eigvalsh(A).max()<1e-12
    assert np.all(A<=1e-12)


def test_rank_selection_is_deterministic():
    rng=np.random.default_rng(99);X=rng.normal(size=(60,8));B=make_B(rng,8,8,2,.3);Y=X@B+rng.normal(size=(60,8))
    a=cv_select(X,Y,99);b=cv_select(X,Y,99)
    assert a==b
