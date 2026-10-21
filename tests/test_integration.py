"""Computational checks aligned with the scientific methods."""
from pathlib import Path
import numpy as np
import pandas as pd
from microdim.cr_theory import interaction, factorized_interaction, qss
from microdim.metrics import participation_dimension, energy_rank, retained_improvement
from microdim.rrr import rrr_predict
from microdim.support_strength import ordinary_nmf, support_latent

ROOT=Path(__file__).resolve().parents[1]

def test_psd_and_factorization_rank():
    rng=np.random.default_rng(321)
    F=rng.uniform(.2,1,(7,3)); G=rng.uniform(.2,1,(5,2)); H=rng.uniform(.1,.9,(3,2))
    C=F@H@G.T; y=np.ones(5); a=np.ones(5)
    A=interaction(C,y,a)
    assert np.allclose(A,A.T,atol=1e-12)
    assert np.linalg.eigvalsh(A).min()>-1e-9
    assert np.linalg.matrix_rank(A,tol=1e-9)<=2
    assert np.allclose(A,factorized_interaction(F,G,H,y,a))


def test_participation_ratio_definition():
    assert np.isclose(participation_dimension(np.eye(4)),4)
    assert energy_rank(np.diag([3,1]),.9)==1
    assert np.isclose(retained_improvement(2,2,4),1)


def test_nmf_support_identity():
    C=np.array([[.2,0,.8],[0,.4,.3],[.1,.2,0.]])
    fixed,_,_=support_latent(C,k=2,seed=7,maxiter=400)
    assert np.array_equal(C>0,fixed>0)
    fit,_=ordinary_nmf(C,2,seed=7,max_iter=600)
    assert np.all(fit>=0) and fit.shape==C.shape


def test_rank_zero_baseline():
    X=np.array([[0.,1.],[1.,0.],[2.,1.],[3.,1.]])
    Y=np.array([[0.,2.],[1.,3.],[3.,2.],[4.,0.]])
    out=rrr_predict(X,Y,X,0)
    assert np.allclose(out,Y.mean(0))


def test_recorded_thresholds():
    ordinary=pd.read_csv(ROOT/'expected'/'hcom2_nmf_threshold.csv')
    fixed=pd.read_csv(ROOT/'expected'/'hcom2_support_strength.csv')
    assert ordinary.loc[(ordinary.retained_improvement>=.95)&(ordinary.p3_pearson>=.772953-.03),'rank'].min()==52
    assert fixed.loc[(fixed.retained_improvement>=.95)&(fixed.p3_pearson>=.772953-.03),'latent_dimension'].min()==6
    assert ordinary.loc[ordinary.seed_pass_fraction>.5,'rank'].min()==57
    assert fixed.loc[fixed.seed_pass_fraction>.5,'latent_dimension'].min()==6
    assert ordinary.loc[ordinary.seed_pass_fraction==1,'rank'].min()==60
    assert fixed.loc[fixed.seed_pass_fraction==1,'latent_dimension'].min()==12
