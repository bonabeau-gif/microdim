import numpy as np
from microdim.support_controls import permute_species, degree_preserving_rewire

def test_species_shuffle_preserves_column_totals():
    x=np.array([[1,0,1,0],[0,1,1,0],[0,0,1,1],[1,1,0,0]])
    y=permute_species(x,seed=7)
    assert y.shape==x.shape
    assert np.array_equal(y.sum(axis=0),x.sum(axis=0))
    assert np.array_equal(np.sort(y.sum(axis=1)),np.sort(x.sum(axis=1)))

def test_rewire_preserves_both_degree_sequences():
    x=np.array([[1,0,1,0],[0,1,1,0],[0,0,1,1],[1,1,0,0]])
    y,n=degree_preserving_rewire(x,1000,seed=5)
    assert n>=0
    assert np.array_equal(y.sum(axis=0),x.sum(axis=0))
    assert np.array_equal(y.sum(axis=1),x.sum(axis=1))
