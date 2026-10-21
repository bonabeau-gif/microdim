"""Null models of binary consumer-by-resource support.

These routines generate counterfactual supports only. Predictive evaluation and
model-parameter fitting are performed by the relevant ecological model kernel.
"""
import numpy as np


def _binary(mask):
    s = np.asarray(mask)
    if s.ndim != 2 or not np.isin(s, [0, 1]).all():
        raise ValueError('Expected a two-dimensional binary support matrix')
    return s.astype(np.uint8, copy=True)


def permute_species(mask, seed=0):
    """Permute row identities while retaining every resource sharing relationship."""
    s = _binary(mask)
    return s[np.random.default_rng(seed).permutation(s.shape[0]), :]


def degree_preserving_rewire(mask, attempts=1000, seed=0):
    """Double-edge swaps conserving the degree of every consumer and resource.

    A swap (i,j),(k,l) -> (i,l),(k,j) is accepted only for i!=k,
    j!=l, and the proposed two new links are absent. This generates a
    support-only null, not an automatically viable resource-kinetic fit.
    """
    s = _binary(mask)
    if attempts < 0:
        raise ValueError('attempts must be nonnegative')
    rng = np.random.default_rng(seed)
    accepted = 0
    for _ in range(attempts):
        ii, jj = np.nonzero(s)
        if len(ii) < 2:
            break
        a, b = rng.choice(len(ii), 2, replace=False)
        i, j, k, l = int(ii[a]), int(jj[a]), int(ii[b]), int(jj[b])
        if i == k or j == l or s[i, l] or s[k, j]:
            continue
        s[i, j] = 0
        s[k, l] = 0
        s[i, l] = 1
        s[k, j] = 1
        accepted += 1
    return s, accepted
