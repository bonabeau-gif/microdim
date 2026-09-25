from __future__ import annotations
import numpy as np


def _safe_half_min_positive(a, default=1e-12):
    x = np.asarray(a, float)
    vals = x[np.isfinite(x) & (x > 0)]
    return float(0.5 * vals.min()) if vals.size else default


def clr_transform(a, pseudocount):
    x = np.asarray(a, float).copy()
    x[~np.isfinite(x)] = 0.0
    x[x <= 0] = pseudocount
    lx = np.log(x)
    return lx - lx.mean(axis=1, keepdims=True)


class FoldPreprocessor:
    """Leakage-safe preprocessing fitted on one outer training fold.

    The manuscript specifies the prevalence/observation filters, feature caps,
    CLR/log transforms, imputation, and training-only standardization. The
    historical exact pseudocount convention was not preserved in the manuscript;
    this implementation uses one half of the minimum positive training value and
    exposes that convention here rather than hiding it.
    """
    def __init__(self, taxa_presence=0.25, metabolite_observation=0.75, cap_fraction=0.25):
        self.taxa_presence = taxa_presence
        self.metabolite_observation = metabolite_observation
        self.cap_fraction = cap_fraction

    def fit(self, Xraw, Yraw):
        Xraw = np.asarray(Xraw, float)
        Yraw = np.asarray(Yraw, float)
        n = Xraw.shape[0]
        cap = max(1, int(np.floor(self.cap_fraction * n)))

        present = np.mean(np.isfinite(Xraw) & (Xraw > 0), axis=0) >= self.taxa_presence
        idx = np.flatnonzero(present)
        if idx.size:
            means = np.nanmean(np.where(Xraw[:, idx] > 0, Xraw[:, idx], np.nan), axis=0)
            means = np.nan_to_num(means, nan=-np.inf)
            idx = idx[np.argsort(means)[::-1][:cap]]
        self.x_idx = idx
        self.x_pc = _safe_half_min_positive(Xraw[:, self.x_idx]) if len(self.x_idx) else 1e-12
        X = clr_transform(Xraw[:, self.x_idx], self.x_pc) if len(self.x_idx) else np.empty((n,0))
        self.x_mean = X.mean(0) if X.shape[1] else np.array([])
        self.x_sd = X.std(0, ddof=0) if X.shape[1] else np.array([])
        self.x_sd[self.x_sd == 0] = 1.0

        observed = np.mean(np.isfinite(Yraw), axis=0) >= self.metabolite_observation
        yidx = np.flatnonzero(observed)
        self.y_pc = _safe_half_min_positive(Yraw[:, yidx]) if len(yidx) else 1e-12
        # Median-impute on the original scale, then log; rank by training log variance.
        med = np.nanmedian(Yraw[:, yidx], axis=0) if len(yidx) else np.array([])
        self.y_median_all = med
        Yimp = Yraw[:, yidx].copy() if len(yidx) else np.empty((n,0))
        if len(yidx):
            inds = np.where(~np.isfinite(Yimp))
            Yimp[inds] = med[inds[1]]
            Yimp[Yimp <= 0] = self.y_pc
            Ylog = np.log(Yimp)
            variances = np.var(Ylog, axis=0)
            keep_local = np.argsort(variances)[::-1][:cap]
            self.y_idx = yidx[keep_local]
            self.y_median = med[keep_local]
            Ylog = Ylog[:, keep_local]
        else:
            self.y_idx = yidx
            self.y_median = np.array([])
            Ylog = np.empty((n,0))
        self.y_mean = Ylog.mean(0) if Ylog.shape[1] else np.array([])
        self.y_sd = Ylog.std(0, ddof=0) if Ylog.shape[1] else np.array([])
        self.y_sd[self.y_sd == 0] = 1.0
        return self

    def transform_X(self, Xraw):
        X = clr_transform(np.asarray(Xraw, float)[:, self.x_idx], self.x_pc)
        return (X - self.x_mean) / self.x_sd

    def transform_Y(self, Yraw):
        Y = np.asarray(Yraw, float)[:, self.y_idx].copy()
        inds = np.where(~np.isfinite(Y))
        if len(inds[0]): Y[inds] = self.y_median[inds[1]]
        Y[Y <= 0] = self.y_pc
        Y = np.log(Y)
        return (Y - self.y_mean) / self.y_sd
