"""Connectivity scoring: four methods, Robust Rank Aggregation, benefit index.

Conventions: Z is the compound x gene consensus matrix (LINCS z-scores, BING genes).
A positive score always means "the drug moves genes the way exercise does".
"""
import numpy as np
from scipy.stats import beta, rankdata


# ---------- preprocessing ----------
def rank_positions(Z: np.ndarray, chunk: int = 2000) -> np.ndarray:
    """Position of each gene in each compound's list sorted high -> low (0 = most up)."""
    n, g = Z.shape
    R = np.empty((n, g), dtype=np.int16)
    ar = np.arange(g, dtype=np.int16)
    for s in range(0, n, chunk):
        order = np.argsort(-Z[s:s + chunk], axis=1)
        rows = np.arange(order.shape[0])[:, None]
        blk = np.empty_like(order, dtype=np.int16)
        blk[rows, order] = ar
        R[s:s + chunk] = blk
    return R


def drug_extremes(Z: np.ndarray, n: int, chunk: int = 4000):
    """Indices of each compound's top-n up and top-n down genes."""
    up = np.empty((Z.shape[0], n), dtype=np.int32)
    dn = np.empty((Z.shape[0], n), dtype=np.int32)
    for s in range(0, Z.shape[0], chunk):
        blk = Z[s:s + chunk]
        up[s:s + chunk] = np.argpartition(-blk, n, axis=1)[:, :n]
        dn[s:s + chunk] = np.argpartition(blk, n, axis=1)[:, :n]
    return up, dn


# ---------- method 1: CMap weighted connectivity score ----------
def enrichment_score(R: np.ndarray, Z: np.ndarray, S: np.ndarray) -> np.ndarray:
    """Weighted Kolmogorov-Smirnov ES of gene set S in every compound's ranked list.

    Computed analytically from hit positions (O(k log k) per compound) instead of
    walking all genes; weights are |z| of the hit genes (CMap weighted, p = 1).
    """
    g, k = R.shape[1], len(S)
    pos = R[:, S].astype(np.int32)
    w = np.abs(Z[:, S]).astype(np.float64)
    o = np.argsort(pos, axis=1)
    pos = np.take_along_axis(pos, o, 1)
    w = np.take_along_axis(w, o, 1)
    W = w.sum(1, keepdims=True)
    W[W == 0] = 1.0
    cw = np.cumsum(w, 1) / W
    miss = (pos - np.arange(k)[None, :]) / (g - k)
    after = cw - miss
    before = cw - w / W - miss
    mx, mn = after.max(1), before.min(1)
    return np.where(np.abs(mx) >= np.abs(mn), mx, mn)


def wtcs(R, Z, up, dn):
    es_u, es_d = enrichment_score(R, Z, up), enrichment_score(R, Z, dn)
    return np.where(np.sign(es_u) != np.sign(es_d), (es_u - es_d) / 2.0, 0.0)


# ---------- method 2: Spearman correlation on landmark genes ----------
def ranked_standardized(M: np.ndarray) -> np.ndarray:
    r = np.argsort(np.argsort(M, axis=-1), axis=-1).astype(np.float32)
    r -= r.mean(axis=-1, keepdims=True)
    r /= np.linalg.norm(r, axis=-1, keepdims=True) + 1e-12
    return r


def spearman(Zl_ranked: np.ndarray, t_landmark: np.ndarray) -> np.ndarray:
    return Zl_ranked @ ranked_standardized(t_landmark)


# ---------- method 3: XSum ----------
def xsum(t_full: np.ndarray, drug_up: np.ndarray, drug_dn: np.ndarray) -> np.ndarray:
    t0 = np.nan_to_num(t_full, nan=0.0)
    return t0[drug_up].sum(1) - t0[drug_dn].sum(1)


# ---------- method 4: GSFM embedding similarity ----------
def cosine_rows(A: np.ndarray, b: np.ndarray) -> np.ndarray:
    return (A @ b) / (np.linalg.norm(A, axis=1) * np.linalg.norm(b) + 1e-12)


def gsfm_score(D_up, D_dn, e_up, e_dn):
    return (cosine_rows(D_up, e_up) + cosine_rows(D_dn, e_dn)
            - cosine_rows(D_up, e_dn) - cosine_rows(D_dn, e_up))


# ---------- aggregation ----------
def norm_ranks(score: np.ndarray) -> np.ndarray:
    """1/n = best (highest score), 1 = worst."""
    return rankdata(-score, method="average") / len(score)


def rra(rank_matrix: np.ndarray) -> np.ndarray:
    """Robust Rank Aggregation rho score (Kolde et al. 2012), Bonferroni over methods."""
    r = np.sort(rank_matrix, axis=1)
    m = r.shape[1]
    k = np.arange(1, m + 1)
    p = beta.cdf(r, k[None, :], (m - k + 1)[None, :])
    return np.minimum(p.min(axis=1) * m, 1.0)


def final_rank(rho: np.ndarray, tiebreak: np.ndarray) -> np.ndarray:
    """1 = best. Order by RRA rho, ties (rho is capped at 1) broken by the tiebreak value
    (mean normalised method rank), so every compound gets a distinct position."""
    order = np.lexsort((tiebreak, rho))
    rank = np.empty(len(rho), dtype=np.int32)
    rank[order] = np.arange(1, len(rho) + 1)
    return rank


def projection(Z_sub: np.ndarray, t_sub: np.ndarray) -> np.ndarray:
    """Agreement of each drug with exercise on a gene subset (exercise-t weighted z-sum)."""
    return Z_sub @ t_sub / (np.linalg.norm(t_sub) + 1e-12)


def zscore(x: np.ndarray) -> np.ndarray:
    return (x - x.mean()) / (x.std() + 1e-12)
