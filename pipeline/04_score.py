"""Step 4: score every LINCS compound against the 18 human exercise signatures.

Methods (positive = drug mimics exercise):
  1. WTCS     CMap weighted connectivity score (weighted KS on query up/down sets)
  2. Spearman rank correlation with the exercise t-statistics, landmark genes
  3. XSum     exercise t summed over each drug's top-100 up minus top-100 down genes
  4. GSFM     cosine similarity of gene-set foundation-model embeddings
                (exercise up/down sets vs drug top-100 up/down sets)
Each score is standardised against the compound's own null from N_PERM random queries
(gene-label permutations of the exercise t-vector), as CMap normalises connectivity per
perturbagen, so compounds with large global transcriptional effects do not dominate.
The standardised scores are combined with Robust Rank Aggregation (RRA); significance comes
from the pooled null of standardised RRA scores, BH-FDR.
Benefit index: agreement on adaptive-module genes minus agreement on stress-module genes.
"""
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from tejas import config as C  # noqa: E402
from tejas import scoring as S  # noqa: E402
from tejas.gsfm_embed import Embedder  # noqa: E402

METHODS = ["wtcs", "spearman", "xsum", "gsfm"]


def cached(path: Path, fn):
    if path.exists():
        return np.load(path)
    arr = fn()
    np.save(path, arr)
    return arr


def bh(p: np.ndarray) -> np.ndarray:
    n = len(p)
    o = np.argsort(p)
    q = p[o] * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.minimum(q, 1.0)
    return out


def query_sets(sub: pd.DataFrame):
    """Top FDR-significant genes per direction; fall back to top-|t| genes if too few."""
    up = sub[(sub.fdr < C.FDR) & (sub.t > 0)].nlargest(C.N_QUERY, "t")
    dn = sub[(sub.fdr < C.FDR) & (sub.t < 0)].nsmallest(C.N_QUERY, "t")
    fallback = len(up) < C.MIN_QUERY or len(dn) < C.MIN_QUERY
    if fallback:
        up = sub[sub.t > 0].nlargest(50, "t")
        dn = sub[sub.t < 0].nsmallest(50, "t")
    return up, dn, fallback


class Scorer:
    def __init__(self):
        t0 = time.time()
        self.genes = pd.read_parquet(C.PROC / "lincs_genes.parquet")
        self.cpd = pd.read_parquet(C.PROC / "compounds.parquet")
        self.Z = np.load(C.PROC / "lincs_consensus_z.npy")
        self.G = self.Z.shape[1]
        self.idx_of = dict(zip(self.genes.gene_id.astype(str), self.genes.idx))
        self.symbols = self.genes.gene_symbol.astype(str).to_numpy()
        self.landmark = (self.genes.feature_space == "landmark").to_numpy()
        cls = pd.read_parquet(C.PROC / "gene_class.parquet")
        gclass = dict(zip(cls.gene_symbol, cls.gene_class))
        self.gene_class = np.array([gclass.get(s, "") for s in self.symbols])

        print("ranking genes within each compound ...", flush=True)
        self.R = cached(C.PROC / "lincs_rank_positions.npy", lambda: S.rank_positions(self.Z))
        up_p, dn_p = C.PROC / "lincs_drug_up100.npy", C.PROC / "lincs_drug_dn100.npy"
        if not (up_p.exists() and dn_p.exists()):
            u, d = S.drug_extremes(self.Z, C.N_DRUG_SET)
            np.save(up_p, u)
            np.save(dn_p, d)
        self.d_up, self.d_dn = np.load(up_p), np.load(dn_p)

        print("loading GSFM and embedding drug gene sets ...", flush=True)
        self.emb = Embedder(C.GSFM_MODEL)

        def embed_drugs(idx):
            return self.emb.encode([list(self.symbols[r]) for r in idx])
        D_up = cached(C.PROC / "gsfm_drug_up.npy", lambda: embed_drugs(self.d_up))
        D_dn = cached(C.PROC / "gsfm_drug_dn.npy", lambda: embed_drugs(self.d_dn))
        self.center = np.vstack([D_up, D_dn]).mean(axis=0)
        self.D_up, self.D_dn = D_up - self.center, D_dn - self.center
        self.gsfm_coverage = self.emb.coverage(self.symbols)
        self._zl = {}
        print(f"ready in {(time.time() - t0) / 60:.1f} min  (GSFM vocab covers "
              f"{self.gsfm_coverage:.1%} of L1000 genes)", flush=True)

    def zl_ranked(self, tissue: str, lm_idx: np.ndarray):
        if tissue not in self._zl:
            self._zl[tissue] = S.ranked_standardized(self.Z[:, lm_idx])
        return self._zl[tissue]

    def embed_sets(self, sets):
        return self.emb.encode(sets) - self.center

    def score(self, tissue, t, up_idx, dn_idx, lm_idx, e_up, e_dn):
        return {
            "wtcs": S.wtcs(self.R, self.Z, up_idx, dn_idx),
            "spearman": S.spearman(self.zl_ranked(tissue, lm_idx), t[lm_idx]),
            "xsum": S.xsum(t, self.d_up, self.d_dn),
            "gsfm": S.gsfm_score(self.D_up, self.D_dn, e_up, e_dn),
        }


def combine(scores: dict, null_stats: dict | None = None):
    """Rank each method's score (optionally standardised per compound) and aggregate with RRA."""
    if null_stats is not None:
        scores = {m: (scores[m] - null_stats[m][0]) / null_stats[m][1] for m in METHODS}
    ranks = np.column_stack([S.norm_ranks(scores[m]) for m in METHODS])
    rho = S.rra(ranks)
    inv_rho = S.rra(1.0 - ranks + 1.0 / len(ranks))
    return ranks, rho, inv_rho


def random_queries(sc, t, n_up, n_dn, rng, n_perm):
    """Gene-label permutations of the exercise t-vector and the query sets they imply."""
    measured = np.flatnonzero(~np.isnan(t))
    perms, sets = [], []
    for _ in range(n_perm):
        tp = t.copy()
        tp[measured] = rng.permutation(t[measured])
        o = measured[np.argsort(-tp[measured])]
        pu, pdn = o[:n_up], o[::-1][:n_dn]
        perms.append((tp, pu, pdn))
        sets += [list(sc.symbols[pu]), list(sc.symbols[pdn])]
    return perms, sc.embed_sets(sets)


def score_signature(sc, tissue, t, up_idx, dn_idx, rng, n_perm=C.N_PERM):
    """Score all compounds, standardise each compound against its own random-query null
    (CMap-style normalisation, so compounds with large global effects do not dominate),
    aggregate with RRA, and return pooled-null permutation p-values for the combined score."""
    lm_idx = np.flatnonzero(sc.landmark & ~np.isnan(t))
    e_up, e_dn = sc.embed_sets([list(sc.symbols[up_idx]), list(sc.symbols[dn_idx])])
    raw = sc.score(tissue, t, up_idx, dn_idx, lm_idx, e_up, e_dn)
    perms, pe = random_queries(sc, t, len(up_idx), len(dn_idx), rng, n_perm)
    null_scores = {m: np.empty((n_perm, sc.Z.shape[0]), dtype=np.float32) for m in METHODS}
    for i, (tp, pu, pdn) in enumerate(perms):
        s = sc.score(tissue, tp, pu, pdn, lm_idx, pe[2 * i], pe[2 * i + 1])
        for m in METHODS:
            null_scores[m][i] = s[m]
    stats = {m: (null_scores[m].mean(0), null_scores[m].std(0) + 1e-6) for m in METHODS}
    ranks, rho, inv_rho = combine(raw, stats)
    null_rho = np.empty((n_perm, sc.Z.shape[0]))
    for i in range(n_perm):
        _, null_rho[i], _ = combine({m: null_scores[m][i] for m in METHODS}, stats)
    null_sorted = np.sort(null_rho.ravel())
    perm_p = (np.searchsorted(null_sorted, rho, side="right") + 1) / (null_sorted.size + 1)
    _, raw_rho, _ = combine(raw)
    z = {m: (raw[m] - stats[m][0]) / stats[m][1] for m in METHODS}
    return {"raw": raw, "z": z, "ranks": ranks, "rho": rho, "inv_rho": inv_rho, "perm_p": perm_p,
            "raw_rho": raw_rho, "stats": stats}


def main():
    t_start = time.time()
    sc = Scorer()
    sig = pd.read_parquet(C.PROC / "signatures_long.parquet")
    sig["entrez"] = sig["entrez"].astype(str)
    sig = sig[sig.entrez.isin(sc.idx_of)].copy()
    sig["gidx"] = sig.entrez.map(sc.idx_of).astype(int)
    rng = np.random.default_rng(C.SEED)

    out, queries = [], {}
    for sid, sub in sig.groupby("signature", sort=False):
        t1 = time.time()
        tissue, modality, timepoint = sid.split("|")
        t = np.full(sc.G, np.nan)
        t[sub.gidx.to_numpy()] = sub.t.to_numpy()
        measured = np.flatnonzero(~np.isnan(t))
        lm_idx = np.flatnonzero(sc.landmark & ~np.isnan(t))

        up, dn, fallback = query_sets(sub)
        up_idx, dn_idx = up.gidx.to_numpy(), dn.gidx.to_numpy()
        r_ = score_signature(sc, tissue, t, up_idx, dn_idx, rng)
        scores, ranks, rho, inv_rho, perm_p = r_["raw"], r_["ranks"], r_["rho"], r_["inv_rho"], r_["perm_p"]
        np.savez_compressed(C.PROC / f"null_stats_{sid.replace('|', '_')}.npz",
                            **{f"{m}_{k}": r_["stats"][m][j] for m in METHODS for j, k in enumerate(("mu", "sd"))})

        # benefit index: adaptive vs stress agreement
        A = measured[sc.gene_class[measured] == "adaptive"]
        St = measured[sc.gene_class[measured] == "stress"]
        a_raw = S.projection(sc.Z[:, A], t[A])
        s_raw = S.projection(sc.Z[:, St], t[St])
        a_z, s_z = S.zscore(a_raw), S.zscore(s_raw)

        df = pd.DataFrame({"signature": sid, "tissue": tissue, "modality": modality, "timepoint": timepoint,
                           "drug_id": sc.cpd.drug_id.to_numpy()})
        for j, m in enumerate(METHODS):
            df[m] = scores[m].astype(np.float32)
            df[f"{m}_z"] = r_["z"][m].astype(np.float32)
            df[f"r_{m}"] = ranks[:, j].astype(np.float32)
        df["rank_unnormalised"] = S.rankdata(r_["raw_rho"], method="min").astype(np.int32)
        df["rho"] = rho
        df["inv_rho"] = inv_rho
        df["perm_p"] = perm_p
        df["perm_fdr"] = bh(perm_p)
        df["adaptive_z"] = a_z.astype(np.float32)
        df["stress_z"] = s_z.astype(np.float32)
        df["benefit"] = (a_z - s_z).astype(np.float32)
        out.append(df)

        queries[sid] = {"up": list(sc.symbols[up_idx]), "down": list(sc.symbols[dn_idx]),
                        "fallback": bool(fallback), "n_measured": int(len(measured)),
                        "n_adaptive_genes": int(len(A)), "n_stress_genes": int(len(St))}
        n_sig = int((df.perm_fdr < 0.05).sum())
        print(f"{sid:32s} up {len(up_idx):3d} dn {len(dn_idx):3d}{' (fallback)' if fallback else ''}  "
              f"FDR<0.05: {n_sig:5d}  top: {', '.join(sc.cpd.name.to_numpy()[np.argsort(rho)[:3]])}  "
              f"[{time.time() - t1:.0f}s]", flush=True)

    res = pd.concat(out, ignore_index=True)
    (C.RESULTS / "query_sets.json").write_text(json.dumps(queries, indent=1))
    informative = [s for s, q in queries.items() if not q["fallback"]]
    finalize(res, informative)
    meta = {"gsfm_vocab_coverage_of_l1000_genes": sc.gsfm_coverage, "n_compounds": int(len(sc.cpd)),
            "informative_signatures": informative,
            "n_perm": C.N_PERM, "minutes": round((time.time() - t_start) / 60, 1)}
    (C.RESULTS / "scoring_meta.json").write_text(json.dumps(meta, indent=1))
    print(f"done in {meta['minutes']} min")


def finalize(res: pd.DataFrame, informative: list):
    """Distinct ranks per signature (rho, ties broken by mean method rank) and the overall
    ranking per exercise type: RRA across informative signatures only (signatures that needed
    the top-|t| fallback have too few FDR-significant genes to trust)."""
    rcols = [f"r_{m}" for m in METHODS]
    parts = []
    for _, g in res.groupby("signature", sort=False):
        g = g.copy()
        mean_r = g[rcols].to_numpy().mean(axis=1)
        g["rank"] = S.final_rank(g.rho.to_numpy(), mean_r)
        g["inv_rank"] = S.final_rank(g.inv_rho.to_numpy(), -mean_r)
        parts.append(g)
    res = pd.concat(parts, ignore_index=True)
    res.to_parquet(C.RESULTS / "scores_long.parquet", index=False)

    overall = []
    for modality, g in res[res.signature.isin(informative)].groupby("modality"):
        wide = g.pivot(index="drug_id", columns="signature", values="rank")
        nr = wide.to_numpy() / len(wide)
        med = np.median(wide.to_numpy(), axis=1)
        o = pd.DataFrame({"modality": modality, "drug_id": wide.index, "rho": S.rra(nr), "median_rank": med})
        o["rank"] = S.final_rank(o.rho.to_numpy(), med)
        overall.append(o)
    pd.concat(overall).to_parquet(C.RESULTS / "overall.parquet", index=False)


if __name__ == "__main__":
    if "--rerank" in sys.argv:
        meta = json.loads((C.RESULTS / "scoring_meta.json").read_text())
        finalize(pd.read_parquet(C.RESULTS / "scores_long.parquet"), meta["informative_signatures"])
        print("ranks recomputed from saved scores")
    else:
        main()
