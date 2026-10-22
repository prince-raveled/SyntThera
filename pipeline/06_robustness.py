"""Step 6: robustness of the ranking and a map of drug space.

* Method ablation: re-aggregate with subsets of the four methods (no rescoring needed) and measure
  positive-control recovery and agreement with the full ensemble. Shows what GSFM adds.
* Bootstrap stability: 30 runs per informative signature in which 20% of the measured genes have
  their exercise t-values randomly reassigned; all four methods are rescored and re-aggregated.
  Reports Jaccard overlap of the top-50 approved drugs with the unperturbed ranking.
* Drug map: UMAP (cosine) of every compound's GSFM up/down embedding, with the exercise
  signatures projected into the same space.
"""
import importlib.util
import json
import sys
import time
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from tejas import config as C  # noqa: E402
from tejas import scoring as S  # noqa: E402

spec = importlib.util.spec_from_file_location("score04", Path(__file__).with_name("04_score.py"))
score04 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(score04)

CONFIGS = {
    "Full ensemble (4 methods)": ["wtcs", "spearman", "xsum", "gsfm"],
    "Without GSFM": ["wtcs", "spearman", "xsum"],
    "GSFM only": ["gsfm"],
    "Without weighted connectivity": ["spearman", "xsum", "gsfm"],
    "Without XSum": ["wtcs", "spearman", "gsfm"],
}
N_BOOT, NOISE, TOPK = 30, 0.20, 50


def jaccard(a, b):
    a, b = set(a), set(b)
    return len(a & b) / max(len(a | b), 1)


def main():
    t0 = time.time()
    res = pd.read_parquet(C.RESULTS / "scores_long.parquet")
    cpd = pd.read_parquet(C.PROC / "compounds.parquet")
    meta = json.loads((C.RESULTS / "scoring_meta.json").read_text())
    informative = meta["informative_signatures"]
    n = len(cpd)
    approved = (cpd.clinical_phase == "Launched").to_numpy()
    drug_ids = cpd.drug_id.to_numpy()

    # ---------- ablation ----------
    abl_rows = []
    for sid in informative:
        rs = res[res.signature == sid].set_index("drug_id").loc[drug_ids]
        mod = sid.split("|")[1]
        controls = [c for c in C.POSITIVE_CONTROLS[mod] if c in set(drug_ids)]
        full_rank = None
        for name, methods in CONFIGS.items():
            ranks = rs[[f"r_{m}" for m in methods]].to_numpy()
            rho = S.rra(ranks) if len(methods) > 1 else ranks[:, 0]
            rank = S.final_rank(rho, ranks.mean(axis=1))
            if full_rank is None:
                full_rank = rank
            top_full = drug_ids[approved][np.argsort(full_rank[approved])[:TOPK]]
            top_cfg = drug_ids[approved][np.argsort(rank[approved])[:TOPK]]
            ctrl_pct = np.mean([rank[drug_ids == c][0] / n * 100 for c in controls])
            abl_rows.append({"signature": sid, "config": name, "control_mean_percentile": ctrl_pct,
                             "spearman_vs_full": spearmanr(rank, full_rank).statistic,
                             "top50_jaccard_vs_full": jaccard(top_full, top_cfg)})
    abl = pd.DataFrame(abl_rows)
    abl_sum = abl.groupby("config", sort=False)[["control_mean_percentile", "spearman_vs_full",
                                                 "top50_jaccard_vs_full"]].mean().round(3).reset_index()
    print("method ablation (mean over informative signatures):")
    print(abl_sum.to_string(index=False))

    # ---------- bootstrap stability ----------
    sc = score04.Scorer()
    sig = pd.read_parquet(C.PROC / "signatures_long.parquet")
    sig["entrez"] = sig["entrez"].astype(str)
    sig = sig[sig.entrez.isin(sc.idx_of)].copy()
    sig["gidx"] = sig.entrez.map(sc.idx_of).astype(int)
    rng = np.random.default_rng(C.SEED)
    boot_rows = []
    for sid in informative:
        sub = sig[sig.signature == sid]
        tissue = sid.split("|")[0]
        t = np.full(sc.G, np.nan)
        t[sub.gidx.to_numpy()] = sub.t.to_numpy()
        measured = np.flatnonzero(~np.isnan(t))
        lm_idx = np.flatnonzero(sc.landmark & ~np.isnan(t))
        up, dn, _ = score04.query_sets(sub)
        k_up, k_dn = len(up), len(dn)
        npz = np.load(C.PROC / f"null_stats_{sid.replace('|', '_')}.npz")
        stats = {m: (npz[f"{m}_mu"], npz[f"{m}_sd"]) for m in score04.METHODS}
        base = res[res.signature == sid].set_index("drug_id").loc[drug_ids]["rank"].to_numpy()
        top_base = set(drug_ids[approved][np.argsort(base[approved])[:TOPK]])
        tops, sets = [], []
        for _ in range(N_BOOT):
            tp = t.copy()
            pick = rng.choice(measured, int(NOISE * len(measured)), replace=False)
            tp[pick] = rng.permutation(tp[pick])
            o = measured[np.argsort(-tp[measured])]
            pu, pdn = o[:k_up], o[::-1][:k_dn]
            sets.append((tp, pu, pdn))
        emb = sc.embed_sets([s for _, pu, pdn in sets for s in (list(sc.symbols[pu]), list(sc.symbols[pdn]))])
        for b, (tp, pu, pdn) in enumerate(sets):
            rks, rho, _ = score04.combine(sc.score(tissue, tp, pu, pdn, lm_idx, emb[2 * b], emb[2 * b + 1]), stats)
            rank = S.final_rank(rho, rks.mean(axis=1))
            tops.append(set(drug_ids[approved][np.argsort(rank[approved])[:TOPK]]))
        vs_base = [jaccard(tb, top_base) for tb in tops]
        pairwise = [jaccard(a, b) for a, b in combinations(tops, 2)]
        shared = TOPK * TOPK / approved.sum()  # expected overlap of two random top-K lists
        boot_rows.append({"signature": sid, "jaccard_vs_unperturbed_median": float(np.median(vs_base)),
                          "jaccard_pairwise_median": float(np.median(pairwise)),
                          "jaccard_random_expected": shared / (2 * TOPK - shared)})
        print(f"bootstrap {sid:30s} top-{TOPK} Jaccard vs unperturbed {np.median(vs_base):.2f}, "
              f"pairwise {np.median(pairwise):.2f}", flush=True)
    boot = pd.DataFrame(boot_rows)

    # ---------- drug map ----------
    X = np.hstack([sc.D_up, sc.D_dn]).astype(np.float32)
    sig_emb, sig_ids = [], []
    for sid in sorted(sig.signature.unique()):
        up, dn, _ = score04.query_sets(sig[sig.signature == sid])
        e_up, e_dn = sc.embed_sets([list(sc.symbols[up.gidx.to_numpy()]), list(sc.symbols[dn.gidx.to_numpy()])])
        sig_emb.append(np.hstack([e_up, e_dn]))
        sig_ids.append(sid)
    try:
        import umap
        reducer = umap.UMAP(n_neighbors=30, min_dist=0.15, metric="cosine", random_state=C.SEED, low_memory=True)
        xy = reducer.fit_transform(X)
        sxy = reducer.transform(np.vstack(sig_emb).astype(np.float32))
        method = "UMAP (cosine, n_neighbors=30)"
    except ImportError:
        Xc = X - X.mean(0)
        _, _, vt = np.linalg.svd(Xc[np.random.default_rng(0).choice(len(Xc), 5000, replace=False)], full_matrices=False)
        xy, sxy = Xc @ vt[:2].T, (np.vstack(sig_emb) - X.mean(0)) @ vt[:2].T
        method = "PCA"
    pd.DataFrame({"drug_id": drug_ids, "x": xy[:, 0], "y": xy[:, 1]}).to_parquet(C.RESULTS / "drug_map.parquet", index=False)

    out = {
        "ablation_summary": abl_sum.to_dict(orient="records"),
        "ablation": abl.round(4).to_dict(orient="records"),
        "bootstrap": boot.round(3).to_dict(orient="records"),
        "bootstrap_params": {"n_boot": N_BOOT, "noise_fraction": NOISE, "top_k": TOPK, "pool": "approved drugs"},
        "map_method": method,
        "map_signatures": [{"id": s, "x": float(p[0]), "y": float(p[1])} for s, p in zip(sig_ids, sxy)],
    }
    (C.RESULTS / "robustness.json").write_text(json.dumps(out, indent=1))
    print(f"done in {(time.time() - t0) / 60:.1f} min (map: {method})")


if __name__ == "__main__":
    main()
