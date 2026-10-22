"""Step 7: rodent comparison arm (the conventional pipeline the proposal argues against).

Rat gastrocnemius after 8 weeks of endurance training (MoTrPAC PASS1B, both sexes,
Stouffer-combined z), converted to human genes keeping 1:1 orthologs only, then scored
with exactly the same four methods + RRA as the human signatures.
Reports how far the drug rankings diverge and how many genes the conversion loses.
Caveat: rat = 8-week training, human = single bout, so species and acute-vs-chronic are confounded.
"""
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyreadr
from scipy.stats import norm, spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from tejas import config as C  # noqa: E402
from tejas import scoring as S  # noqa: E402

spec = importlib.util.spec_from_file_location("score04", Path(__file__).with_name("04_score.py"))
score04 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(score04)

RAT = C.RAW / "motrpac_rat"
HUMAN_REF = ["muscle|endurance|3.5-4h", "muscle|endurance|24h"]


def bh(p):
    return score04.bh(np.asarray(p, dtype=float))


def rat_signature() -> pd.DataFrame:
    da = next(iter(pyreadr.read_r(RAT / "TRNSCRPT_SKMGN_DA.rda").values()))
    da = da[da.comparison_group == "8w"]
    wide = da.pivot_table(index="feature_ID", columns="sex", values="zscore")
    wide = wide.dropna()
    z = (wide["male"] + wide["female"]) / np.sqrt(2)
    p = 2 * norm.sf(np.abs(z))
    return pd.DataFrame({"rat_ensembl": z.index, "z": z.to_numpy(), "fdr": bh(p)})


def orthologs() -> pd.DataFrame:
    o = next(iter(pyreadr.read_r(RAT / "RAT_TO_HUMAN_GENE.rda").values()))
    o = o[["RAT_ENSEMBL_ID", "HUMAN_ORTHOLOG_NCBI_GENE_ID"]].dropna()
    o = o.assign(rat_ensembl=o.RAT_ENSEMBL_ID.astype(str).str.split(";")).explode("rat_ensembl")
    o["entrez"] = o.HUMAN_ORTHOLOG_NCBI_GENE_ID.astype(str).str.replace(r"\.0$", "", regex=True)
    o = o[["rat_ensembl", "entrez"]].drop_duplicates()
    one_h = o.groupby("rat_ensembl").entrez.transform("nunique") == 1
    one_r = o.groupby("entrez").rat_ensembl.transform("nunique") == 1
    return o[one_h & one_r]


def main():
    sc = score04.Scorer()
    rat = rat_signature()
    orth = orthologs()
    n_rat = len(rat)
    n_rat_sig = int((rat.fdr < C.FDR).sum())
    m = rat.merge(orth, on="rat_ensembl", how="inner")
    m = m[m.entrez.isin(sc.idx_of)]
    lost_sig = 1 - (m.fdr < C.FDR).sum() / max(n_rat_sig, 1)
    print(f"rat genes {n_rat:,} (FDR<0.05: {n_rat_sig:,}); after 1:1 orthologs in L1000 genes: {len(m):,} "
          f"({(m.fdr < C.FDR).sum():,} significant) -> lost {lost_sig:.1%} of significant genes")

    m = m.assign(gidx=m.entrez.map(sc.idx_of).astype(int), t=m.z, gene_symbol=sc.symbols[m.entrez.map(sc.idx_of).astype(int)])
    up, dn, fallback = score04.query_sets(m)
    t = np.full(sc.G, np.nan)
    t[m.gidx.to_numpy()] = m.z.to_numpy()
    up_idx, dn_idx = up.gidx.to_numpy(), dn.gidx.to_numpy()
    # identical procedure to the human signatures, including per-compound null standardisation
    r_ = score04.score_signature(sc, "rat_muscle", t, up_idx, dn_idx, np.random.default_rng(C.SEED))
    rho_rat = r_["rho"]

    res = pd.read_parquet(C.RESULTS / "scores_long.parquet")
    cpd = sc.cpd.assign(rat_rho=rho_rat, rat_rank=S.final_rank(rho_rat, r_["ranks"].mean(axis=1)))
    cpd.loc[:, ["drug_id", "rat_rho", "rat_rank"]].to_parquet(C.RESULTS / "rat_scores.parquet", index=False)
    approved = cpd.clinical_phase == "Launched"

    comps = {}
    for ref in HUMAN_REF:
        h = res[res.signature == ref].set_index("drug_id").loc[cpd.drug_id]
        rs = spearmanr(-np.log10(h.rho.clip(lower=1e-300)), -np.log10(cpd.rat_rho.clip(lower=1e-300))).statistic
        top_h, top_r = set(h.nsmallest(100, "rank").index), set(cpd.nsmallest(100, "rat_rank").drug_id)
        ah = h[approved.to_numpy()].nsmallest(100, "rank")
        ar = cpd[approved].nsmallest(100, "rat_rank")
        comps[ref] = {"spearman": float(rs), "top100_overlap": len(top_h & top_r),
                      "top100_overlap_approved": len(set(ah.index) & set(ar.drug_id)),
                      "top_human": list(cpd.set_index("drug_id").loc[ah.index[:10], "name"]),
                      "top_rat": list(ar.name[:10])}
        print(f"vs {ref}: Spearman {rs:.3f}, top-100 overlap {comps[ref]['top100_overlap']}, "
              f"approved top-100 overlap {comps[ref]['top100_overlap_approved']}")

    ctrl = {}
    for name in C.POSITIVE_CONTROLS["endurance"]:
        if name in set(cpd.drug_id):
            hr = int(res[(res.signature == HUMAN_REF[0]) & (res.drug_id == name)]["rank"].iloc[0])
            rr = int(cpd.loc[cpd.drug_id == name, "rat_rank"].iloc[0])
            ctrl[name] = {"human_pct": round(hr / len(cpd) * 100, 1), "rat_pct": round(rr / len(cpd) * 100, 1)}
    print("endurance controls (percentile, lower = better):", ctrl)

    main_ref = comps[HUMAN_REF[0]]
    out = {
        "title": "Rat gastrocnemius, 8 weeks endurance training vs human muscle, 3.5–4 h after endurance exercise",
        "spearman": main_ref["spearman"], "top100_overlap": main_ref["top100_overlap"],
        "top100_overlap_approved": main_ref["top100_overlap_approved"],
        "genes_lost_frac": float(lost_sig), "n_rat_sig": n_rat_sig, "query_fallback": bool(fallback),
        "vs_24h": comps[HUMAN_REF[1]], "controls": ctrl,
        "top_human": main_ref["top_human"], "top_rat": main_ref["top_rat"],
        "text": (f"Converting the rat signature to human genes kept only clean 1:1 orthologs present in L1000, "
                 f"losing {lost_sig:.0%} of the {n_rat_sig:,} significant rat genes. Across all {len(cpd):,} compounds the "
                 f"two rankings correlate at r = {main_ref['spearman']:.2f}, and {main_ref['top100_overlap']} of the top 100 "
                 f"compounds are shared (r = {comps[HUMAN_REF[1]]['spearman']:.2f} against the human 24 h signature)."),
        "caveat": ("The rat data are 8 weeks of training and the human data a single bout, so this difference mixes "
                   "species with acute-versus-chronic exercise. It measures how much the drug list depends on the "
                   "starting signature, which is the assumption earlier rodent-based screens left untested."),
    }
    (C.RESULTS / "rat_comparison.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
