"""Step 1: build the 18 human exercise signatures from MoTrPAC differential analysis.

Each signature is a "difference-in-changes" contrast from the MoTrPAC dream model:
(exercise group post - pre) - (non-exercise control post - pre), for
3 tissues x 2 modalities (endurance, resistance) x 3 post-exercise timepoints.
Genes are mapped Ensembl -> Entrez so they can be joined to LINCS L1000.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from tejas import config as C  # noqa: E402


def contrast_name(group: str, tp: str) -> str:
    return (f"group_timepoint{group}.{tp} - group_timepoint{group}.pre_exercise - "
            f"group_timepointADUControl.{tp} + group_timepointADUControl.pre_exercise")


def load_tissue(tissue: str, code: str) -> pd.DataFrame:
    da_path = C.MOTRPAC / f"human-precovid-sed-adu_{code}_transcript-rna-seq_da_dream-acute_v2.1.txt"
    fx_path = C.MOTRPAC / f"human-precovid-sed-adu_{code}_transcript-rna-seq_metadata_features_v2.0.txt"
    wanted = {contrast_name(g, tp): (mod, tpk)
              for mod, g in C.MODALITIES.items() for tpk, tp in C.TIMEPOINTS.items()}
    da = pd.read_csv(da_path, sep="\t", usecols=["feature_id", "logFC", "t", "adj_p_value", "z.std", "contrast"])
    da = da[da["contrast"].isin(wanted)].copy()
    missing = set(wanted) - set(da["contrast"].unique())
    if missing:
        raise SystemExit(f"{tissue}: missing contrasts {missing}")
    da["modality"] = da["contrast"].map(lambda c: wanted[c][0])
    da["timepoint"] = da["contrast"].map(lambda c: wanted[c][1])

    fx = pd.read_csv(fx_path, sep="\t", dtype={"entrez_gene": "string"})
    fx = fx.dropna(subset=["entrez_gene"])
    fx["entrez"] = fx["entrez_gene"].str.split(r"[;,|]").str[0].str.strip()
    fx = fx[fx["entrez"].str.fullmatch(r"\d+")]
    da = da.merge(fx[["feature_id", "entrez", "gene_symbol"]], on="feature_id", how="inner")
    da["gene_symbol"] = da["gene_symbol"].fillna(da["entrez"])
    # several Ensembl IDs can map to one Entrez gene: keep the strongest
    da["abs_t"] = da["t"].abs()
    da = (da.sort_values("abs_t", ascending=False)
            .drop_duplicates(["modality", "timepoint", "entrez"])
            .drop(columns=["abs_t", "contrast"]))
    da.insert(0, "tissue", tissue)
    return da


def main():
    gi = pd.read_csv(C.LINCS / "geneinfo_beta.txt", sep="\t", dtype={"gene_id": "string"})
    space = dict(zip(gi["gene_id"], gi["feature_space"]))

    parts = [load_tissue(t, code) for t, code in C.TISSUES.items()]
    sig = pd.concat(parts, ignore_index=True)
    sig["entrez"] = sig["entrez"].astype("string")
    sig["l1000_space"] = sig["entrez"].map(space).fillna("not in L1000")
    sig["signature"] = [C.signature_id(t, m, p) for t, m, p in zip(sig.tissue, sig.modality, sig.timepoint)]
    sig = sig.rename(columns={"z.std": "z_std", "adj_p_value": "fdr"})
    sig.to_parquet(C.PROC / "signatures_long.parquet", index=False)

    rows = []
    for (t, m, p), g in sig.groupby(["tissue", "modality", "timepoint"], sort=False):
        sigf = g[g.fdr < C.FDR]
        rows.append({
            "signature": C.signature_id(t, m, p), "tissue": t, "modality": m, "timepoint": p,
            "timepoint_label": C.timepoint_label(t, p),
            "n_genes": len(g), "n_in_l1000_bing": int(g.l1000_space.isin(["landmark", "best inferred"]).sum()),
            "n_fdr_up": int((sigf.t > 0).sum()), "n_fdr_down": int((sigf.t < 0).sum()),
            "top_up": ", ".join(g.nlargest(5, "t").gene_symbol.astype(str)),
            "top_down": ", ".join(g.nsmallest(5, "t").gene_symbol.astype(str)),
        })
    summ = pd.DataFrame(rows)
    summ.to_csv(C.PROC / "signature_summary.csv", index=False)
    with pd.option_context("display.width", 200, "display.max_columns", 20):
        print(summ.drop(columns=["signature"]).to_string(index=False))


if __name__ == "__main__":
    main()
