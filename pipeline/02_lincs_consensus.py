"""Step 2: read the 33 GB LINCS matrix once and build one consensus signature per compound.

* Integrity check of the GCTX file (size, HDF5 layout, first/last/random signatures readable).
* Keep exemplar compound signatures (one per compound x cell line; 136,460 of 720,216).
* Keep the L1000 "BING" gene space (978 landmark + 9,196 best-inferred genes).
* Collapse each compound's cell-line signatures with MODZ (CMap's weighted average:
  weights = summed Spearman correlation to the compound's other signatures, landmark genes).
* Annotate compounds with LINCS compoundinfo and Drug Repurposing Hub clinical phase.
"""
import sys
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from tejas import config as C  # noqa: E402

BLOCK = 4096


def decode(arr):
    return np.array([x.decode().strip() if isinstance(x, bytes) else str(x).strip() for x in arr])


def integrity_check(f, ds):
    size = C.GCTX.stat().st_size
    print(f"file size {size:,} bytes (expected {C.GCTX_BYTES:,}) -> {'OK' if size == C.GCTX_BYTES else 'MISMATCH'}")
    print(f"matrix shape {ds.shape}, dtype {ds.dtype}, chunks {ds.chunks}, compression {ds.compression}")
    rng = np.random.default_rng(C.SEED)
    idx = sorted({0, ds.shape[0] - 1, *rng.integers(0, ds.shape[0], 20).tolist()})
    bad = 0
    for i in idx:
        v = ds[i, :]
        if not np.all(np.isfinite(v)) or float(np.std(v)) == 0.0:
            bad += 1
            print(f"  row {i}: suspicious (finite={np.all(np.isfinite(v))}, sd={np.std(v):.3f})")
    print(f"read test: {len(idx)} signatures incl. first and last -> {len(idx) - bad} OK")
    if bad:
        raise SystemExit("GCTX integrity check failed")


def modz(block: np.ndarray, lm_cols: np.ndarray) -> np.ndarray:
    """MODZ consensus of k signatures (k x G)."""
    k = block.shape[0]
    if k == 1:
        return block[0]
    if k == 2:
        return block.mean(axis=0)
    r = np.argsort(np.argsort(block[:, lm_cols], axis=1), axis=1).astype(np.float32)
    rho = np.corrcoef(r)
    np.fill_diagonal(rho, 0.0)
    w = np.clip(rho.sum(axis=1), 0.01, None)
    w /= w.sum()
    return w @ block


def join_unique(s) -> str:
    return "; ".join(sorted({str(x) for x in s.dropna() if str(x) not in ("", "nan")}))


def drug_names() -> pd.Series:
    """pert_id -> drug name. One drug can have several LINCS pert_ids (vendors, batches, salts)."""
    ci = pd.read_csv(C.LINCS / "compoundinfo_beta.txt", sep="\t")
    return ci.groupby("pert_id")["cmap_name"].first()


def annotate(cpd: pd.DataFrame, si: pd.DataFrame) -> pd.DataFrame:
    ci = pd.read_csv(C.LINCS / "compoundinfo_beta.txt", sep="\t")
    ci["drug_id"] = ci["pert_id"].map(drug_names())
    lincs = ci.groupby("drug_id").agg(lincs_moa=("moa", join_unique), lincs_target=("target", join_unique))
    cpd = cpd.merge(lincs, left_on="drug_id", right_index=True, how="left")

    def read_hub(name):
        path = C.HUB / name
        with open(path, encoding="utf-8", errors="replace") as fh:
            skip = 0
            for line in fh:
                if not line.startswith("!"):
                    break
                skip += 1
        return pd.read_csv(path, sep="\t", skiprows=skip, encoding="utf-8", encoding_errors="replace")

    drugs = read_hub("repurposing_drugs_20200324.txt")
    drugs["key"] = drugs["pert_iname"].str.lower()
    samples = read_hub("repurposing_samples_20200324.txt")
    samples["pert_id"] = samples["broad_id"].str[:13]
    by_brd = samples.drop_duplicates("pert_id").set_index("pert_id")["pert_iname"].str.lower()

    # match each drug to the Hub by name first, then through any of its Broad IDs
    hub_keys = set(drugs["key"])
    pert_lists = si.groupby("drug_id")["pert_id"].agg(lambda s: sorted(set(s)))

    def hub_key(row):
        k = str(row.drug_id).lower()
        if k in hub_keys:
            return k
        for p in pert_lists.get(row.drug_id, []):
            if p in by_brd.index and by_brd[p] in hub_keys:
                return by_brd[p]
        return None

    cpd["key"] = cpd.apply(hub_key, axis=1)
    hub = drugs.drop_duplicates("key").rename(columns={
        "pert_iname": "hub_name", "moa": "hub_moa", "target": "hub_target"})
    cpd = cpd.merge(hub, on="key", how="left").drop(columns=["key"])
    cpd["clinical_phase"] = cpd["clinical_phase"].fillna("Not in Hub")
    cpd["name"] = cpd["hub_name"].fillna(cpd["drug_id"])
    cpd["moa"] = cpd["hub_moa"].fillna(cpd["lincs_moa"]).replace("", np.nan)
    cpd["target"] = cpd["hub_target"].fillna(cpd["lincs_target"]).replace("", np.nan)
    return cpd


def main():
    t0 = time.time()
    gi = pd.read_csv(C.LINCS / "geneinfo_beta.txt", sep="\t", dtype={"gene_id": "string"})
    si = pd.read_csv(C.LINCS / "siginfo_beta.txt", sep="\t", low_memory=False,
                     usecols=["sig_id", "pert_id", "pert_type", "cell_iname", "is_exemplar_sig",
                              "is_hiq", "tas", "pert_idose", "pert_itime"])
    si = si[(si.pert_type == "trt_cp") & (si.is_exemplar_sig == 1)].copy()
    si["drug_id"] = si["pert_id"].map(drug_names()).fillna(si["pert_id"])

    with h5py.File(C.GCTX, "r") as f:
        ds = f["0/DATA/0/matrix"]
        col_ids = decode(f["0/META/COL/id"][:])
        row_ids = decode(f["0/META/ROW/id"][:])
        assert ds.shape == (len(col_ids), len(row_ids)), "unexpected GCTX orientation"
        integrity_check(f, ds)

        genes = pd.DataFrame({"gene_id": row_ids, "col": np.arange(len(row_ids))})
        genes = genes.merge(gi[["gene_id", "gene_symbol", "feature_space"]], on="gene_id")
        genes = genes[genes.feature_space.isin(["landmark", "best inferred"])].sort_values("col")
        gcols = genes["col"].to_numpy()
        genes = genes.reset_index(drop=True)
        genes["idx"] = np.arange(len(genes))
        print(f"gene universe (landmark + best inferred): {len(genes):,}")

        pos = pd.Series(np.arange(len(col_ids)), index=col_ids)
        si = si[si.sig_id.isin(pos.index)].copy()
        si["row"] = pos.loc[si.sig_id].to_numpy()
        si = si.sort_values("row").reset_index(drop=True)
        rows = si["row"].to_numpy()
        print(f"exemplar compound signatures: {len(si):,} from {si.pert_id.nunique():,} LINCS IDs "
              f"= {si.drug_id.nunique():,} distinct drugs")

        tmp = C.PROC / "_exemplar_z_f16.npy"
        ex = np.lib.format.open_memmap(tmp, mode="w+", dtype=np.float16, shape=(len(si), len(genes)))
        n_blocks = int(np.ceil(ds.shape[0] / BLOCK))
        done = 0
        for b in range(n_blocks):
            lo, hi = b * BLOCK, min((b + 1) * BLOCK, ds.shape[0])
            a, z = np.searchsorted(rows, lo), np.searchsorted(rows, hi)
            if a == z:
                continue
            blk = ds[lo:hi, :]
            ex[a:z] = blk[rows[a:z] - lo][:, gcols].astype(np.float16)
            done = z
            if b % 20 == 0 or b == n_blocks - 1:
                el = time.time() - t0
                print(f"  block {b + 1}/{n_blocks}  signatures kept {done:,}/{len(si):,}  {el / 60:.1f} min", flush=True)
        ex.flush()

    # MODZ consensus per drug (pooling all of its LINCS IDs and cell lines)
    lm_cols = genes.index[genes.feature_space == "landmark"].to_numpy()
    si["grp"] = si.groupby("drug_id", sort=True).ngroup()
    order = np.argsort(si["grp"].to_numpy(), kind="stable")
    grp_sorted = si["grp"].to_numpy()[order]
    starts = np.r_[0, np.flatnonzero(np.diff(grp_sorted)) + 1]
    ends = np.r_[starts[1:], len(order)]
    n_cpd = len(starts)
    cons = np.zeros((n_cpd, len(genes)), dtype=np.float32)
    for i, (s, e) in enumerate(zip(starts, ends)):
        cons[i] = modz(ex[np.sort(order[s:e])].astype(np.float32), lm_cols)
        if i % 5000 == 0:
            print(f"  MODZ {i:,}/{n_cpd:,}", flush=True)
    np.save(C.PROC / "lincs_consensus_z.npy", cons)

    g = si.groupby("grp")
    cpd = pd.DataFrame({
        "drug_id": g["drug_id"].first(),
        "pert_ids": g["pert_id"].agg(join_unique),
        "n_sigs": g.size(),
        "n_hiq": g["is_hiq"].sum(),
        "n_cells": g["cell_iname"].nunique(),
        "mean_tas": g["tas"].mean().round(3),
        "cells": g["cell_iname"].agg(lambda s: ", ".join(sorted(s.unique()))),
    }).sort_index().reset_index(drop=True)
    cpd = annotate(cpd, si)
    cpd.to_parquet(C.PROC / "compounds.parquet", index=False)
    genes.drop(columns=["col"]).to_parquet(C.PROC / "lincs_genes.parquet", index=False)
    del ex  # release the memory map so Windows lets us delete the temp file
    tmp.unlink(missing_ok=True)

    print(f"consensus matrix {cons.shape} saved; compounds by phase:")
    print(cpd["clinical_phase"].value_counts().to_string())
    print(f"done in {(time.time() - t0) / 60:.1f} min")


if __name__ == "__main__":
    main()
