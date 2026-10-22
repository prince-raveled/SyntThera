"""Step 5: validation of the screen.

* Positive-control recovery: where known exercise mimetics land among ~29,700 compounds,
  with an empirical p-value (random drug sets of the same size).
* Negative controls: mTOR inhibitors should oppose the resistance signature.
* Permutation significance: how many compounds pass FDR < 0.05 against random queries.
* Agreement: between the four scoring methods, and between tissues / timepoints.
Writes results/validation.json, CSV tables and PNG figures for slides.
"""
import json
import sys
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from tejas import config as C  # noqa: E402

FIG = C.RESULTS / "figures"
FIG.mkdir(exist_ok=True)
METHODS = ["wtcs", "spearman", "xsum", "gsfm"]


def control_table(res, n, names, modality, direction):
    rank_col = "rank" if direction == "mimic" else "inv_rank"
    sub = res[(res.modality == modality) & res.drug_id.isin(names)]
    out = sub[["signature", "tissue", "timepoint", "drug_id", rank_col, "perm_fdr", "benefit"]].copy()
    out["percentile"] = out[rank_col] / n * 100
    out["direction"] = direction
    return out.rename(columns={rank_col: "rank"})


def empirical_p(res_sig: pd.DataFrame, k: int, observed_mean_pct: float, n: int, rng, draws=20000):
    pct = res_sig["rank"].to_numpy() / n * 100
    # k random compounds per draw (duplicates are negligible for k << n)
    sims = pct[rng.integers(0, len(pct), size=(draws, k))].mean(axis=1)
    return float((np.sum(sims <= observed_mean_pct) + 1) / (draws + 1))


def main():
    res = pd.read_parquet(C.RESULTS / "scores_long.parquet")
    cpd = pd.read_parquet(C.PROC / "compounds.parquet")
    meta = json.loads((C.RESULTS / "scoring_meta.json").read_text())
    informative = meta["informative_signatures"]
    n = len(cpd)
    rng = np.random.default_rng(C.SEED)
    launched = set(cpd.loc[cpd.clinical_phase == "Launched", "drug_id"])

    # ---- controls ----
    tabs = []
    for mod, names in C.POSITIVE_CONTROLS.items():
        tabs.append(control_table(res, n, names, mod, "mimic"))
    for mod, names in C.NEGATIVE_CONTROLS.items():
        tabs.append(control_table(res, n, names, mod, "oppose"))
    ctrl = pd.concat(tabs, ignore_index=True)
    ctrl["informative"] = ctrl.signature.isin(informative)
    ctrl.to_csv(C.RESULTS / "validation_controls.csv", index=False)

    summary = []
    for sig_id in informative:
        tissue, mod, tp = sig_id.split("|")
        rs = res[res.signature == sig_id]
        for direction, names in (("mimic", C.POSITIVE_CONTROLS.get(mod, [])),
                                 ("oppose", C.NEGATIVE_CONTROLS.get(mod, []))):
            if not names:
                continue
            g = ctrl[(ctrl.signature == sig_id) & (ctrl.direction == direction)]
            if g.empty:
                continue
            col = "rank" if direction == "mimic" else "inv_rank"
            ref = rs[[col]].rename(columns={col: "rank"})
            m = g.percentile.mean()
            summary.append({"signature": sig_id, "tissue": tissue, "modality": mod, "timepoint": tp,
                            "direction": direction, "n_controls": len(g), "mean_percentile": round(m, 1),
                            "best_control": g.sort_values("rank").drug_id.iloc[0],
                            "best_percentile": round(g.percentile.min(), 1),
                            "p_vs_random": empirical_p(ref, len(g), m, n, rng)})
    ctrl_summary = pd.DataFrame(summary)
    print("control recovery (informative signatures):")
    print(ctrl_summary.to_string(index=False))

    # ---- permutation significance ----
    perm = (res.assign(sig=res.perm_fdr < 0.05, sig_launched=(res.perm_fdr < 0.05) & res.drug_id.isin(launched))
              .groupby("signature")[["sig", "sig_launched"]].sum().astype(int).reset_index()
              .rename(columns={"sig": "n_fdr05", "sig_launched": "n_fdr05_launched"}))
    perm["informative"] = perm.signature.isin(informative)
    print("\ncompounds passing permutation FDR < 0.05:")
    print(perm.to_string(index=False))

    # ---- method agreement (informative signatures) ----
    mats = []
    for s in informative:
        rs = res[res.signature == s]
        mats.append(spearmanr(rs[METHODS].to_numpy()).statistic)
    method_corr = np.mean(mats, axis=0)

    # ---- signature agreement (combined score, all compounds) ----
    res["score"] = -np.log10(res.rho.clip(lower=1e-300))
    wide = res.pivot(index="drug_id", columns="signature", values="score")
    order = sorted(wide.columns, key=lambda s: (s.split("|")[1], s.split("|")[0], ["early", "3.5-4h", "24h"].index(s.split("|")[2])))
    wide = wide[order]
    sig_corr = spearmanr(wide.to_numpy()).statistic

    # ---- figures ----
    fig, ax = plt.subplots(figsize=(8, 3.6))
    pos = ctrl[ctrl.informative & (ctrl.direction == "mimic")]
    for i, (name, g) in enumerate(pos.groupby("drug_id", sort=False)):
        ax.scatter(g.percentile, [i] * len(g), s=36, alpha=0.8)
    ax.set_yticks(range(pos.drug_id.nunique()))
    ax.set_yticklabels(list(dict.fromkeys(pos.drug_id)))
    ax.axvline(50, color="grey", lw=0.8, ls="--")
    ax.set_xlim(0, 100)
    ax.set_xlabel("percentile among all compounds (lower = stronger exercise match)")
    ax.set_title("Known exercise mimetics across informative human signatures")
    fig.tight_layout()
    fig.savefig(FIG / "positive_controls.png", dpi=160)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.5, 6.5))
    im = ax.imshow(sig_corr, cmap="RdBu_r", vmin=-1, vmax=1)
    labels = [s.replace("|", " · ") for s in order]
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(labels, rotation=90, fontsize=7)
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(labels, fontsize=7)
    fig.colorbar(im, ax=ax, shrink=0.7, label="Spearman r of drug scores")
    ax.set_title("Agreement between exercise signatures")
    fig.tight_layout()
    fig.savefig(FIG / "signature_agreement.png", dpi=160)
    plt.close(fig)

    out = {
        "n_compounds": n,
        "informative_signatures": informative,
        "controls_summary": ctrl_summary.to_dict(orient="records"),
        "controls": ctrl.round(4).to_dict(orient="records"),
        "permutation": perm.to_dict(orient="records"),
        "method_agreement": {"methods": METHODS, "spearman": np.round(method_corr, 3).tolist()},
        "signature_agreement": {"order": order, "spearman": np.round(sig_corr, 3).tolist()},
    }
    (C.RESULTS / "validation.json").write_text(json.dumps(out, indent=1))
    print("\nmean method agreement (Spearman):")
    print(pd.DataFrame(method_corr, index=METHODS, columns=METHODS).round(2).to_string())


if __name__ == "__main__":
    main()
