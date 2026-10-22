"""Step 6: export results into the static website (website/index.html, data embedded).

The page is a single self-contained HTML file: it opens from disk, GitHub Pages, or any host.
"""
import json
import sys
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from tejas import config as C  # noqa: E402
from tejas.candidates import MATCH_TOP, MIN_CELLS, candidate_table, flag_cytotoxic  # noqa: E402

TOP_PER_SIG = 100
TP_ORDER = ["early", "3.5-4h", "24h"]


def r(x, nd=2):
    return None if x is None or (isinstance(x, float) and not np.isfinite(x)) else round(float(x), nd)


def sig_fdr(x):
    return float(f"{x:.2g}")


def main():
    res = pd.read_parquet(C.RESULTS / "scores_long.parquet")
    cpd = pd.read_parquet(C.PROC / "compounds.parquet")
    summ = pd.read_csv(C.PROC / "signature_summary.csv")
    sigl = pd.read_parquet(C.PROC / "signatures_long.parquet")
    queries = json.loads((C.RESULTS / "query_sets.json").read_text())
    smeta = json.loads((C.RESULTS / "scoring_meta.json").read_text())
    val = json.loads((C.RESULTS / "validation.json").read_text())
    overall = pd.read_parquet(C.RESULTS / "overall.parquet")
    mods = pd.read_parquet(C.PROC / "gene_class.parquet")
    gclass = dict(zip(mods.gene_symbol, mods.gene_class))
    informative = set(smeta["informative_signatures"])
    n = len(cpd)

    # ---- which drugs go into the page ----
    keep = set(cpd.loc[cpd.clinical_phase == "Launched", "drug_id"])
    keep |= set(res[res["rank"] <= TOP_PER_SIG].drug_id)
    keep |= set(overall[overall["rank"] <= TOP_PER_SIG].drug_id)
    for names in list(C.POSITIVE_CONTROLS.values()) + list(C.NEGATIVE_CONTROLS.values()):
        keep |= set(names) & set(cpd.drug_id)
    d = cpd[cpd.drug_id.isin(keep)].reset_index(drop=True)
    d["cyt"] = flag_cytotoxic(d.moa, d.target)
    pos = {k: i for i, k in enumerate(d.drug_id)}
    controls = {name: kind for kind, groups in (("positive", C.POSITIVE_CONTROLS), ("negative", C.NEGATIVE_CONTROLS))
                for names in groups.values() for name in names}

    def txt(v, limit=160):
        if v is None or (isinstance(v, float) and np.isnan(v)) or str(v) in ("", "nan", "None"):
            return ""
        s = str(v).replace("|", ", ")
        return s if len(s) <= limit else s[:limit - 1] + "…"

    drugs = [{"id": row.drug_id, "name": row["name"], "phase": row.clinical_phase, "moa": txt(row.moa),
              "target": txt(row.target, 120), "disease": txt(row.get("disease_area"), 80),
              "cells": int(row.n_cells), "sigs": int(row.n_sigs), "tas": r(row.mean_tas),
              "control": controls.get(row.drug_id, ""), "cyt": bool(row.cyt),
              "use": txt(row.get("indication"), 90)}
             for _, row in d.iterrows()]

    # ---- signatures ----
    sig_order = [C.signature_id(t, m, p) for m in C.MODALITIES for t in C.TISSUES for p in TP_ORDER]
    sigs, scores = [], {}
    perm = {p["signature"]: p for p in val["permutation"]}
    for sid in sig_order:
        t, m, p = sid.split("|")
        s = summ[summ.signature == sid].iloc[0]
        q = queries[sid]
        g = sigl[sigl.signature == sid]
        tv = dict(zip(g.gene_symbol.astype(str), g.t))
        g = g.assign(cls=g.gene_symbol.map(gclass).fillna(""))
        module_shift = {c: r(g.loc[g.cls == c, "t"].mean()) for c in ("adaptive", "stress")}
        sigs.append({
            "id": sid, "tissue": t, "modality": m, "timepoint": p, "label": C.timepoint_label(t, p),
            "informative": sid in informative, "fallback": q["fallback"],
            "n_genes": int(s.n_genes), "n_up": int(s.n_fdr_up), "n_down": int(s.n_fdr_down),
            "n_drugs_fdr": int(perm[sid]["n_fdr05"]), "n_launched_fdr": int(perm[sid]["n_fdr05_launched"]),
            "up": [[gname, r(tv.get(gname), 1)] for gname in q["up"][:30]],
            "down": [[gname, r(tv.get(gname), 1)] for gname in q["down"][:30]],
            "query_sizes": [len(q["up"]), len(q["down"])],
            "module_t": module_shift,
        })
        rs = res[(res.signature == sid) & res.drug_id.isin(keep)].copy()
        rs["i"] = rs.drug_id.map(pos)
        rs = rs.sort_values("i")
        scores[sid] = {
            "rank": rs["rank"].astype(int).tolist(),
            "inv": rs["inv_rank"].astype(int).tolist(),
            "fdr": [sig_fdr(x) for x in rs.perm_fdr],
            "ben": [r(x) for x in rs.benefit],
            "ad": [r(x) for x in rs.adaptive_z],
            "st": [r(x) for x in rs.stress_z],
            "m": np.round(rs[["r_wtcs", "r_spearman", "r_xsum", "r_gsfm"]].to_numpy() * 100, 1).tolist(),
        }

    ov = {}
    for m, g in overall[overall.drug_id.isin(keep)].groupby("modality"):
        g = g.assign(i=g.drug_id.map(pos)).sort_values("i")
        ov[m] = {"rank": g["rank"].astype(int).tolist(), "median_rank": g.median_rank.round(0).astype(int).tolist()}

    # ---- human response context: canonical gene trajectories and volcano plots (muscle) ----
    traj_genes = ["PPARGC1A", "NR4A3", "VEGFA", "PDK4", "HK2", "ESRRA", "TFAM", "ATF3", "MYC", "IL6"]
    traj = {"tissue": "muscle", "times": ["pre", "early", "3.5-4h", "24h"], "genes": [], "lfc": {}, "fdr": {}}
    musc = sigl[sigl.tissue == "muscle"]
    for gname in traj_genes:
        gg = musc[musc.gene_symbol.astype(str) == gname]
        if gg.empty:
            continue
        traj["genes"].append(gname)
        for mod in C.MODALITIES:
            row = gg[gg.modality == mod].set_index("timepoint")
            traj["lfc"].setdefault(mod, {})[gname] = [0.0] + [r(row.logFC.get(p)) for p in TP_ORDER]
            traj["fdr"].setdefault(mod, {})[gname] = [1.0] + [sig_fdr(row.fdr.get(p, 1.0)) for p in TP_ORDER]
    volcano = {}
    for mod in C.MODALITIES:
        sid_v = C.signature_id("muscle", mod, "3.5-4h")
        gv = sigl[sigl.signature == sid_v]
        q = -np.log10(gv.fdr.clip(lower=1e-30))
        lab = gv.assign(q=q.to_numpy()).nlargest(14, "q")
        volcano[mod] = {"lfc": np.round(gv.logFC.to_numpy(), 2).tolist(), "q": np.round(q.to_numpy(), 2).tolist(),
                        "labels": [[r(a), r(b), str(c)] for a, b, c in zip(lab.logFC, lab.q, lab.gene_symbol)]}

    # ---- hero animation: strongest muscle genes after endurance exercise, across the time course ----
    hs = sigl[(sigl.signature == C.signature_id("muscle", "endurance", "3.5-4h")) & (sigl.fdr < C.FDR)]
    hero_genes = list(hs.nlargest(9, "t").gene_symbol.astype(str)) + list(hs.nsmallest(4, "t").gene_symbol.astype(str))
    me = musc[musc.modality == "endurance"]
    hero = {"genes": hero_genes, "lfc": {}}
    for p in TP_ORDER:
        sub = me[me.timepoint == p].drop_duplicates("gene_symbol")
        lf = dict(zip(sub.gene_symbol.astype(str), sub.logFC))
        hero["lfc"][p] = [r(lf.get(gname, 0.0)) for gname in hero_genes]

    # ---- drug map and robustness ----
    rob_path, map_path = C.RESULTS / "robustness.json", C.RESULTS / "drug_map.parquet"
    rob = json.loads(rob_path.read_text()) if rob_path.exists() else None
    dmap = None
    if map_path.exists():
        mp = pd.read_parquet(map_path).set_index("drug_id")
        page = mp.loc[d.drug_id]
        rest = mp.drop(index=d.drug_id).sample(n=min(7000, len(mp) - len(d)), random_state=C.SEED)
        dmap = {"method": rob["map_method"] if rob else "", "x": np.round(page.x.to_numpy(), 3).tolist(),
                "y": np.round(page.y.to_numpy(), 3).tolist(),
                "bg": np.round(rest[["x", "y"]].to_numpy(), 3).ravel().tolist(),
                "sigs": rob["map_signatures"] if rob else []}

    outputs = []
    for f in sorted(list(C.RESULTS.glob("*.tsv")) + list(C.RESULTS.glob("*.jsonld")) + list(C.RESULTS.glob("*.parquet"))):
        outputs.append({"file": f"results/{f.name}", "mb": round(f.stat().st_size / 1e6, 2)})

    cands = {"rule": {"match_top": MATCH_TOP, "min_cells": MIN_CELLS}}
    for mod in C.MODALITIES:
        ct = candidate_table(res, cpd, overall, sorted(informative), mod)
        ct = ct[ct.approved & ct.strong_match]
        cands[mod] = [{"i": pos[row.drug_id], "c": None if np.isnan(row.candidate_rank) else int(row.candidate_rank),
                       "m": int(row.match_rank), "b": r(row.mean_benefit), "a": r(row.mean_adaptive_z),
                       "s": r(row.mean_stress_z), "f": int(row.n_signatures_fdr05), "cyt": bool(row.flag_cytotoxic),
                       "low": bool(row.flag_low_evidence), "best": row.best_signature}
                      for _, row in ct.iterrows() if row.drug_id in pos]

    rat_path = C.RESULTS / "rat_comparison.json"
    data = {
        "meta": {
            "generated": date.today().isoformat(), "n_compounds": n,
            "n_launched": int((cpd.clinical_phase == "Launched").sum()),
            "n_drugs_in_page": len(drugs), "n_perm": smeta["n_perm"],
            "gsfm_coverage": r(smeta["gsfm_vocab_coverage_of_l1000_genes"], 3),
            "releases": {"motrpac": C.MOTRPAC_RELEASE, "lincs": C.LINCS_RELEASE, "hub": C.HUB_RELEASE,
                         "gsfm": C.GSFM_MODEL},
            "modules": json.loads((C.PROC / "gene_modules_meta.json").read_text()),
            "n_query": C.N_QUERY,
        },
        "sigs": sigs, "drugs": drugs, "scores": scores, "overall": ov,
        "validation": {k: val[k] for k in ("controls_summary", "controls", "permutation",
                                            "method_agreement", "signature_agreement")},
        "rat": json.loads(rat_path.read_text()) if rat_path.exists() else None,
        "traj": traj, "volcano": volcano, "map": dmap, "robust": rob, "outputs": outputs, "candidates": cands,
        "hero": hero,
    }
    blob = json.dumps(data, separators=(",", ":"), ensure_ascii=False).replace("</", "<\\/")
    tpl = (C.ROOT / "website" / "template.html").read_text(encoding="utf-8")
    html = tpl.replace("/*__DATA__*/null", blob)
    # content-only copy for publishing (the host adds doctype/head), plus a full standalone page
    (C.ROOT / "website" / "_artifact.html").write_text(html.replace("<!--BODY-->\n", ""), encoding="utf-8")
    head, body = html.split("<!--BODY-->", 1)
    standalone = ('<!doctype html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n'
                  '<meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover">\n'
                  f"{head.strip()}\n</head>\n<body>\n{body.strip()}\n</body>\n</html>\n")
    html = standalone
    (C.ROOT / "website" / "index.html").write_text(html, encoding="utf-8")
    (C.WEB_DATA / "results.json").write_text(blob, encoding="utf-8")
    print(f"website/index.html written: {len(html) / 1e6:.2f} MB, {len(drugs)} drugs, {len(sigs)} signatures")


if __name__ == "__main__":
    main()
