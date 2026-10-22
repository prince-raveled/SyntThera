"""Step 8: FAIR outputs.

* results/ranking_<modality>.tsv          every compound, overall and per-signature ranks
* results/top25_approved_annotated.tsv     top approved drugs per exercise type with annotation
* results/exercise_mimetic_atlas.jsonld    schema.org / Bioschemas JSON-LD: provenance (isBasedOn),
  licences, variables, ranked compounds (InChIKey, Broad IDs, PubChem link) and the query gene sets
  (NCBI Gene identifiers)
"""
import json
import sys
from datetime import date
from pathlib import Path
from urllib.parse import quote

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from tejas import config as C  # noqa: E402
from tejas.candidates import MATCH_TOP, MIN_CELLS, candidate_table  # noqa: E402

TEAM = ["Jahnavi Sai G", "Prince Kumar", "Apurva V", "Priyanka G", "Vishnupriya K"]
PRIMARY = {"endurance": "muscle|endurance|3.5-4h", "resistance": "muscle|resistance|3.5-4h"}


def main():
    res = pd.read_parquet(C.RESULTS / "scores_long.parquet")
    cpd = pd.read_parquet(C.PROC / "compounds.parquet")
    overall = pd.read_parquet(C.RESULTS / "overall.parquet")
    meta = json.loads((C.RESULTS / "scoring_meta.json").read_text())
    queries = json.loads((C.RESULTS / "query_sets.json").read_text())
    sigl = pd.read_parquet(C.PROC / "signatures_long.parquet")
    informative = meta["informative_signatures"]

    ci = pd.read_csv(C.LINCS / "compoundinfo_beta.txt", sep="\t")
    ci = ci.groupby("cmap_name").agg(inchi_key=("inchi_key", "first"), smiles=("canonical_smiles", "first")).reset_index()
    cpd = cpd.merge(ci, left_on="drug_id", right_on="cmap_name", how="left").drop(columns=["cmap_name"])
    base_cols = ["drug_id", "name", "clinical_phase", "moa", "target", "n_cells", "n_sigs", "pert_ids", "inchi_key"]

    tables = {}
    for mod in C.MODALITIES:
        o = overall[overall.modality == mod][["drug_id", "rank", "rho", "median_rank"]] \
            .rename(columns={"rank": "overall_rank", "rho": "overall_rra_rho", "median_rank": "median_signature_rank"})
        wide = res[(res.modality == mod) & res.signature.isin(informative)] \
            .pivot(index="drug_id", columns="signature", values="rank")
        wide.columns = [f"rank_{c.replace('|', '_')}" for c in wide.columns]
        prim = res[res.signature == PRIMARY[mod]].set_index("drug_id")[["perm_fdr", "benefit", "adaptive_z", "stress_z"]]
        prim.columns = [f"primary_{c}" for c in prim.columns]
        t = cpd[base_cols].merge(o, on="drug_id").merge(prim, left_on="drug_id", right_index=True) \
            .merge(wide, left_on="drug_id", right_index=True).sort_values("overall_rank")
        t.to_csv(C.RESULTS / f"ranking_{mod}.tsv", sep="\t", index=False)
        tables[mod] = t
        print(f"ranking_{mod}.tsv: {len(t):,} compounds")

    cands = {}
    for mod in C.MODALITIES:
        ct = candidate_table(res, cpd.drop(columns=["inchi_key", "smiles"]), overall, informative, mod)
        ct = ct.merge(cpd[["drug_id", "inchi_key", "pert_ids"]], on="drug_id", how="left")
        ct[ct.approved & ct.strong_match].to_csv(C.RESULTS / f"candidates_{mod}.tsv", sep="\t", index=False)
        cands[mod] = ct
        n_flag = int((ct.approved & ct.strong_match & ct.flag_cytotoxic).sum())
        print(f"{mod}: {int(ct.candidate.sum())} candidates; {n_flag} strong-match approved drugs flagged cytotoxic; "
              f"top 5: {', '.join(ct[ct.candidate].head(5).name)}")
    top = pd.concat([cands[m][cands[m].candidate].head(25) for m in C.MODALITIES])
    top.to_csv(C.RESULTS / "top25_candidates_annotated.tsv", sep="\t", index=False)

    def compound_node(row, mod):
        props = [("exerciseType", mod), ("exerciseMatchRank", int(row.match_rank)),
                 ("exerciseMatchRRArho", float(f"{row.match_rho:.3g}")),
                 ("medianSignatureRank", int(row.median_rank)),
                 ("meanBenefitIndex", round(float(row.mean_benefit), 3)),
                 ("meanAdaptiveZ", round(float(row.mean_adaptive_z), 3)),
                 ("meanStressZ", round(float(row.mean_stress_z), 3)),
                 ("informativeSignaturesFDR05", int(row.n_signatures_fdr05)),
                 ("bestSignature", row.best_signature),
                 ("clinicalPhase", row.clinical_phase), ("l1000CellLines", int(row.n_cells)),
                 ("flagCytotoxicMechanism", bool(row.flag_cytotoxic)), ("flagLowEvidence", bool(row.flag_low_evidence))]
        if row.candidate:
            props.insert(1, ("candidateRank", int(row.candidate_rank)))
        node = {"@type": "MolecularEntity", "name": row["name"],
                "identifier": [f"broad:{p.strip()}" for p in str(row.pert_ids).split(";") if p.strip()],
                "additionalProperty": [{"@type": "PropertyValue", "name": k, "value": v} for k, v in props]}
        if isinstance(row.inchi_key, str) and len(row.inchi_key) == 27:
            node["inChIKey"] = row.inchi_key
            node["sameAs"] = f"https://pubchem.ncbi.nlm.nih.gov/#query={quote(row.inchi_key)}"
        if isinstance(row.moa, str):
            node["description"] = f"Mechanism: {row.moa}"
        return node

    entrez = dict(zip(sigl.gene_symbol.astype(str), sigl.entrez.astype(str)))

    def gene_set_node(sid):
        q = queries[sid]
        return {"@type": "Dataset", "name": f"Query gene set {sid.replace('|', ' / ')}",
                "description": ("Top FDR-significant genes (difference-in-changes vs control). "
                                + ("Exploratory: top-|t| fallback." if q["fallback"] else "")).strip(),
                "about": [{"@type": "Gene", "name": g, "identifier": f"ncbigene:{entrez.get(g, '')}", "additionalType": d}
                          for d in ("up", "down") for g in q["up" if d == "up" else "down"]]}

    variables = [
        ("rra_rho", "Robust Rank Aggregation rho over four connectivity methods (lower = stronger exercise match)"),
        ("rank", "Rank among all screened compounds for one signature"),
        ("perm_fdr", "Benjamini-Hochberg FDR from a pooled null of random-query permutations"),
        ("benefit", "Adaptive-module agreement minus stress-module agreement (z units)"),
        ("wtcs", "CMap weighted connectivity score"), ("spearman", "Spearman correlation, landmark genes"),
        ("xsum", "Extreme-sum score"), ("gsfm", "GSFM embedding cosine score"),
    ]
    doc = {
        "@context": {"@vocab": "https://schema.org/", "bioschemas": "https://bioschemas.org/",
                     "broad": "https://clue.io/command?q=", "ncbigene": "https://www.ncbi.nlm.nih.gov/gene/"},
        "@type": "Dataset",
        "@id": "#exercise-mimetic-atlas",
        "name": "Exercise Mimetic Atlas: human-signature exercise-mimetic drug rankings",
        "description": ("Rankings of LINCS L1000 compounds by similarity to human acute-exercise transcriptomic "
                        "signatures (MoTrPAC, 175 sedentary adults; muscle, blood, adipose; endurance and resistance; "
                        "difference-in-changes vs non-exercise control). Four connectivity methods including the GSFM "
                        "gene-set foundation model, combined with Robust Rank Aggregation; permutation FDR; "
                        "adaptive-vs-stress benefit index; rat comparison arm."),
        "version": "1.0.0", "dateCreated": date.today().isoformat(),
        "license": "https://creativecommons.org/licenses/by/4.0/",
        "creator": [{"@type": "Person", "name": n, "affiliation": {"@type": "Organization", "name": "SASTRA Deemed University"}} for n in TEAM],
        "keywords": ["exercise mimetics", "drug repurposing", "connectivity map", "LINCS L1000", "MoTrPAC",
                     "gene set foundation model", "robust rank aggregation", "FAIR"],
        "measurementTechnique": ["connectivity mapping", "robust rank aggregation", "permutation testing"],
        "isBasedOn": [
            {"@type": "Dataset", "name": f"MoTrPAC {C.MOTRPAC_RELEASE} (human pre-suspension acute exercise, analysis collection)",
             "url": "https://motrpac-data.org/data-download", "license": "https://creativecommons.org/licenses/by/4.0/"},
            {"@type": "Dataset", "name": C.LINCS_RELEASE, "url": "https://clue.io/releases/data-dashboard",
             "citation": "https://doi.org/10.1016/j.cell.2017.10.049"},
            {"@type": "Dataset", "name": C.HUB_RELEASE, "url": "https://clue.io/repurposing"},
            {"@type": "SoftwareApplication", "name": "GSFM", "url": "https://github.com/MaayanLab/gsfm",
             "identifier": "https://doi.org/10.1016/j.patter.2026.101565", "softwareVersion": C.GSFM_MODEL},
            {"@type": "Dataset", "name": "MoTrPAC rat endurance training (MotrpacRatTraining6moData)",
             "url": "https://github.com/MoTrPAC/MotrpacRatTraining6moData",
             "citation": "https://doi.org/10.1038/s41586-023-06877-w"},
            {"@type": "Dataset", "name": "MSigDB Hallmark 2020 and Reactome 2022 via Enrichr", "url": "https://maayanlab.cloud/Enrichr/"},
        ],
        "variableMeasured": [{"@type": "PropertyValue", "name": k, "description": v} for k, v in variables],
        "hasPart": (
            [{"@type": "Dataset", "name": f"Approved {m} exercise-mimetic candidates and flagged strong matches",
              "description": (f"Approved drugs in the top {MATCH_TOP:.0%} of exercise match (RRA across informative "
                              f"signatures). Candidates exclude cytotoxic/generic-stress mechanisms and drugs profiled in "
                              f"fewer than {MIN_CELLS} cell lines, and are ranked by mean benefit index."),
              "hasPart": [compound_node(r, m) for _, r in cands[m][cands[m].approved & cands[m].strong_match].iterrows()]}
             for m in C.MODALITIES]
            + [gene_set_node(s) for s in informative]
            + [{"@type": "DataDownload", "name": f.name, "contentUrl": f"results/{f.name}", "encodingFormat": "text/tab-separated-values"}
               for f in sorted(C.RESULTS.glob("*.tsv"))]
        ),
    }
    (C.RESULTS / "exercise_mimetic_atlas.jsonld").write_text(json.dumps(doc, indent=1, ensure_ascii=False), encoding="utf-8")
    n_ik = int(top.inchi_key.notna().sum())
    print(f"JSON-LD written; top-25 approved table has InChIKeys for {n_ik}/{len(top)} drugs")


if __name__ == "__main__":
    main()
