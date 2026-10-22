"""Step 3: define adaptive vs stress gene modules (MSigDB Hallmark 2020 + Reactome 2022 via Enrichr).

Adaptive = adaptations worth copying (mitochondria, oxidative metabolism, muscle, vessels).
Stress   = inflammation / damage / cell-stress programmes that accompany a hard bout.
Genes in both, plus immediate-early genes (FOS, JUN, EGR, NR4A, ATF3...), are 'ambiguous'
and are left out of both scores: NR4A3 is the most exercise-responsive gene in MetaMEx and
is a mediator of adaptation, not damage.
"""
import json
import re
import sys
from pathlib import Path

import gseapy
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from tejas import config as C  # noqa: E402

LIBRARIES = ["MSigDB_Hallmark_2020", "Reactome_2022"]

ADAPTIVE = {
    "Oxidative phosphorylation": r"^Oxidative Phosphorylation$",
    "Fatty acid metabolism": r"^Fatty Acid Metabolism$",
    "Peroxisome": r"^P?peroxisome$",   # Enrichr's Hallmark 2020 spells it "Pperoxisome"
    "Myogenesis": r"^Myogenesis$",
    "Angiogenesis": r"^Angiogenesis$",
    "mTORC1 signalling": r"^mTORC1 Signaling$",
    "Mitochondrial biogenesis": r"^Mitochondrial Biogenesis R-HSA-1592230$",
    "PGC-1a mitochondrial programme": r"^Transcriptional Activation Of Mitochondrial Biogenesis R-HSA-2151201$",
}
STRESS = {
    "TNFa / NF-kB": r"^TNF-alpha Signaling via NF-kB$",
    "Inflammatory response": r"^Inflammatory Response$",
    "IL-6 / JAK / STAT3": r"^IL-6/JAK/STAT3 Signaling$",
    "Interferon-gamma response": r"^Interferon Gamma Response$",
    "Complement": r"^Complement$",
    "Apoptosis": r"^Apoptosis$",
    "p53 pathway": r"^p53 Pathway$",
    "Unfolded protein response": r"^Unfolded Protein Response$",
}
AMBIGUOUS_SETS = {
    "Hypoxia": r"^Hypoxia$",
    "Reactive oxygen species": r"^Reactive Oxygen Species Pathway$",
}
IMMEDIATE_EARLY = ["FOS", "FOSB", "FOSL1", "FOSL2", "JUN", "JUNB", "JUND", "EGR1", "EGR2", "EGR3",
                   "ATF3", "NR4A1", "NR4A2", "NR4A3", "IER2", "IER3", "IER5", "MYC", "DUSP1"]


def fetch(lib: str) -> dict:
    path = C.GENESETS / f"{lib}.json"
    if path.exists():
        return json.loads(path.read_text())
    gs = gseapy.get_library(name=lib, organism="Human")
    path.write_text(json.dumps(gs))
    return gs


def pick(libs: dict, patterns: dict) -> dict:
    out = {}
    for label, pat in patterns.items():
        hits = [(lib, name) for lib, gs in libs.items() for name in gs if re.match(pat, name, re.I)]
        if len(hits) != 1:
            raise SystemExit(f"module '{label}': expected 1 gene set for /{pat}/, found {hits}")
        lib, name = hits[0]
        out[label] = {"source": f"{lib}: {name}", "genes": sorted(set(libs[lib][name]))}
    return out


def main():
    libs = {lib: fetch(lib) for lib in LIBRARIES}
    for lib, gs in libs.items():
        print(f"{lib}: {len(gs)} gene sets")
    adaptive, stress, amb = pick(libs, ADAPTIVE), pick(libs, STRESS), pick(libs, AMBIGUOUS_SETS)

    a_genes = {g for m in adaptive.values() for g in m["genes"]}
    s_genes = {g for m in stress.values() for g in m["genes"]}
    amb_genes = {g for m in amb.values() for g in m["genes"]} | (a_genes & s_genes) | set(IMMEDIATE_EARLY)

    rows = []
    for cls, mods in (("adaptive", adaptive), ("stress", stress), ("ambiguous", amb)):
        for label, m in mods.items():
            for g in m["genes"]:
                rows.append({"gene_symbol": g, "module": label, "module_class": cls})
    for g in IMMEDIATE_EARLY:
        rows.append({"gene_symbol": g, "module": "Immediate-early genes", "module_class": "ambiguous"})
    long = pd.DataFrame(rows).drop_duplicates()
    long.to_parquet(C.PROC / "gene_modules_long.parquet", index=False)

    gene_class = {g: "adaptive" for g in a_genes - amb_genes}
    gene_class.update({g: "stress" for g in s_genes - amb_genes})
    gene_class.update({g: "ambiguous" for g in amb_genes})
    pd.DataFrame({"gene_symbol": list(gene_class), "gene_class": list(gene_class.values())}) \
        .to_parquet(C.PROC / "gene_class.parquet", index=False)

    meta = {"adaptive": {k: v["source"] for k, v in adaptive.items()},
            "stress": {k: v["source"] for k, v in stress.items()},
            "ambiguous": {**{k: v["source"] for k, v in amb.items()},
                          "Immediate-early genes": "curated: " + ", ".join(IMMEDIATE_EARLY)}}
    (C.PROC / "gene_modules_meta.json").write_text(json.dumps(meta, indent=2))
    counts = pd.Series(gene_class).value_counts()
    print("genes per class (after removing overlaps):")
    print(counts.to_string())


if __name__ == "__main__":
    main()
