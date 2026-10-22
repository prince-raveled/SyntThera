# SynThera (TEJAS019)

**Mimicking exercise through molecular discovery.** Connecting molecular signatures to drug discovery.

A computational exercise-mimetic screen built directly on **human** data: MoTrPAC's randomised,
control-matched acute-exercise signatures (muscle, blood, adipose; endurance and resistance) are
matched against ~29,700 compounds in LINCS L1000 with four scoring methods, including the GSFM
gene-set foundation model, combined by Robust Rank Aggregation. A rat pipeline runs alongside as
a control arm to measure the species gap.

TEJAS 2026 AI Innovation Challenge · SASTRA Deemed University · SDG 9 · AI for Drug Discovery & Personalized Medicine

**Website:** open `website/index.html` in a browser (self-contained, works offline).

## Layout

```
pipeline/            run in order (01 → 09)
  01_build_signatures.py   18 human signatures from MoTrPAC DA tables
  02_lincs_consensus.py    reads the 33 GB L1000 matrix once; MODZ consensus per drug
  03_gene_modules.py       adaptive vs stress gene modules (Hallmark, Reactome)
  04_score.py              4 methods, per-compound null standardisation, RRA, permutation FDR, benefit index
                           (--rerank recomputes ranks from saved scores without rescoring)
  05_validate.py           positive/negative controls, permutation summary, agreement, figures
  06_robustness.py         method ablation, bootstrap stability, GSFM drug map (UMAP)
  07_rat_comparison.py     rodent control arm on identical scoring
  08_fair_export.py        TSV tables, candidate lists, schema.org/Bioschemas JSON-LD
  09_export_website.py     embeds results into website/index.html
src/tejas/           shared code (config, scoring methods, GSFM wrapper)
data/raw/            downloaded inputs (see TEJAS019_dataset_plan.md)
data/processed/      signatures, consensus matrix, caches
results/             scores_long.parquet, overall.parquet, validation.json, figures/, logs
website/             template.html + generated index.html
docs/                proposal
```

## Run

```bash
.venv/Scripts/python.exe pipeline/01_build_signatures.py
.venv/Scripts/python.exe pipeline/02_lincs_consensus.py
.venv/Scripts/python.exe pipeline/03_gene_modules.py
.venv/Scripts/python.exe pipeline/04_score.py
.venv/Scripts/python.exe pipeline/05_validate.py
.venv/Scripts/python.exe pipeline/06_robustness.py
.venv/Scripts/python.exe pipeline/07_rat_comparison.py
.venv/Scripts/python.exe pipeline/08_fair_export.py
.venv/Scripts/python.exe pipeline/09_export_website.py
```

## Candidate rule (`src/tejas/candidates.py`)

1. Approved drug (Drug Repurposing Hub, clinical phase "Launched").
2. Top 5% of exercise match across the informative signatures (RRA of per-signature ranks).
3. Ranked by mean benefit index (adaptive-module minus stress-module agreement).
4. Flagged and excluded: cytotoxic or generic-stress mechanisms (tubulin, proteasome, topoisomerase,
   cardiac glycosides, PKC activators, ...) and drugs profiled in fewer than 3 L1000 cell lines.

A plain connectivity screen on acute exercise data is dominated by cytotoxic drugs that copy the
stress of a hard bout in cancer cells; the benefit index is what separates them.

Step 2 needs the 33 GB `level5_beta_trt_cp_n720216x12328.gctx` at `E:\lincs\` (path in `src/tejas/config.py`).
Steps 2 and 4 take a few minutes and ~40 minutes respectively on a 6-core laptop.

The environment was created with `uv` (Python 3.12):

```bash
uv venv .venv --python 3.12
uv pip install --python .venv/Scripts/python.exe numpy pandas scipy h5py pyarrow matplotlib gseapy statsmodels tqdm pyreadr "gsfm @ git+https://github.com/MaayanLab/gsfm.git"
```

## Key outputs

| File | Content |
|---|---|
| `results/scores_long.parquet` | every compound × 18 signatures: raw and null-standardised method scores, ranks, RRA rho, rank, permutation p/FDR, adaptive/stress z, benefit index |
| `results/overall.parquet` | overall rank per exercise type across informative signatures |
| `results/candidates_<exercise>.tsv` | approved strong matches with benefit index, flags and candidate rank |
| `results/ranking_<exercise>.tsv` | all compounds, overall and per-signature ranks |
| `results/exercise_mimetic_atlas.jsonld` | FAIR JSON-LD (schema.org + Bioschemas) |
| `results/validation.json` | control recovery, permutation counts, method and signature agreement |
| `results/robustness.json` | method ablation, bootstrap stability, drug-map coordinates of signatures |
| `results/rat_comparison.json` | human vs rat ranking divergence, genes lost in ortholog conversion |
| `results/figures/*.png` | slide-ready figures |

## Data and licences

- MoTrPAC human pre-suspension acute exercise, Analysis collection c2.0 (CC BY 4.0). Cite: MoTrPAC Bioinformatics Center, Public Release human-precovid-sed-adu-c2.0.
- MoTrPAC rat endurance training, MotrpacRatTraining6moData (CC BY 4.0).
- CMap LINCS 2020 Level 5 compound signatures and Drug Repurposing Hub 2020-03-24 (CLUE terms of use).
- MSigDB Hallmark 2020 and Reactome 2022 gene sets via Enrichr.
- GSFM: Clarke DJB, Marino GB, Ma'ayan A. *Patterns* 2026 (BSD-3-Clause).

Research prototype for hypothesis generation; not medical advice.
