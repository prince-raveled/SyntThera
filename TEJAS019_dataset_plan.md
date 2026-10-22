# TEJAS019: Datasets and Build Plan

**Project:** Human-signature exercise-mimetic drug screen (MoTrPAC human data × LINCS L1000)
**Sub-theme:** AI for Drug Discovery & Personalized Medicine · SDG 9
**Plan date:** 12 Oct 2026 · Data checked live on motrpac-data.org the same day

---

## 1. TL;DR: what you need to start

| # | Dataset | Role in the project | Access | Size |
|---|---------|---------------------|--------|------|
| 1 | **MoTrPAC Human Acute Exercise (pre-suspension), Analysis c2.0**: RNA-seq differential-analysis tables for muscle, adipose and blood | **The exercise signature** (core of the idea) | Public, CC BY 4.0, no application | ~605 MB |
| 2 | **LINCS L1000 / CMap 2020 release**: Level 5 compound signatures + metadata | The drug library to screen against | Free (clue.io) | Very large matrix; slice it, don't load it whole |
| 3 | **Broad Drug Repurposing Hub** annotations | Keep only approved ("Launched") drugs; add MOA/targets | Free | < 10 MB |
| 4 | **MSigDB Hallmark + Reactome gene sets** | Split signature into *adaptive* vs *stress* modules | Free | < 50 MB |
| 5 | **GSFM pretrained model** (Ma'ayan Lab) | AI embedding-based scoring | Free, BSD-3 | Model download via `pip` |
| 6 | **MetaMEx** (human muscle exercise meta-analysis) | Independent human validation + training (chronic) signatures | Free | Small |
| 7 | *(Control arm only)* **MoTrPAC Rat Endurance Training, Analysis c2.0** + rat→human ortholog table | Rodent comparison pipeline | Public, CC BY 4.0 | ~hundreds of MB |

**Day 1:** download #1 (Section 4 has the exact clicks) and the four L1000 metadata files from #2. That is enough to build the human signature and start matching genes to L1000.

---

## 2. What the human MoTrPAC data really is (read before building)

These facts come from the live Data Hub and its release notes. Several of them should change wording in the proposal.

### 2.1 It is an ACUTE (single-bout) dataset, not a training dataset
- The human study randomized **206** sedentary adults (EE 80, RE 81, CON 45). COVID-19 suspended the trial. **176** started the baseline acute exercise test, and **175** are in the analysis release. Your "175" is correct.
- Only **55** finished the 12-week training, so the **public release covers the single acute bout only**.
- The website's four pillars, mapped to what you can download for **humans** today:

| Website pillar | Human data available now? |
|---|---|
| Multi-species: **Humans** / Rodents | ✅ Humans (your focus). Rodents only as the control arm |
| Multi-modal: **Endurance** / **Resistance** | ✅ Both (EE-vs-CON and RE-vs-CON contrasts) |
| Multi-tissue: **Blood / Skeletal muscle / Adipose** | ✅ All three |
| Multi-temporal: **24-h acute profiling** | ✅ Yes, this is the release |
| Multi-temporal: **Weeks of exercise training** | ❌ **Not for humans.** Public training data exist only for rats |

**Consequence:** your human signature is the *response to one bout*, not *training adaptation*. The proposal's "beneficial adaptation vs. unwanted stress" idea still works, but:
- treat the **24 h** timepoint as closest to "adaptation" and the early timepoints as mostly "stress / immediate response";
- use **MetaMEx human training signatures** (dataset 6) as the reference for *chronic adaptation*.

### 2.2 Only summary-level results are open, and that is all you need
| Level | Content | Access |
|---|---|---|
| Level 3 **DA tables** | logFC, t, p, FDR per gene × timepoint × contrast | **Public**, CC BY 4.0 |
| Level 4 enrichment (CAMERA pathway results, clusters) | Pre-computed pathway enrichment | Public, but **only via the R package**, not as Data Hub downloads |
| Level 2 per-sample normalized data + phenotypes | Individual-level data | **Restricted**: dbGaP `phs004523.v1.p1` + committee approval |

You do **not** need dbGaP. Drug-signature matching runs on DA tables.

### 2.3 How the DA tables were made (state this in methods)
- Model: `dream` (variancePartition) linear mixed model; `voomWithDreamWeights` for RNA-seq.
- **Primary contrast = "difference-in-changes":** (exercise group post − pre) minus (control group post − pre). It removes time-of-day, fasting and repeated-biopsy effects, which makes it a clean signature.
- Covariates: clinical site, age, sex, BMI; random effect = participant. BH-FDR 0.05.
- Contrast categories in the files/package: `EE-CON`, `RE-CON` (use these), plus `EE-EE`, `RE-RE`, `EE-RE`, `CON-CON`.

### 2.4 Sampling timepoints (these define your signatures)
| Tissue | Timepoints | RNA-seq genes / samples |
|---|---|---|
| **Skeletal muscle** (vastus lateralis), primary | pre, **15 min**, **3.5 h**, **24 h** post | 16,509 / 485 |
| Adipose | pre, 45 min, 4 h, 24 h post | 20,090 / 406 |
| Blood | EE/CON: pre, during 20 & 40 min, 10 min, 30 min, 3.5 h, 24 h. RE: no during-exercise draws | 20,526 / 951 |

Timepoint labels in the data: `pre_exercise, during_20_min, during_40_min, post_10_min, post_15_30_45_min, post_3.5_4_hr, post_24_hr`.

### 2.5 Use version c2.0 only
c2.0 (released **30 Sep 2026**) fixes a bug in c1.3 where transcriptomics batch correction was applied to the **wrong samples**. Any c1.3 transcriptomic result is superseded. Record "MoTrPAC human-precovid-sed-adu **c2.0**" in every output.

### 2.6 Download early
The Data Hub footer says: *"This repository is under review for potential modification in compliance with Administration directives."* Download now, keep the files, and record checksums (`sha256sum`).

---

## 3. Dataset details

### 3.1 MoTrPAC human DA files (dataset #1)
From **Analysis c2.0 → Transcriptomics** (exact names and sizes from the file browser):

| File | Size | Need? |
|---|---|---|
| `human-precovid-sed-adu_t06-muscle_transcript-rna-seq_da_dream-acute_v2.1.txt` | 144.63 MB | **Must** |
| `human-precovid-sed-adu_t06-muscle_transcript-rna-seq_metadata_features_v2.0.txt` | 1.03 MB | **Must** (gene IDs) |
| `human-precovid-sed-adu_t11-adipose_transcript-rna-seq_da_dream-acute_v2.1.txt` | 175.47 MB | Should |
| `human-precovid-sed-adu_t11-adipose_transcript-rna-seq_metadata_features_v2.0.txt` | 1.25 MB | Should |
| `human-precovid-sed-adu_t04-blood-rna_transcript-rna-seq_da_dream-acute_v2.1.txt` | 281.41 MB | Should |
| `human-precovid-sed-adu_t04-blood-rna_transcript-rna-seq_metadata_features_v2.0.txt` | 1.27 MB | Should |
| `human-precovid-sed-adu_all_transcript-rna-seq_qc-report_v1.2.html` | 8.05 MB | Nice to read |
| `motrpac-mappings-human-feature-to-gene_v2.2.txt` | 182.72 MB | Optional (in c2.0 each `metadata_features` file is already a self-contained gene map) |

Optional later: muscle proteomics DA (`t10-muscle`, `prot-pr`) to check that top genes also move at the protein level.

**DA columns** (per the R package docs): `feature_id, logFC, t, p_value, adj_p_value, AveExpr, z.std, CI.L_calculated, CI.R_calculated`, plus contrast/timepoint fields. Inspect the header first:

```python
import pandas as pd
da = pd.read_csv("human-precovid-sed-adu_t06-muscle_transcript-rna-seq_da_dream-acute_v2.1.txt", sep="\t")
print(da.columns.tolist()); print(da.head())
```

**R alternative** (R ≥ 4.4). This route also gives you the Level-4 CAMERA pathway results:
```r
BiocManager::install(version = "3.20")          # R 4.4 → 3.20, R 4.5 → 3.22, R 4.6 → 3.23
pak::pak("MoTrPAC/MotrpacHumanPreSuspensionAnalysis")
library(MotrpacHumanPreSuspensionAnalysis)
ome_available_list(); tissue_available_list()   # check exact argument values
da <- load_differential_analysis(selected_omes = "transcript-rna-seq", selected_tissues = "muscle")
CAMERA_RESULTS                                   # pre-computed pathway enrichment
```

### 3.2 LINCS L1000 (dataset #2)
- Source: CMap/LINCS **2020 release** via the clue.io data dashboard (`https://clue.io/releases/data-dashboard`). Many projects pull the same files from `https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/`; check that it still works.
- Files:
  - `level5_beta_trt_cp_n720216x12328.gctx`: 720,216 compound signatures × 12,328 genes (~34k compounds). Read **slices only** with `cmapPy` (`parse(..., cid=..., rid=...)`); never load the whole file on Colab.
  - `siginfo_beta.txt`: signature → compound, dose, time, cell line, quality flags
  - `compoundinfo_beta.txt`: compound names, MOA, targets
  - `geneinfo_beta.txt`: gene IDs, **landmark vs inferred** flag
  - `cellinfo_beta.txt`: cell line lineage
- **No-download alternative for the matrix:** the LINCS README says the same data are on **Google BigQuery** (project `cmap-big-table`, dataset `cmap_lincs_public_views`), where you can query just the signatures you need. The free BigQuery sandbox needs a Google account.
- **Faster first pass (no big download):** SigCom LINCS (`maayanlab.cloud/sigcom-lincs`) accepts up/down gene sets and returns ranked drugs. Use it to sanity-check results on day 2–3.
- Caveats to handle:
  - Only **978 landmark genes** are measured; the rest are inferred. Do correlation scoring on landmark genes.
  - Cell lines are mostly **cancer lines**. Check `cellinfo_beta.txt` for any muscle-lineage cells and expect few or none. Prefer compounds with consistent signatures across cell lines.
  - Filter to reproducible signatures (use the quality fields in `siginfo_beta.txt`).

### 3.3 Drug Repurposing Hub (dataset #3)
Broad Institute (`clue.io/repurposing`). Join L1000 `pert_id` (e.g. `BRD-K12345678`) to the first 13 characters of the Hub's `broad_id`. Keep `clinical_phase == "Launched"` for the "approved drugs" claim.

### 3.4 Gene sets (dataset #4), for the adaptive vs stress split
Get them from MSigDB through `msigdbr` (R) or `gseapy` (Python). Proposed starting split; the team should finalize and justify it:

| Adaptive (worth copying) | Stress / damage (avoid) | Ambiguous (report separately) |
|---|---|---|
| OXIDATIVE_PHOSPHORYLATION, FATTY_ACID_METABOLISM, PEROXISOME, MYOGENESIS, ANGIOGENESIS, MTORC1_SIGNALING (resistance), mitochondrial biogenesis (Reactome/GO) | TNFA_SIGNALING_VIA_NFKB, INFLAMMATORY_RESPONSE, IL6_JAK_STAT3_SIGNALING, INTERFERON_GAMMA_RESPONSE, COMPLEMENT, APOPTOSIS, P53_PATHWAY, UNFOLDED_PROTEIN_RESPONSE | HYPOXIA, REACTIVE_OXYGEN_SPECIES, immediate-early genes (FOS, JUN, EGR1, ATF3, **NR4A3**) |

NR4A3 was the most exercise-responsive gene in MetaMEx, so don't throw immediate-early genes into "stress" automatically. Also cross-check against the MoTrPAC `CAMERA_RESULTS`.

### 3.5 GSFM (dataset #5)
```bash
pip install "gsfm@git+https://github.com/MaayanLab/gsfm.git"
```
```python
import torch
from gsfm import Vocab, GSFM   # see repo README for exact import path
vocab = Vocab.from_pretrained("maayanlab/gsfm-rummagene")   # also: gsfm-rummageo, gsfm-rummage
model = GSFM.from_pretrained("maayanlab/gsfm-rummagene").eval()
emb = model.encode(torch.tensor(vocab(gene_symbols))[None, :])
```
Input is **human gene symbols**, so the human-first design fits it directly. Citation: Clarke, Marino & Ma'ayan, *Patterns* 2026.

### 3.6 MetaMEx (dataset #6)
Pillon et al., *Nat Commun* 2020. 66 human muscle datasets covering acute aerobic, acute resistance, aerobic training, resistance training and **inactivity**. Get it from `metamex.eu` or GitHub `NicoPillon/MetaMEx` (`/data`). Uses:
- **Replication:** drug rankings from MoTrPAC muscle should correlate with rankings from MetaMEx acute signatures.
- **Chronic adaptation reference:** MetaMEx training signatures fill the gap left by the missing human training data (§2.1).
- **Bonus angle:** drugs that *reverse* the inactivity signature are another way to define a mimetic.

### 3.7 Rodent control arm (dataset #7, secondary)
Rat endurance training, **Analysis c2.0** (public). Exact names:
- `pass1b-06_t55-gastrocnemius_transcript-rna-seq_timewise-dea_v2.0.txt`
- `pass1b-06_t56-vastus-lateralis_transcript-rna-seq_timewise-dea_v2.0.txt`
- `pass1b-06_t70-white-adipose_transcript-rna-seq_timewise-dea_v2.0.txt`
- `pass1b-06_t30-blood-rna_transcript-rna-seq_timewise-dea_v2.0.txt`

Notes:
- `timewise-dea` is per sex and per week (1, 2, 4, 8 wk vs sedentary). `training-dea` is an overall test with no direction, so it can't serve as a signature.
- Ortholog mapping: HGNC HCOP, Ensembl BioMart, or the `babelgene` R package. Report how many genes are lost; that number is a result.
- ⚠ **Confound:** human data are *acute*, rat data are *8-week training*. A human-vs-rat ranking difference therefore mixes **species** with **acute vs chronic**. Either compare rat training with **MetaMEx human training**, or state this limitation explicitly.
- ⚠ **Do not use** "Acute Exercise in Young Adult Rats (Limited)". Its data use agreement forbids using the data in abstracts or presentations during an embargo, and a hackathon pitch counts.

---

## 4. Step-by-step on the website

The page you pasted (motrpac.org) is the consortium homepage. Its **"Access Data"** link goes to the **Data Hub: https://motrpac-data.org**.

1. Open **motrpac-data.org** → **Browse Data** → **Download Datasets**.
2. On the card **"Acute Exercise in Human Sedentary Adults (Pre-Suspension)"**, the **Analysis c2.0 · Public Release** row → **Browse Files**.
   (The Quant-ID and Phenotype rows say *Restricted, apply through dbGaP*. Skip them.)
3. Left filter → **Omics: Transcriptomics**. Ten files appear.
4. Tick the 3 `_da_dream-acute_v2.1.txt` files, the 3 `_metadata_features_v2.0.txt` files, and the QC report → **Download selected files**. If it asks you to sign in, use **Login** (top right). The download step itself was not tested in this check.
5. Read the **Data Release Notes** PDF (link on the Download page) and the Knowledge Center page
   `motrpac-data.org/knowledge-center/data-collections/analysis-results/analysis-human-precovid-sed-adu`.
6. Explore genes visually at **data-viz.motrpac-data.org/precawg** to see what the signature looks like before coding.
7. (Control arm) Same steps on the **"Endurance Training in Young Adult Rats"** card → Analysis c2.0 → Transcriptomics → the 4 `timewise-dea` files in §3.7.
8. **Cite:** "MoTrPAC Bioinformatics Center (BIC). Public Release: human-precovid-sed-adu-c2.0", license CC BY 4.0. Notify MoTrPAC-ACC@aging.ufl.edu if you publish.

---

## 5. Build plan (pipeline)

### Phase 0: Setup
- Google Colab / Kaggle (free). Python: `pandas, numpy, scipy, cmapPy, gseapy, torch, gsfm`. R (optional): `MotrpacHumanPreSuspensionAnalysis, msigdbr, RobustRankAggreg`.
- Repo layout:
```
D:\Tejas\                         (created 12 Oct 2026)
├── TEJAS019_dataset_plan.md
├── docs\                         # proposal (TEJAS019 final.docx / .pdf)
├── data\raw\motrpac_human\       # the 7 MoTrPAC transcriptomics files
├── data\raw\lincs\               # *_beta.txt metadata + README_LINCS2020.txt
├── data\raw\repurposing_hub\     # Drug Repurposing Hub drug/sample files
├── data\raw\genesets\  data\raw\metamex\  data\raw\motrpac_rat\
├── data\processed\               # signatures, L1000 slices
├── notebooks\                    # 01_signature  02_lincs_slice  03_scoring  04_validation
├── src\                          # reusable functions
└── results\                      # ranked drug tables, figures
E:\lincs\                         # the 33 GB level5_beta_trt_cp gctx (too big for D:)
```

### Phase 1: Human exercise signatures
1. Load muscle DA → keep the `EE-CON` and `RE-CON` contrasts (difference-in-changes).
2. One signature per **tissue × modality × timepoint**. **Primary: muscle, 3.5 h.** Sensitivity: 15 min and 24 h.
3. Signed score per gene = `t` (or `z.std`). Map `feature_id` (Ensembl) → Entrez/symbol via `metadata_features` → intersect with L1000 `geneinfo_beta`.
4. Gene-set form: top ~100–150 up and ~100–150 down by `t` among FDR < 0.05 genes.

### Phase 2: Adaptive vs stress split
Assign signature genes to modules (§3.4) → build an **adaptive sub-signature** and a **stress sub-signature** for each primary signature.

### Phase 3: Score ~34k compounds with four methods (as promised in the proposal)
| # | Method | Input |
|---|---|---|
| 1 | CMap weighted connectivity score (WTCS, KS-based) | up/down sets vs. drug ranked list |
| 2 | Spearman correlation | full signed vector, landmark genes |
| 3 | XSum (extreme-sum) | top-N drug genes vs. signature |
| 4 | **GSFM embedding similarity**: cos(ex_up, drug_up) + cos(ex_down, drug_down) − cross terms | gene sets → GSFM embeddings |

- Combine with **Robust Rank Aggregation**: R `RobustRankAggreg`, or a ~30-line Python implementation.
- **Benefit index** = score vs adaptive sub-signature − score vs stress sub-signature (or the ratio, as written in the proposal).
- **Score ALL ~34k compounds first, validate (Phase 4), and only then filter to Launched drugs.** AICAR (`acadesine`, Phase 3), GW-501516 (Phase 2) and torin-1/2 (Preclinical) are not Launched, so filtering early would remove the key positive and negative controls.
- Final list = Launched drugs (Repurposing Hub). Checked 12 Oct 2026: 2,427 Launched drugs; **1,404** have L1000 signatures (185,089 sigs), **1,076** have ≥1 high-quality `is_hiq` signature (28,367 sigs). Join: first 13 characters of the Hub's `broad_id` = L1000 `pert_id`.
- Report endurance and resistance lists separately.

### Phase 4: Validation (needed before trusting any hit)
- **Permutation:** 1,000 random gene sets of equal size → empirical p for every drug.
- **Positive controls:** AICAR, metformin, GW501516 (PPARδ), resveratrol; for resistance, β2-agonists (formoterol, clenbuterol, salbutamol). Metric: percentile rank / AUROC vs. all compounds.
  - Checked in `compoundinfo_beta.txt` (12 Oct 2026): AICAR is listed as **`acadesine`** (BRD-A81585303, BRD-A95696820); metformin BRD-K79602928; GW-501516 BRD-K14880289; resveratrol, formoterol, clenbuterol and salbutamol are present. **SR9009 is not in L1000.** Drop it.
  - Metformin may give weak L1000 signatures in cancer cell lines at the doses tested. If it ranks low, check that before blaming the pipeline.
- **Negative control:** mTOR inhibitors (`sirolimus` = rapamycin, `everolimus`, `torin-1`, `torin-2`, all present in L1000) are expected to score *opposite* to the resistance signature, because resistance exercise activates mTORC1.
- **Replication:** rerun with MetaMEx human signatures; the top-k overlap should beat chance.
- **Cross-tissue:** muscle vs adipose vs blood rankings. Hits that agree across tissues are stronger.

### Phase 5: Rodent comparison arm
Same pipeline on rat signatures (ortholog-mapped). Report Spearman correlation of rankings, rank-biased overlap of the top 50/100, genes lost in mapping, and positive-control recovery in each arm. State the acute-vs-chronic confound (§3.7).

### Phase 6: Outputs (FAIR)
- Ranked candidate tables (CSV) + code on GitHub, with a frozen release on Zenodo (DOI).
- Record all versions: MoTrPAC c2.0, LINCS 2020, MSigDB version, GSFM model ID.
- **Prototype demo:** a small Streamlit/Gradio app on Hugging Face Spaces. Pick tissue / modality / timepoint → ranked drugs with adaptive-vs-stress scores and a human-vs-rat rank plot.

---

## 6. Suggested wording fixes for the proposal

| Current text | Issue | Safer wording |
|---|---|---|
| "175 sedentary adults randomised to endurance, resistance or control" | Correct, but implies training | "...profiled across a single acute exercise bout (training data not yet public)" |
| "comparable human data did not exist" | MetaMEx (2020; 66 datasets, >1,100 people) exists | "no large, randomised, control-matched, multi-tissue human dataset existed" |
| "Every published computational screen derives its exercise signature from rats or mice" | No systematic review done yet (one 2025–26 mouse CMap study found, consistent with the claim) | "To our knowledge, published screens rely on rodent signatures" plus a quick literature check |
| "beneficial adaptation" | Acute data show early responses, not adaptation | "early adaptive signalling (24 h) vs acute stress", with MetaMEx training as the adaptation reference |
| "first large public human exercise dataset was released in 2026" | Fine: c1.3 Jan 2026, c2.0 30 Sep 2026 | Cite "MoTrPAC c2.0 (Sep 2026)" |

---

## 7. Sources
- MoTrPAC Data Hub, Download Datasets: https://motrpac-data.org/data-download
- Release notes (human pre-suspension, c2.0): https://motrpac-data.org/knowledge-center/data-collections/analysis-results/analysis-human-precovid-sed-adu
- Human preprint: *MoTrPAC: Initial Insights into the Dynamic Human Responses to Exercise* (bioRxiv 2026): https://www.biorxiv.org/content/10.64898/2026.03.02.705347v1.full
- R package: https://github.com/MoTrPAC/MotrpacHumanPreSuspensionAnalysis
- Rat training data: https://github.com/MoTrPAC/MotrpacRatTraining6moData · *Nature* 2024: https://www.nature.com/articles/s41586-023-06877-w
- LINCS/CMap: https://clue.io/releases/data-dashboard · signatureSearchData vignette: https://bioconductor.org/packages/release/data/experiment/vignettes/signatureSearchData/inst/doc/signatureSearchData.html
- GSFM: https://github.com/MaayanLab/gsfm · https://gsfm.maayanlab.cloud
- MetaMEx: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6981202/ · https://github.com/NicoPillon/MetaMEx
- Mouse multi-omics CMap mimetic study: https://pubmed.ncbi.nlm.nih.gov/42033849/
