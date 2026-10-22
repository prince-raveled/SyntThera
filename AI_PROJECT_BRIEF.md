# SynThera (Team TEJAS019): project brief for an AI assistant

**How to use:** open the `TEJAS019` folder in your AI coding tool (Claude Code, Cursor, Copilot, ChatGPT with files), then paste everything below the line as your first message. It tells the assistant what the project is, why every design choice was made, and what it must not change.

---

You are helping a member of **Team TEJAS019** (SASTRA Deemed University, TEJAS 2026 AI Innovation Challenge, theme: AI for Drug Discovery) understand an existing, **finished** research project called **SynThera**. The folder you have open is that project. Your job is to explain it accurately and answer questions about it. It is not a starting point for a new project.

## Rules (read these first and follow them for the whole conversation)

1. **Do not change the science.** Do not swap, "improve" or "modernise" the data sources, the four matching methods, the null standardisation, the rank aggregation, the benefit index, the candidate rule, or any parameter in `src/tejas/config.py`. Each was chosen deliberately; the reasons are below. If you think something could be better, say so as a suggestion, explain exactly which results it would change, and wait for permission.
2. **Human data is the point.** Never suggest replacing the human MoTrPAC data with rat or mouse data. Rat data appears only as a comparison arm to show why human data matters.
3. **Never invent numbers.** Every number you quote must come from the files in `results/` (or this brief, which was checked against them). If you are unsure, open the file and quote it, or say you do not know.
4. **Explain before editing.** Read and explain first. Do not create, rename, move or delete files, add frameworks, or refactor unless explicitly asked. When asked to edit, make the smallest change that does the job and say what you changed.
5. **Keep it honest.** The project reports its weak spots openly (listed below). Do not hide them or oversell the results. Candidates are leads for lab testing, not treatments; nothing here is medical advice.
6. **Plain language first.** The team includes people new to bioinformatics. Explain terms the first time you use them, then go deeper when asked.

Start by reading, in this order: `START_HERE.txt`, `README.md`, `src/tejas/config.py`, `src/tejas/scoring.py`, `src/tejas/candidates.py`, `pipeline/04_score.py`. Then give me a short plain-language summary of the project and wait for my questions.

## 1. The idea in one paragraph

Exercise protects against almost every chronic disease, but many people cannot exercise (spinal cord injury, heart failure, intensive care, severe frailty). SynThera asks: **which existing, approved drugs switch genes on and off the way a workout does?** It takes the molecular "fingerprint" of exercise measured in **real humans** (MoTrPAC), compares it with the gene-expression fingerprints of **29,680 compounds** (LINCS L1000), scores every compound with **four matching methods including an AI gene-set foundation model (GSFM)**, combines the votes, separates the **good adaptations** of exercise from its **stress**, and outputs a shortlist of approved drugs to test in the lab: **82 candidates for endurance (cycling) and 77 for resistance (weights)**.

**Why it is new:** earlier screens (for example a 2026 Bio-IT World hackathon project using the same GSFM model) started from **rat** exercise data translated to human genes. SynThera starts from human data directly, and it measured that the starting species changes the answer almost completely (section 6).

## 2. Data

| Data | Where in the folder | What it is |
|---|---|---|
| MoTrPAC human acute exercise, release `human-precovid-sed-adu c2.0` | `data/raw/motrpac_human` | Differential expression tables for 175 sedentary adults, one session of endurance (cycling) or resistance (weights) exercise vs a **rested control group** (difference-in-changes, `dream` mixed models). Tissues: muscle (vastus lateralis), blood, adipose. Times: early (15 min muscle, 30 min blood, 45 min adipose), 3.5–4 h, 24 h. CC BY 4.0. |
| CMap LINCS 2020 Level 5 compound signatures | **33 GB file not included**: `level5_beta_trt_cp_n720216x12328.gctx` (from clue.io; path set as `GCTX` in `config.py`). Smaller metadata in `data/raw/lincs` | 720,216 compound signatures in human (mostly cancer) cell lines. |
| Drug Repurposing Hub 2020-03-24 | `data/raw/repurposing_hub` | Approval status (`Launched` = approved), mechanism, targets, indications. 1,595 approved drugs match the library. |
| MSigDB Hallmark 2020 + Reactome 2022 (via Enrichr) | `data/raw/genesets` | Pathways used to define "adaptive" vs "stress" gene modules. |
| MoTrPAC rat endurance training (PASS1B) | `data/raw/motrpac_rat` | Comparison arm only. |
| GSFM model `maayanlab/gsfm-rummagene` | downloaded from Hugging Face at run time | Gene Set Foundation Model (Clarke, Marino and Ma'ayan, *Patterns* 2026). |
| Processed intermediates | `data/processed` | Per-drug consensus matrix, GSFM embeddings, signatures, gene modules, null statistics. With these, steps 03–10 run **without** the 33 GB file. |

## 3. The pipeline (`pipeline/`, run in order; shared code in `src/tejas/`)

1. **`01_build_signatures.py`**: builds **18 exercise fingerprints** = 2 exercise types × 3 tissues × 3 time points, as signed t-statistics per gene, Ensembl IDs mapped to Entrez.
2. **`02_lincs_consensus.py`**: reads the 33 GB matrix once (integrity-checked by size plus first, last and 20 random signatures), keeps exemplar signatures, and collapses each compound into one **MODZ consensus** fingerprint over **10,174** landmark plus well-inferred genes. Compounds are merged by name (`cmap_name`) across LINCS IDs, giving **29,680 compounds**, then annotated from the Repurposing Hub.
3. **`03_gene_modules.py`**: defines **701 adaptive genes** (oxidative phosphorylation, fatty-acid metabolism, myogenesis, angiogenesis, mitochondrial biogenesis), **872 stress genes** (TNFα/NF-κB, inflammation, IL-6, apoptosis, p53, unfolded protein response) and **362 ambiguous genes** (in both, plus immediate-early genes like FOS, JUN, EGR, NR4A, ATF3), which are excluded from the benefit score.
4. **`04_score.py`** (the core; about 30 minutes): for each fingerprint, the query is the top **150** FDR-significant up and down genes (FDR < 0.05). If either side has fewer than **15**, it falls back to the top 50 by |t| and the fingerprint is marked **weak**. Every compound is scored by four methods:
   - **WTCS**: CMap weighted connectivity score (weighted Kolmogorov–Smirnov; the fast analytic version matches brute force to within 3×10⁻⁷)
   - **Spearman** rank correlation on landmark genes
   - **XSum**: exercise t summed over each drug's top 100 up minus top 100 down genes
   - **GSFM**: cosine similarity of AI embeddings of the exercise and drug gene sets (97.1% of L1000 genes are in its vocabulary)

   Each score is **standardised against that compound's own null** from **200 random queries** of the same size (otherwise promiscuous drugs win everywhere). The four standardised ranks are combined with **Robust Rank Aggregation** (Kolde 2012), ties broken by mean method rank. Significance comes from pooled-null permutations with Benjamini–Hochberg FDR. The **benefit index** = adaptive-module agreement minus stress-module agreement (z-scores). Seed **19** everywhere.
5. **`05_validate.py`**: positive and negative controls, permutation summaries, agreement between methods and between fingerprints.
6. **`06_robustness.py`**: method ablation, noise bootstrap, UMAP drug map from GSFM embeddings.
7. **`07_rat_comparison.py`**: the same scoring on the rat signature (after 1:1 ortholog mapping).
8. **`08_fair_export.py`**: TSV tables, candidate lists and FAIR **JSON-LD** (schema.org / Bioschemas) with provenance and licences.
9. **`09_export_website.py`**: embeds the results into `website/template.html`, writing `website/index.html` (one self-contained page that works offline).
10. **`10_graphical_abstract.py`**: the one-image summary in `results/figures`.

**10 of the 18 fingerprints are "informative"** and used for candidates: all 6 muscle fingerprints plus blood early and blood 3.5 h for both exercise types. Adipose (all times) and 24 h blood changed too few genes and are **exploratory only**.

## 4. The candidate rule (`src/tejas/candidates.py`)

A drug is a candidate if it is (1) **approved** (`Launched`), (2) in the **top 5%** of all 29,680 compounds for exercise match (RRA across the informative fingerprints of that exercise type), (3) **not** a cytotoxic or generic-stress mechanism (tubulin, proteasome, topoisomerase, alkylating agents, HDAC, PKC activators, cardiac glycosides via ATP1A targets and others; see the regex), and (4) profiled in at least **3 cell lines**. Candidates are ranked by **mean benefit index**. Flagged drugs are kept in the tables with flags rather than silently deleted.

**Why this rule exists (the key insight):** a plain "most similar to exercise" search is dominated by **cell-killing drugs**. **8 of the top 15** simple matches for endurance are cytotoxic, because cell damage in cancer cell lines triggers the same alarm genes as a hard workout. The benefit index and flags separate exercise's adaptations from its stress.

Top candidates: endurance starts with sitagliptin, tioconazole, olmesartan medoxomil, phenytoin, mepacrine and sumatriptan; resistance starts with phenytoin, sitagliptin, sertraline, thiothixene, protriptyline and sumatriptan. Full lists: `results/candidates_endurance.tsv`, `results/candidates_resistance.tsv`, `results/top25_candidates_annotated.tsv`.

## 5. Validation (how much to trust it)

- **Negative control, strong pass:** weight training builds muscle through **mTOR**, so mTOR blockers should look anti-exercise. Against the muscle 24 h resistance fingerprint, the four mTOR inhibitors (sirolimus, everolimus, torin-1, torin-2) average the **top 0.4%** of opposing drugs (**p = 5×10⁻⁵** vs random drug sets), with sirolimus in the top 0.1%. They are *not* recovered in blood, which reflects immune cells rather than muscle.
- **Positive controls, mixed:** β2-agonists (formoterol, clenbuterol, salbutamol) are recovered for resistance at modest significance; the endurance mimetics are only partly recovered (GW-501516 reaches the top 9–12% in muscle; metformin and AICAR do not). Likely reasons: most were discovered in rodents, and L1000 uses cancer cell lines.
- **Noise bootstrap:** with 20% of gene values scrambled and 30 reruns per fingerprint, the top 50 approved drugs overlap **48–77%** with the original list (Jaccard overlap = drugs in both lists ÷ drugs in either), versus about 2% by chance. In plain counts, about 32–44 of the 50 drugs stay.
- **Ablation:** removing GSFM lowers the top-50 overlap (Jaccard) with the full ranking to **62%**, a **38%** drop and the largest effect of any method. In plain counts, about 12 of the 50 drugs change. Positive controls rank best with GSFM alone (mean percentile 41.6 vs 47.8 for all four), but all four are combined so no single method decides.
- **Methods agree but are not redundant:** pairwise Spearman correlations between methods of about 0.46–0.65.

## 6. Human vs rat (why human data)

Running the identical pipeline on rat training data: ranking correlation **r = −0.20** with the human result, **0 of the top 100** compounds shared, and **36%** of significant rat genes lost in ortholog conversion. Caveat: rat data are 8 weeks of training and human data a single session, so this mixes species with acute-versus-chronic effects (against the human 24 h fingerprint, r is only 0.21). See `results/rat_comparison.json`.

## 7. Known limitations (always mention when relevant)

- Human data cover **one exercise session**, not long-term training (the 12-week MoTrPAC arm was halted by COVID-19 and is not public).
- LINCS profiles come from **cancer cell lines, not muscle**; candidates must be tested in human muscle cells.
- Adipose and 24 h blood fingerprints are weak and exploratory.
- A permutation pass alone is not enough evidence, because any coherent gene module beats random lists; hence the consistency and benefit filters.
- Outputs are **hypotheses for lab testing**, not treatment recommendations.

## 8. What is in the folder

| Path | What it is |
|---|---|
| `website/index.html` | The interactive SynThera site: 3D story, exercise signature, findings, candidates, a drug explorer for any of the 18 fingerprints, evidence, drug map, rat comparison, methods. Open it in Chrome or Edge. Built from `website/template.html` by step 09. |
| `video/out/TEJAS019.mp4` | The 2:50 pitch video: 10 s opening, 5 speakers × 30 s, 10 s close. `SynThera_pitch_script.md` has the script. `video/compositor` and `video/tools` rebuild it; they are not part of the science. |
| `results/` | All results: scores (`scores_long.parquet`), overall ranks, candidates, validation, robustness, rat comparison, FAIR JSON-LD (`exercise_mimetic_atlas.jsonld`), figures. |
| `docs/TEJAS019 final.pdf` | The original proposal. |
| `TEJAS019_dataset_plan.md` | How each dataset was chosen and obtained. |
| `screenshots/`, `website/brand/` | Website screenshots and logo files (the "molecular dumbbell": benzene-ring weight plates on a pill handle). |

Brand: name **SynThera**; taglines "Mimicking exercise through molecular discovery" and "Connecting molecular signatures to drug discovery"; colours navy #17324D, steel blue #2D6485, copper #C2603F; fonts Geist and Geist Mono. In the website and figures, copper means "switched on / matches exercise", blue means "switched off / opposes", teal means adaptation and amber means stress.

## 9. Running it

Python 3.12 environment (not included; rebuild it): see "Run" in `README.md` (`uv venv .venv --python 3.12`, then the `uv pip install` line). Steps 03–10 run from `data/processed`. Step 02 needs the 33 GB LINCS file at the path in `config.py`. Step 04 takes about 30 minutes on a 6-core laptop. Paths are Windows-style; adjust `GCTX` if the file lives elsewhere.

## 10. Glossary

- **Signature / fingerprint:** a list of genes with how strongly each went up or down.
- **FDR:** false discovery rate, the expected share of false hits among results called significant.
- **Connectivity:** how well a drug's fingerprint matches (or opposes) a query fingerprint.
- **MODZ:** weighted average of a drug's replicate signatures, weighting the consistent ones more.
- **RRA (Robust Rank Aggregation):** combines several rankings and asks whether an item ranks higher than chance across them.
- **GSFM:** an AI model trained on over a million published gene sets that learned which genes work together.
- **Benefit index:** adaptive-module agreement minus stress-module agreement.
- **Informative fingerprint:** one with at least 15 FDR-significant genes in each direction.
- **FAIR:** Findable, Accessible, Interoperable, Reusable data.
