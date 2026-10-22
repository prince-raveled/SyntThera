"""Paths and analysis parameters shared by every pipeline step."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "data" / "raw"
PROC = ROOT / "data" / "processed"
RESULTS = ROOT / "results"
WEB_DATA = ROOT / "website" / "data"

MOTRPAC = RAW / "motrpac_human"
LINCS = RAW / "lincs"
HUB = RAW / "repurposing_hub"
GENESETS = RAW / "genesets"

# 33 GB LINCS 2020 Level-5 compound matrix (kept on E: because of its size)
GCTX = Path(r"E:\lincs\level5_beta_trt_cp_n720216x12328.gctx")
GCTX_BYTES = 35_518_405_386

MOTRPAC_RELEASE = "human-precovid-sed-adu c2.0"
LINCS_RELEASE = "CMap LINCS 2020 (beta), level5_beta_trt_cp"
HUB_RELEASE = "Drug Repurposing Hub 2020-03-24"

TISSUES = {"muscle": "t06-muscle", "adipose": "t11-adipose", "blood": "t04-blood-rna"}
MODALITIES = {"endurance": "ADUEndur", "resistance": "ADUResist"}
# MoTrPAC timepoint labels; "early" is 15 min (muscle), 30 min (blood) or 45 min (adipose)
TIMEPOINTS = {"early": "post_15_30_45_min", "3.5-4h": "post_3.5_4_hr", "24h": "post_24_hr"}
EARLY_LABEL = {"muscle": "15 min", "adipose": "45 min", "blood": "30 min"}

PRIMARY_TISSUE = "muscle"
PRIMARY_TIMEPOINT = "3.5-4h"

# Query gene sets
FDR = 0.05
N_QUERY = 150          # max genes per direction (CMap query convention: 10-150)
MIN_QUERY = 15         # below this many FDR-significant genes, fall back to top-|t| genes
N_DRUG_SET = 100       # top up/down genes per drug for XSum and GSFM

GSFM_MODEL = "maayanlab/gsfm-rummagene"
N_PERM = 200           # random-query permutations per signature (pooled null)
SEED = 19

# Validation controls (names as in LINCS compoundinfo `cmap_name`)
POSITIVE_CONTROLS = {
    "endurance": ["acadesine", "metformin", "GW-501516", "resveratrol"],
    "resistance": ["formoterol", "clenbuterol", "salbutamol"],
}
NEGATIVE_CONTROLS = {"resistance": ["sirolimus", "everolimus", "torin-1", "torin-2"]}


def signature_id(tissue: str, modality: str, timepoint: str) -> str:
    return f"{tissue}|{modality}|{timepoint}"


def timepoint_label(tissue: str, timepoint: str) -> str:
    return EARLY_LABEL[tissue] if timepoint == "early" else timepoint


for _d in (PROC, RESULTS, WEB_DATA, GENESETS):
    _d.mkdir(parents=True, exist_ok=True)
