"""Candidate prioritisation rule (used by the FAIR export and the website).

1. Approved drug (Drug Repurposing Hub clinical_phase == "Launched").
2. Exercise match: overall rank in the top MATCH_TOP fraction of all screened compounds
   (RRA across informative signatures).
3. Ranked by mean benefit index (adaptive-module minus stress-module agreement) across the
   informative signatures of that exercise type.
4. Flagged rather than silently dropped: cytotoxic / generic-stress mechanisms, which match
   the stress component of acute exercise in cancer cell lines, and drugs profiled in fewer
   than MIN_CELLS cell lines. Flagged drugs are excluded from the candidate list.
"""
import re

import numpy as np
import pandas as pd

MATCH_TOP = 0.05
MIN_CELLS = 3
CYTOTOXIC_MOA = re.compile(
    r"tubulin|microtubule|proteasome|topoisomerase|alkylating|ribonucleotide reductase|thymidylate synthase|"
    r"protein synthesis inhibitor|HDAC inhibitor|DNA synthesis|antimetabolite|PKC activator|"
    r"protein phosphatase inhibitor|RNA polymerase inhibitor|dihydrofolate reductase|CDK inhibitor|HSP inhibitor|"
    r"DNA intercalat|aurora kinase|PLK inhibitor|exportin",
    re.I,
)
CARDIAC_GLYCOSIDE_TARGET = re.compile(r"\bATP1A\d", re.I)


def flag_cytotoxic(moa: pd.Series, target: pd.Series) -> pd.Series:
    m = moa.fillna("").astype(str).str.contains(CYTOTOXIC_MOA)
    t = target.fillna("").astype(str).str.contains(CARDIAC_GLYCOSIDE_TARGET)
    return m | t


def candidate_table(res: pd.DataFrame, cpd: pd.DataFrame, overall: pd.DataFrame,
                    informative: list, modality: str) -> pd.DataFrame:
    n = len(cpd)
    r = res[(res.modality == modality) & res.signature.isin(informative)]
    best = r.loc[r.groupby("drug_id")["rank"].idxmin(), ["drug_id", "signature", "rank"]] \
        .rename(columns={"signature": "best_signature", "rank": "best_signature_rank"})
    agg = r.groupby("drug_id").agg(mean_benefit=("benefit", "mean"), mean_adaptive_z=("adaptive_z", "mean"),
                                   mean_stress_z=("stress_z", "mean"),
                                   n_signatures_fdr05=("perm_fdr", lambda s: int((s < 0.05).sum())))
    t = (overall[overall.modality == modality][["drug_id", "rank", "rho", "median_rank"]]
         .rename(columns={"rank": "match_rank", "rho": "match_rho"})
         .merge(cpd[["drug_id", "name", "clinical_phase", "moa", "target", "n_cells", "n_sigs"]], on="drug_id")
         .merge(agg, left_on="drug_id", right_index=True)
         .merge(best, on="drug_id"))
    t["approved"] = t.clinical_phase == "Launched"
    t["strong_match"] = t.match_rank <= MATCH_TOP * n
    t["flag_cytotoxic"] = flag_cytotoxic(t.moa, t.target)
    t["flag_low_evidence"] = t.n_cells < MIN_CELLS
    t["candidate"] = t.approved & t.strong_match & ~t.flag_cytotoxic & ~t.flag_low_evidence
    t["candidate_rank"] = np.nan
    cand = t[t.candidate].sort_values(["mean_benefit", "match_rank"], ascending=[False, True])
    t.loc[cand.index, "candidate_rank"] = np.arange(1, len(cand) + 1)
    t.insert(0, "exercise", modality)
    return t.sort_values(["candidate_rank", "match_rank"], na_position="last").reset_index(drop=True)
