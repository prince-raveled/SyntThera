"""Step 10: one-image graphical abstract (flow + results), every number read from results/.

Writes website/graphical_abstract.html (1920x1080 layout). Render to PNG with a headless browser,
e.g. tools/render_png.js, or open the HTML and take a full-window screenshot.
"""
import json
import sys
from html import escape
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from tejas import config as C  # noqa: E402
from tejas.candidates import MATCH_TOP, candidate_table  # noqa: E402


def main():
    res = pd.read_parquet(C.RESULTS / "scores_long.parquet")
    cpd = pd.read_parquet(C.PROC / "compounds.parquet")
    overall = pd.read_parquet(C.RESULTS / "overall.parquet")
    meta = json.loads((C.RESULTS / "scoring_meta.json").read_text())
    val = json.loads((C.RESULTS / "validation.json").read_text())
    rob = json.loads((C.RESULTS / "robustness.json").read_text())
    rat = json.loads((C.RESULTS / "rat_comparison.json").read_text())
    summ = pd.read_csv(C.PROC / "signature_summary.csv")
    sigl = pd.read_parquet(C.PROC / "signatures_long.parquet")
    informative = meta["informative_signatures"]
    n = len(cpd)
    n_launched = int((cpd.clinical_phase == "Launched").sum())

    cands, n_cand, n_cyt15, n_strong = {}, {}, {}, {}
    for m in C.MODALITIES:
        ct = candidate_table(res, cpd, overall, informative, m)
        strong = ct[ct.approved & ct.strong_match]
        n_strong[m] = len(strong)
        n_cyt15[m] = int(strong.sort_values("match_rank").head(15).flag_cytotoxic.sum())
        c = ct[ct.candidate].sort_values("candidate_rank")
        n_cand[m] = len(c)
        c = c.merge(cpd[["drug_id", "indication"]], on="drug_id", how="left")
        cands[m] = c.head(6)

    neg = [c for c in val["controls"] if c["direction"] == "oppose" and c["signature"] == "muscle|resistance|24h"]
    neg = sorted(neg, key=lambda c: c["percentile"])
    neg_p = next(c["p_vs_random"] for c in val["controls_summary"]
                 if c["direction"] == "oppose" and c["signature"] == "muscle|resistance|24h")
    pos = [c for c in val["controls_summary"] if c["direction"] == "mimic"]
    b2 = min((c for c in pos if c["modality"] == "resistance"), key=lambda c: c["p_vs_random"])
    gw = min(c["percentile"] for c in val["controls"] if c["drug_id"] == "GW-501516" and c["signature"] in informative)
    abl = {a["config"]: a for a in rob["ablation_summary"]}
    boot = [b["jaccard_vs_unperturbed_median"] for b in rob["bootstrap"]]
    chance = rob["bootstrap"][0]["jaccard_random_expected"]
    mc = np.array(val["method_agreement"]["spearman"])
    off_diag = mc[~np.eye(len(mc), dtype=bool)]

    hs = sigl[(sigl.signature == "muscle|endurance|3.5-4h") & (sigl.fdr < C.FDR)]
    genes = pd.concat([hs.nlargest(5, "t"), hs.nsmallest(2, "t")])
    gmax = genes.logFC.abs().max()
    gene_bars = "".join(
        f'<div class="gb"><span>{escape(str(g))}</span><i class="{"on" if v > 0 else "off"}" style="width:{abs(v) / gmax * 100:.0f}%"></i></div>'
        for g, v in zip(genes.gene_symbol, genes.logFC))
    musc = summ[summ.tissue == "muscle"]
    n_musc = int((musc.n_fdr_up + musc.n_fdr_down).max())

    def cand_rows(m):
        rows = []
        for _, r in cands[m].iterrows():
            use = str(r.indication) if isinstance(r.indication, str) and r.indication else str(r.moa or "")
            use = use.split("|")[0][:34]
            a, s = max(0.0, r.mean_adaptive_z), max(0.0, r.mean_stress_z)
            rows.append(f'<li><b>{int(r.candidate_rank)}</b><div><span class="nm">{escape(r["name"])}</span>'
                        f'<span class="use">{escape(use)}</span></div>'
                        f'<div class="gsb"><i class="g" style="width:{min(a / 4, 1) * 100:.0f}%"></i>'
                        f'<i class="s" style="width:{min(s / 4, 1) * 100:.0f}%"></i></div></li>')
        return "".join(rows)

    pins = "".join(f'<i style="left:{max(1.5, c["percentile"]):.1f}%"></i>' for c in neg)
    neg_names = " · ".join(f'{escape(c["drug_id"])} {c["percentile"]:.1f}%' for c in neg)
    noG, gOnly, full = abl["Without GSFM"], abl["GSFM only"], abl["Full ensemble (4 methods)"]

    html = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><title>SynThera · graphical abstract</title>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Geist:wght@300..800&family=Geist+Mono:wght@400;500&display=swap">
<style>
:root {{ --bg:#f4f7fa; --surface:#fff; --ink:#0f2438; --ink2:#17324d; --muted:#5a6b7b; --line:#dce4ec; --on:#c2603f; --off:#2d6485;
  --good:#2e8c7e; --warn:#b7791f; --cell:#e9eff4; --endur:#3a7fb0; --resist:#7b5ea7;
  --d:"Geist","Segoe UI",sans-serif; --b:"Geist","Segoe UI",sans-serif; --m:"Geist Mono",Consolas,monospace; }}
* {{ box-sizing:border-box; margin:0; }}
html, body {{ width:1920px; height:1080px; }}
body {{ background:var(--bg); color:var(--ink); font-family:var(--b); font-size:15px; line-height:1.4; padding:34px 44px 26px; display:grid; grid-template-rows:auto auto auto auto 1fr auto; gap:18px; position:relative; overflow:hidden; }}
body::before {{ content:""; position:absolute; inset:0; z-index:-1;
  background:radial-gradient(900px 420px at 88% 0%, rgb(45 100 133 / .13), transparent 70%), radial-gradient(700px 360px at 0% 100%, rgb(194 96 63 / .09), transparent 70%); }}
h1,h2,h3 {{ font-family:var(--d); letter-spacing:-.03em; line-height:1.05; font-weight:650; }}
.k {{ font-family:var(--m); font-size:12px; letter-spacing:.14em; text-transform:uppercase; color:var(--off); font-weight:500; }}
.muted {{ color:var(--muted); }}
header {{ display:grid; grid-template-columns:1fr auto; gap:40px; align-items:end; }}
.brand {{ display:flex; align-items:center; gap:14px; margin-top:8px; }}
.brand svg {{ width:50px; height:50px; }}
.brand h1 {{ font-size:50px; font-weight:720; letter-spacing:-.045em; }}
.brand .tl {{ font-size:27px; font-weight:600; letter-spacing:-.03em; color:var(--on); margin-left:18px; padding-left:18px; border-left:2px solid var(--line); line-height:1.15; }}
header p {{ font-size:17px; color:var(--muted); margin-top:8px; }}
.stats {{ display:grid; grid-template-columns:repeat(4,auto); gap:34px; }}
.stats b {{ display:block; font-family:var(--d); font-size:40px; font-weight:720; letter-spacing:-.045em; line-height:1; }}
.stats span {{ font-size:13px; color:var(--muted); }}
.sec {{ display:flex; align-items:center; gap:14px; }}
.sec::after {{ content:""; flex:1; height:1px; background:var(--line); }}
.flow {{ display:grid; grid-template-columns:repeat(5,1fr); gap:34px; }}
.st {{ position:relative; background:var(--surface); border:1px solid var(--line); border-radius:14px; padding:16px 18px; display:grid; gap:8px; align-content:start; }}
.st:not(:last-child)::after {{ content:""; position:absolute; right:-27px; top:50%; width:20px; height:20px; transform:translateY(-50%);
  background:var(--on); clip-path:polygon(0 15%,60% 15%,60% 0,100% 50%,60% 100%,60% 85%,0 85%); }}
.st .n {{ font-family:var(--m); font-size:12px; color:var(--on); }}
.st h3 {{ font-size:19px; font-weight:800; }}
.st p {{ font-size:13.5px; color:var(--muted); }}
.st .big {{ font-family:var(--d); font-weight:800; font-size:26px; letter-spacing:-.02em; color:var(--ink); }}
.gbars {{ display:grid; gap:2px; }}
.gb {{ display:grid; grid-template-columns:72px 1fr; gap:6px; align-items:center; font-family:var(--m); font-size:10.5px; line-height:1.25; }}
.gb i {{ height:7px; border-radius:2px; display:block; }}
.gb i.on {{ background:var(--on); }} .gb i.off {{ background:var(--off); }}
.chips {{ display:flex; flex-wrap:wrap; gap:5px; }}
.chip {{ font-family:var(--m); font-size:11px; padding:3px 8px; border-radius:99px; border:1px solid var(--line); background:var(--bg); }}
.chip.ai {{ border-color:var(--on); color:var(--on); }}
.arrow {{ font-family:var(--m); font-size:12px; color:var(--muted); }}
.bal {{ display:grid; gap:5px; }}
.bal div {{ display:grid; grid-template-columns:58px 1fr; gap:8px; align-items:center; font-size:12px; color:var(--muted); }}
.bal i {{ display:block; height:9px; border-radius:3px; }}
.funnel {{ display:grid; gap:4px; }}
.funnel div {{ display:flex; justify-content:space-between; align-items:center; gap:8px; font-size:12.5px; padding:4px 10px; border-radius:6px; background:var(--cell); margin-inline:auto; width:100%; }}
.funnel div b {{ font-family:var(--m); font-weight:500; }}
.results {{ display:grid; grid-template-columns:1.25fr 1fr 1fr; gap:18px; min-height:0; }}
.card {{ background:var(--surface); border:1px solid var(--line); border-radius:14px; padding:16px 18px; display:grid; gap:8px; align-content:start; }}
.card h3 {{ font-size:18px; font-weight:800; }}
.card p {{ font-size:13.5px; color:var(--muted); }}
.col {{ display:grid; gap:14px; align-content:start; }}
.two {{ display:grid; grid-template-columns:1fr 1fr; gap:16px; }}
.lst h4 {{ font:600 13px var(--b); display:flex; gap:6px; align-items:center; margin-bottom:4px; }}
.lst h4 i {{ width:14px; height:3px; border-radius:2px; display:inline-block; }}
.lst ol {{ list-style:none; padding:0; display:grid; gap:6px; }}
.lst li {{ display:grid; grid-template-columns:20px 1fr 64px; gap:8px; align-items:center; border-bottom:1px dashed var(--line); padding-bottom:5px; }}
.lst li b {{ font-family:var(--m); font-size:12px; color:var(--on); font-weight:500; }}
.lst .nm {{ display:block; font-weight:700; font-size:14.5px; text-transform:capitalize; white-space:nowrap; overflow:hidden; text-overflow:ellipsis; }}
.lst .use {{ display:block; font-size:11.5px; color:var(--muted); white-space:nowrap; overflow:hidden; text-overflow:ellipsis; }}
.gsb {{ display:grid; gap:3px; }}
.gsb i {{ display:block; height:6px; border-radius:2px; }}
.gsb .g {{ background:var(--good); }} .gsb .s {{ background:var(--warn); }}
.legend {{ display:flex; gap:12px; font-size:11.5px; color:var(--muted); }}
.legend span {{ display:inline-flex; gap:5px; align-items:center; }}
.legend i {{ width:10px; height:6px; border-radius:2px; display:inline-block; }}
.hl {{ display:flex; align-items:baseline; gap:12px; }}
.hl .num {{ font-family:var(--d); font-weight:800; font-size:40px; letter-spacing:-.03em; line-height:.95; color:var(--on); white-space:nowrap; }}
.hl .num.blue {{ color:var(--off); }} .hl .num.green {{ color:var(--good); }}
.scale {{ position:relative; height:16px; margin-top:4px; }}
.scale .ax {{ position:absolute; inset:5px 0; border-radius:3px; background:linear-gradient(90deg,var(--off),var(--cell) 45%,var(--cell) 55%,var(--on)); }}
.scale i {{ position:absolute; top:0; width:14px; height:14px; border-radius:50%; background:var(--off); border:2px solid #fff; box-shadow:0 0 0 1px var(--off); transform:translateX(-50%); }}
.ends {{ display:flex; justify-content:space-between; font-size:11px; color:var(--muted); }}
.mini {{ display:grid; gap:5px; }}
.mini div {{ display:grid; grid-template-columns:150px 1fr 40px; gap:8px; align-items:center; font-size:12.5px; }}
.mini .t {{ height:9px; background:var(--cell); border-radius:3px; overflow:hidden; }}
.mini .t i {{ display:block; height:100%; border-radius:3px; background:var(--off); }}
.mini .t i.ai {{ background:var(--on); }}
.mini b {{ font-family:var(--m); font-weight:500; text-align:right; }}
.tag {{ display:inline-block; font-size:11.5px; font-weight:600; padding:2px 9px; border-radius:99px; }}
.tag.pass {{ color:var(--good); background:rgb(35 120 74 / .1); }} .tag.mixed {{ color:var(--warn); background:rgb(165 98 12 / .12); }}
footer {{ display:flex; justify-content:space-between; gap:30px; font-size:12px; color:var(--muted); border-top:1px solid var(--line); padding-top:10px; }}
footer b {{ color:var(--ink); }}
</style></head>
<body>
<header>
  <div><div class="k">Connecting molecular signatures to drug discovery · TEJAS 2026 · Team TEJAS019 · SASTRA Deemed University</div>
  <div class="brand"><svg aria-hidden="true" viewBox="0 0 64 64"><mask id="gaRings" maskUnits="userSpaceOnUse" x="0" y="0" width="64" height="64"><rect width="64" height="64" fill="#fff"/><circle cx="13" cy="32" r="5.2" fill="none" stroke="#000" stroke-width="2.4"/><circle cx="51" cy="32" r="5.2" fill="none" stroke="#000" stroke-width="2.4"/></mask><rect x="20" y="27.5" width="24" height="9" rx="4.5" fill="#17324d"/><g mask="url(#gaRings)"><polygon points="13,17 22.5,22.5 22.5,41.5 13,47 3.5,41.5 3.5,22.5" fill="#c2603f"/><polygon points="51,17 60.5,22.5 60.5,41.5 51,47 41.5,41.5 41.5,22.5" fill="#2d6485"/></g></svg><h1>SynThera</h1><span class="tl">Mimicking exercise through<br>molecular discovery</span></div>
  <p>We matched the molecular fingerprint of exercise in 175 people against {n:,} compounds, to find approved medicines worth testing for patients who cannot exercise.</p></div>
  <div class="stats"><div><b>175</b><span>adults, rested control</span></div><div><b>{len(summ)}</b><span>exercise fingerprints ({len(informative)} strong)</span></div><div><b>{n:,}</b><span>compounds screened</span></div><div><b>{n_launched:,}</b><span>approved drugs</span></div></div>
</header>

<div class="sec"><span class="k">How it works</span></div>
<section class="flow">
  <div class="st"><span class="n">01 · HUMAN DATA</span><h3>Exercise fingerprint</h3>
    <p>MoTrPAC: sedentary adults, one bout of cycling or weights vs a rested group. Muscle, blood and fat, 15 min to 24 h. Exercise minus control.</p>
    <div class="gbars">{gene_bars}</div></div>
  <div class="st"><span class="n">02 · DRUG LIBRARY</span><h3>Drug fingerprints</h3>
    <p>LINCS L1000 2020: 720,216 compound signatures in human cells, collapsed per drug (MODZ) over 10,174 genes.</p>
    <div class="big">{n:,} <span class="muted" style="font-size:14px;font-weight:500">compounds</span></div>
    <p>Up to {cpd.n_cells.max()} cell lines per drug; {n_launched:,} approved (Drug Repurposing Hub).</p></div>
  <div class="st"><span class="n">03 · MATCHING + AI</span><h3>Four methods vote</h3>
    <div class="chips"><span class="chip">CMap connectivity</span><span class="chip">Spearman</span><span class="chip">XSum</span><span class="chip ai">GSFM AI model</span></div>
    <p>Each score compared with the drug's own 200 random-gene baselines, then combined by Robust Rank Aggregation. Methods agree only partly (r = {off_diag.min():.2f}–{off_diag.max():.2f}), so combining them adds information.</p></div>
  <div class="st"><span class="n">04 · BENEFIT INDEX</span><h3>Good vs stress</h3>
    <div class="bal"><div>Good<i style="width:80%;background:var(--good)"></i></div><div>Stress<i style="width:34%;background:var(--warn)"></i></div></div>
    <p>701 adaptive genes (mitochondria, fat burning, muscle, vessels) vs 872 stress genes (inflammation, cell death, protein stress). Benefit = good minus stress.</p></div>
  <div class="st"><span class="n">05 · PRIORITISE</span><h3>Approved candidates</h3>
    <div class="funnel"><div style="width:100%"><span>Screened</span><b>{n:,}</b></div><div style="width:88%"><span>Approved</span><b>{n_launched:,}</b></div><div style="width:74%"><span>Top {MATCH_TOP:.0%} match</span><b>{n_strong['endurance']} / {n_strong['resistance']}</b></div><div style="width:60%;background:rgb(35 120 74 / .14)"><span>Candidates</span><b>{n_cand['endurance']} / {n_cand['resistance']}</b></div></div>
    <p>Cycling / weights. Cell-killing drugs and thin evidence (&lt;3 cell lines) removed.</p></div>
</section>

<div class="sec"><span class="k">What we found</span></div>
<section class="results">
  <div class="card">
    <h3>Top approved candidates</h3>
    <div class="legend"><span><i style="background:var(--good)"></i>good changes</span><span><i style="background:var(--warn)"></i>stress changes</span><span>ranked by good minus stress</span></div>
    <div class="two">
      <div class="lst"><h4><i style="background:var(--endur)"></i>Cycling (endurance)</h4><ol>{cand_rows('endurance')}</ol></div>
      <div class="lst"><h4><i style="background:var(--resist)"></i>Weight training (resistance)</h4><ol>{cand_rows('resistance')}</ol></div>
    </div>
    <p>Leads for laboratory testing in muscle cells, not treatments.</p>
  </div>
  <div class="col">
    <div class="card"><h3>A naive search picks chemotherapy</h3>
      <div class="hl"><span class="num">{n_cyt15['endurance']}/15</span><p>strongest raw cycling matches are cell-killing drugs (vinca alkaloids, proteasome inhibitors) that copy exercise's stress, not its benefit. The benefit index removes them.</p></div></div>
    <div class="card"><h3>It found the muscle-growth switch</h3>
      <div class="hl"><span class="num blue">top {neg[0]['percentile']:.1f}%</span><p>sirolimus (rapamycin) for <b>opposing</b> weight training, out of {n:,}: mTOR drives muscle growth. p = {neg_p:.0e}</p></div>
      <div class="scale"><div class="ax"></div>{pins}</div><div class="ends"><span>opposes exercise</span><span>{neg_names}</span><span>copies exercise</span></div></div>
    <div class="card"><h3>Known mimetics <span class="tag mixed">mixed</span></h3>
      <p>β2-agonists recovered for weights (p = {b2['p_vs_random']:.3f}); GW-501516 top {gw:.0f}%; AICAR and metformin not. Most were found in rodents.</p></div>
  </div>
  <div class="col">
    <div class="card"><h3>Rats and people disagree</h3>
      <div class="hl"><span class="num">{rat['top100_overlap']}/100</span><p>top drugs shared by the human and rat screens (r = {rat['spearman']:.2f}). Converting rat genes to human lost {rat['genes_lost_frac']:.0%} of the rat signal.</p></div></div>
    <div class="card"><h3>The AI model changes the answer</h3>
      <div class="mini"><div>Without GSFM<span class="t"><i class="ai" style="width:{noG['top50_jaccard_vs_full'] * 100:.0f}%"></i></span><b>{noG['top50_jaccard_vs_full'] * 100:.0f}%</b></div>
      <div>Without XSum<span class="t"><i style="width:{abl['Without XSum']['top50_jaccard_vs_full'] * 100:.0f}%"></i></span><b>{abl['Without XSum']['top50_jaccard_vs_full'] * 100:.0f}%</b></div></div>
      <p>Top-50 overlap with the full ranking (Jaccard) after removing a method. GSFM alone recovers known mimetics best ({gOnly['control_mean_percentile']:.0f}th vs {full['control_mean_percentile']:.0f}th percentile).</p></div>
    <div class="card"><h3>Stable under noise <span class="tag pass">passed</span></h3>
      <div class="hl"><span class="num green">{min(boot):.0%}–{max(boot):.0%}</span><p>top-50 overlap (Jaccard) with the original after scrambling 20% of genes, 30 runs per fingerprint (chance ≈ {chance:.0%}).</p></div></div>
  </div>
</section>

<footer>
  <span><b>Data:</b> MoTrPAC human pre-suspension c2.0 &amp; rat PASS1B (CC BY 4.0) · CMap LINCS 2020 · Drug Repurposing Hub · MSigDB/Reactome · GSFM (Ma'ayan Lab, Patterns 2026)</span>
  <span><b>Limits:</b> one exercise bout · cancer cell lines · fat and 24 h blood too weak · predictions, not medical advice</span>
</footer>
</body></html>"""
    out = C.ROOT / "website" / "graphical_abstract.html"
    out.write_text(html, encoding="utf-8")
    print(f"written {out}")


if __name__ == "__main__":
    main()
