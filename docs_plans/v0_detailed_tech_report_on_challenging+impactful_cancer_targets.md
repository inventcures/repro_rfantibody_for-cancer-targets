# Tech Report Plan: De Novo VHH Antibody Design for Challenging Cancer Targets Using RFAntibody

> **Output**: `docs_plans/v0_detailed_tech_report_on_challenging+impactful_cancer_targets.tex` + PDF
> **Status**: WAITING — all 10 campaigns must complete on Modal before writing
> **Data source**: Modal persistent Volume `rfab-cancer-drivers-results`

---

## Report Structure (LaTeX, scientific_report.sty)

### 1. Title & Abstract
- De novo VHH nanobody design targeting 10 cancer-associated antigens using RFAntibody
- Targets span GBM, PDAC, neuroblastoma, MPNST, pediatric glioma
- 100K total backbones → 500K sequences → filtered candidates

### 2. Introduction
- Unmet need in these 5 cancer indications (from v1_specs.md)
- RFAntibody pipeline (Watson et al., Nature 2025): RFdiffusion → ProteinMPNN → RF2
- Why VHH format (tissue penetration, BBB crossing, bispecific building blocks)
- Why these 10 targets (structural readiness + clinical impact ranking)

### 3. Methods
- Target selection rationale (structural readiness matrix)
- PDB template selection + epitope/hotspot identification via contact analysis
- Campaign configuration (CDR loop ranges, truncation, filtering thresholds)
- Computational pipeline: RFAntibody 3-stage on Modal A100-80GB × 10 parallel
- Filtering criteria: pAE < threshold, RMSD < threshold, ddG < threshold
- Ranking: composite score

### 4. Results (per campaign, populated from Volume data)
For each of the 10 campaigns:
- Stage timing breakdown (RFdiffusion, ProteinMPNN, RF2)
- Design funnel: total → filtered → ranked
- Score distributions (pAE, RMSD, ddG) — histograms + summary stats
- Top candidates table (rank, tag, pAE, RMSD, ddG, composite)
- Pass rate comparison across campaigns

### 5. Cross-Campaign Analysis
- Comparative pass rates across 10 targets
- Which targets produced highest-confidence binders
- Correlation between target properties (size, epitope area, hydrophobicity) and success rate
- Compute cost breakdown per campaign

### 6. Discussion
- Structural factors affecting design success
- Limitations: single-epitope per target, no experimental validation yet
- Comparison to Watson et al. benchmarks
- Potential for bispecific combinations (e.g., MSLN-nterm + MSLN-cterm)

### 7. Conclusion & Next Steps
- Experimental validation priorities (yeast display, SPR)
- Most promising candidates per indication

---

## Data Required (from Volume after completion)

```
For each campaign in {b7h3_vhh, cd47_vhh, ceacam5_vhh, egfr_cetuximab_vhh,
                       egfrviii_vhh, epha2_vhh, gpc2_vhh, her2_domIV_vhh,
                       msln_cterm_vhh, msln_nterm_vhh}:

  cancer_drivers/{campaign}_v1/logs/{campaign}_result.json
    → stage_times, elapsed_seconds, design_count, success

  cancer_drivers/{campaign}_v1/analysis/ranked_candidates.csv
    → rank, tag, pae, rmsd, ddg, composite_score

  cancer_drivers/{campaign}_v1/analysis/report.html
    → filtering summary, score distributions

  .batch_state.json
    → aggregate timing and status for all campaigns
```

## Visualization Guidelines

**ALL figures MUST follow Saloni Dattani's data visualization guidelines** (saved in
`docs_plans/saloni_dataviz_guidelines.md`). Key rules enforced:

- Horizontal text everywhere — no rotated axis labels
- Direct labeling instead of legends (labels adjacent to data elements)
- Colorblind-safe palette (Paul Tol's bright scheme)
- Every figure standalone: title, subtitle, units, data source ON the figure
- No 3D, no pie charts, no dual-axis with mismatched scales
- Small multiples (faceting) for multi-campaign comparisons
- "Clearer, not simpler" — complexity is fine if guided
- Check every figure against the 6 questions:
  1. Is the chart type meaningful? 2. Can it be clearer? 3. Is complexity guided?
  4. Does it work standalone? 5. Is the presentation justifiable? 6. Is it reproducible?

### Planned Figures

| # | Figure | Type | Saloni Principle |
|---|--------|------|-----------------|
| 1 | Graphical abstract | Schematic (scientific-schematics) | Standalone summary |
| 2 | Pipeline overview | Flow diagram (scientific-schematics) | Guide through complexity |
| 3 | Target structural readiness matrix | Heatmap with direct labels | Meaningful chart type |
| 4 | Score distributions (pAE, RMSD, ddG) | Histograms, small multiples per campaign | Faceting for clarity |
| 5 | Design funnel per campaign | Horizontal waterfall/bar | Meaningful for attrition data |
| 6 | Cross-campaign pass rates | Horizontal bar, sorted by rate | Logical ordering |
| 7 | Per-stage timing breakdown | Stacked horizontal bar | Composition + total |
| 8 | Composite score vs individual metrics | Scatter with direct annotations | Show correlations |
| 9 | Target size vs success rate | Annotated scatter | Guide interpretation |
| 10 | Top candidates summary | Clean formatted table | Tables for exact values |

## Generation Plan

1. Download all result data from Volume: `modal volume get rfab-cancer-drivers-results ...`
2. Parse result JSONs and CSVs into summary DataFrames
3. Generate figures following `docs_plans/saloni_dataviz_guidelines.md` strictly
4. Write LaTeX using scientific_report.sty with scientific-writing skill
5. Compile with xelatex → PDF
