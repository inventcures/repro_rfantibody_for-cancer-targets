# Deep Research: Confidence Score Failure, Data Infrastructure, and the rfab-harness Pipeline

> **Date**: 2026-03-11 | **Author**: Ashish Makani
> **Corpus**: 7 works (2 Greiff-group papers + 5 project websites)
> **Scope**: Thematic synthesis of structural plausibility limitations, antibody-antigen data infrastructure, and implications for computational antibody design pipelines

---

## 1. Executive Summary

Two papers from the Greiff group (University of Oslo / NaturalAntibody) published in early 2026 fundamentally challenge the filtering assumptions underlying our rfab-harness pipeline and most computational antibody discovery workflows:

1. **Smorodyina et al. (2026)** demonstrate that confidence scores (ipTM, pTM) from AlphaFold3, Boltz-2, and Chai-1 **fail to discriminate cognate VHH-antigen pairings from shuffled (non-binding) controls** (PR-AUC: AF3 AP=0.187, Chai-1 AP=0.067, Boltz-2 AP=0.026; random baseline ~0.011). The ipTM metric exhibits a "locked-in" effect where structural quality improvements from sampling do not translate to confidence score changes (delta-DockQ vs delta-ipTM correlations r=-0.03 to -0.04 across all tools).

2. **Czerwinski et al. (2025)** introduce the Antigen-Specific Antibody Database (ASD), aggregating 1,097,946 antibody-antigen interactions from 15 sources across 9,575 unique antigens. The database includes a VHH/nanobody-specific subset (1,258 structures, 390 antigens) directly relevant to our campaign designs, but suffers from HER2 dominance (524K/1.1M records from a single "buzz" dataset) and heterogeneous affinity measurement (48% "fuzzy" categorical).

These findings converge with our prior work across 5 project analyses to expose a systemic gap: **the field's standard generate-then-filter-by-confidence workflow selects for structural self-consistency, not binding specificity**. Our rfab-harness pipeline (pAE < 10, RMSD < 2.0 A, composite ranking by 0.4*pAE + 0.3*RMSD + 0.3*ddG) is directly implicated, as all three metrics measure structural plausibility rather than interaction specificity.

---

## 2. Corpus Table

| ID | Work | Authors | Year | Central Claim | Evidence Type |
|----|------|---------|------|---------------|---------------|
| P1 | Structural Plausibility Without Binding Specificity | Smorodyina, Ali, Brumat, Salicari, Miklavc, Kappassov, Fu, Sormanni, de Marco, Greiff | 2026 | Confidence scores from AF3/Boltz-2/Chai-1 fail to discriminate real from shuffled VHH-antigen pairings | Experimental benchmark (106 complexes, 1.8M predictions) |
| P2 | ASD: Antigen-Specific Antibody Database | Czerwinski, Dudzic, Wojtowicz, Jaszczyszyn, Bielska, Wrobel, Demharter, Spreafico, Greiff, Krawczyk | 2025 | Unified repository of 1.1M antibody-antigen interactions enables model development | Database (15 sources, 25 datasets) |
| P3 | Dogfooding rfab-harness on Cancer Targets (preprint) | Makani | 2026 | 10-target VHH campaign produces 135/4,085 passing designs with 60-fold pass rate variation | Computational (RFAntibody pipeline) |
| P4 | rfab-harness Platform (website) | Makani | 2026 | Open-source campaign orchestration for RFAntibody with automated analysis | Software tool |
| P5 | FLAb Analyses (website) | Makani | 2026 | Training data representation biases model performance on antibody structure prediction | Computational analysis |
| P6 | DADB Benchmark (website) | Makani | 2026 | Multi-metric gatekeeper benchmark for de novo antibody designs across structure, sequence, developability | Benchmark framework |
| P7 | 5-Platform Comparison (website) | Makani | 2026 | Cross-platform comparison reveals 1-2% experimental binding rates; only Latent-X2 reports immunogenicity | Comparative analysis |

---

## 3. Thematic Synthesis

### Theme A: Confidence Score Limitations --- The Core Challenge

The central finding of Smorodyina et al. (P1) is devastating for confidence-based filtering workflows: across 561,800 predictions per tool (106 real VHH-antigen complexes x 106 antigens x 50 replicates), none of the three state-of-the-art structure prediction tools could reliably separate cognate from shuffled pairings using their internal confidence scores.

**Quantitative evidence:**
- PR-AUC (Average Precision): AF3 = 0.187, Chai-1 = 0.067, Boltz-2 = 0.026 (random baseline ~0.011)
- Cross-tool score correlations: r=0.13 (Boltz-2/AF3), r=0.14 (Chai-1/Boltz-2), r=0.18 (Chai-1/AF3)
- ipTM score overlap: Boltz-2 assigns uniformly high ipTM (~0.85-0.91) to both real and shuffled complexes; AF3 shows bimodal distributions but with substantial overlap; Chai-1 yields uniformly low scores for both classes
- Epitope recall for shuffled complexes frequently matched or exceeded that of real complexes (Fig. 4A, P1)

**The "locked-in" hypothesis:** Smorodyina et al. provide mechanistic insight into why sampling does not rescue confidence-based filtering. While increased stochastic sampling (N=1 to N=100) improved best-case DockQ substantially (median improvement: AF3 +0.43, Boltz-2 +0.20, Chai-1 +0.15), the corresponding ipTM changes showed near-zero correlation with structural improvements (delta-DockQ vs delta-ipTM: r=-0.027 AF3, r=-0.040 Boltz-2, r=-0.019 Chai-1). Confidence scores are "locked in" to the initial MSA/embedding representation and do not update to reflect refinement achieved through diffusion sampling. This is an architectural limitation: ipTM derives from trunk network pairwise representations computed once from the MSA, remaining fixed regardless of downstream structure module iterations.

**Calibration failure modes differ by tool:**
- AF3: Best calibrated (DockQ vs ipTM r=0.736 at best-sample), but still 2% overconfident failures (Q2) and 28% underconfident successes (Q4)
- Boltz-2: Systematically overconfident (18% Q2 failures), calibration degrades catastrophically at single-sample (r=0.493)
- Chai-1: Systematically underconfident (24% Q4), confidence-based filtering would exclude many structurally accurate predictions

**Direct implication for rfab-harness:** Our pipeline uses pAE (predicted aligned error) as the primary filter (`filter.py:40-81`, threshold pAE < 10.0). While pAE is related to but distinct from ipTM (both derive from pairwise error predictions), Smorodyina et al.'s findings apply directly: pAE measures the model's self-assessed confidence in interface geometry, which P1 demonstrates is not calibrated to binding specificity. A shuffled (non-binding) VHH-antigen pair can achieve low pAE simply by adopting a geometrically plausible (but biologically meaningless) docking pose.

### Theme B: Data Scarcity and the ASD Opportunity

The ASD database (P2) addresses the data infrastructure gap that limits computational antibody model development. Key characteristics:

**Scale and composition:**
- 1,097,946 unique antibody-antigen interactions across 9,575 unique antigens
- 15 distinct sources compiled into 25 individual datasets
- Data types span structural (PDB-derived), sequence-based (GenBank), literature-extracted, and patent-derived records
- Schema includes complete antibody heavy/light sequences, antigen sequences, affinity measurements, UniProt/PDB metadata, RIOT-numbered CDR annotations, and confidence labels

**VHH-specific subset directly relevant to our work:**
- `structures-nanobodies`: 1,258 records, 390 unique antigens, 552 unique VHH antibodies (confidence: very_high)
- `structures-antibodies`: 2,711 records, 1,083 unique antigens, 1,327 unique antibodies (confidence: very_high)
- Combined structural subset: 3,969 records with experimental binding confirmation

**Critical limitations:**
- **HER2 dominance**: The "buzz" dataset alone contributes 524,346 records targeting a single antigen (HER2), representing ~48% of total ASD. Any model trained on raw ASD without stratification would be severely biased toward HER2 binding characteristics.
- **Affinity heterogeneity**: 48% of all entries use "fuzzy" affinity (categorical h/m/l from the buzz dataset). Remaining entries use diverse measurement types: KD (nM), IC50 (ug/ml), boolean binding, delta_g. Cross-dataset affinity comparison requires careful normalization.
- **Confidence stratification**: Four tiers (very_high, high, medium, low). Structural data is very_high; patent-derived data is medium due to chain-mapping uncertainties. This stratification is crucial for downstream use.

**Opportunity for rfab-harness validation:** The 1,258 VHH structures in ASD overlap significantly with the 106 complexes used in Smorodyina et al.'s benchmark (both drawn from PDB-deposited nanobody-antigen complexes). ASD's broader antigen coverage (390 unique antigens for VHH structures) could enable expanded negative-control generation for our 10 cancer targets. Specifically, ASD could provide:
1. Known binders to our target antigens (positive controls)
2. Known binders to unrelated antigens (shuffled negative controls, following P1's methodology)
3. Affinity-stratified training data for developing specificity-aware scoring

### Theme C: Developability and Multi-Metric Gaps

The convergence of P1's confidence score failure with our prior analyses (P5, P6, P7) exposes the narrow scope of current computational evaluation:

**FLAb training data bias (P5):** Models systematically perform better on targets present in their training data. In P1, Smorodyina et al. explicitly tracked train/test splits for each tool (AF3: 30 train / 76 test; Boltz-2: 64 train / 42 test; Chai-1: 25 train / 81 test) and found that discrimination failure persists across both splits. Train-leakage complexes showed slightly higher median epitope recall, consistent with partial memorization, but shuffled complexes spanned similar recall ranges regardless of leakage status. This confirms P5's finding that training data bias exists but does not explain the fundamental discrimination failure.

**DADB benchmark limitations (P6):** DADB's structure component (25% of total score) currently assesses structural plausibility using metrics from the same class (pAE, RMSD, pLDDT) that P1 demonstrates are non-discriminative. The gatekeeper concept --- requiring candidates to pass each metric category independently --- is sound architecturally, but the structure gate needs a specificity discrimination sub-score to avoid passing structurally plausible non-binders.

**5-platform comparison (P7):** The 1-2% experimental binding rate across platforms establishes the ground-truth hit rate for de novo antibody design. Our 3.3% computational pass rate (P3) exceeds this, consistent with P1's finding that confidence-based filtering inflates apparent success by passing structurally plausible non-binders. The fact that only Latent-X2 among 5 platforms reports immunogenicity data (P7) highlights that developability assessment remains an afterthought --- Smorodyina et al. confirm that confidence scores encode none of the drug-relevant properties (induced fit, off-target propensity, developability, conformational selection).

### Theme D: Benchmarking and Validation Standards

Smorodyina et al. establish a new gold standard for antibody-antigen prediction benchmarking:

**The shuffled-decoy benchmark design:**
- 106 experimentally solved VHH-antigen complexes from PDB (91 unique PDB IDs)
- All-vs-all pairing creates 11,236 shuffled non-cognate pairs per tool
- 50 stochastic replicates per pairing (total: ~1.8M structures across 3 tools)
- Train/test split tracking for each tool's training set cutoff date

**Key methodological innovations:**
1. **PR-AUC over ROC-AUC**: Under extreme class imbalance (~89 positives out of ~8,000+ evaluated pairs), PR-AUC is the appropriate metric. ROC-AUC would be misleadingly high due to the massive number of true negatives.
2. **Saturation sampling analysis**: Testing N=1 to N=100 samples reveals that most quality gain comes from N=10-25 (efficiency frontier), and that "just sample more" is expensive without improving discrimination.
3. **Cross-seed consensus**: Seed selection introduces DockQ variance of 0.04-0.05 (median), with 10-15% of systems strongly seed-sensitive. This suggests cross-seed consensus as a potential discrimination signal.

**ASD as a validation oracle:** ASD's 1,258 VHH-antigen structural records (P2) could serve as the expansion set for P1's benchmark. The 106-complex benchmark is powerful but small; ASD could enable scaling to ~500+ systems while maintaining the shuffled-decoy methodology. The antigen diversity in ASD (390 unique antigens for VHH structures vs. ~91 unique PDB antigens in P1) would also test whether discrimination failure generalizes across antigen families.

---

## 4. Impact on rfab-harness

### 4.1 Specific Code-Level Consequences

**`harness/analysis/filter.py:40-81` --- pAE/RMSD filtering:**
The `apply_filters()` function applies three thresholds:
- `pae_threshold: float = 10.0` --- pAE measures the model's predicted error in interface alignment. P1 shows this metric has near-zero discriminative power for binding specificity (PR-AUC 0.187 for the best tool, AF3).
- `rmsd_threshold: float = 2.0` --- CDR RMSD measures geometric deviation from the template-docked pose. P1's RMSD distributions (Fig. 4E) show similar values for real and shuffled complexes, indicating structural consistency does not imply correct pairing.
- `ddg_threshold: float | None = -20.0` --- Rosetta ddG is energy-based and partially orthogonal to confidence scores, but was not evaluated in P1's benchmark.

**`harness/analysis/rank.py:14-45` --- composite scoring:**
The `rank_candidates()` function computes: `0.4 * pae_norm + 0.3 * rmsd_norm + 0.3 * ddg_norm`. All three components measure structural self-consistency. None measures binding specificity, target selectivity, or developability. A shuffled (non-binding) VHH-antigen pair that adopts a geometrically clean pose would rank identically to a true binder.

### 4.2 Proposed Additions

1. **Shuffled controls per campaign**: Generate 100+ shuffled VHH-antigen pairings per target using ASD's nanobody pool. Score with the same pipeline. Any passing design that does not significantly outperform shuffled controls on at least one metric should be flagged.

2. **Cross-seed consensus scoring**: Run RF2 scoring with 5+ independent seeds per design. Designs with high variance across seeds are likely occupying shallow energy minima. P1 shows seed choice introduces DockQ variance of 0.04-0.05 median, with 10-15% of systems strongly seed-sensitive.

3. **ASD lookup for target coverage**: For each of our 10 target antigens, query ASD for existing known binders. This provides (a) positive controls for pipeline validation and (b) sequence/structural references for CDR comparison.

4. **Epitope recovery scoring**: Quantify what fraction of the intended epitope residues are actually contacted by the designed VHH in the RF2-predicted complex. P1 shows that high ipTM (and by extension low pAE) does not guarantee correct epitope placement.

5. **Specificity discrimination sub-score for DADB**: Add a component that measures whether a design's score distribution meaningfully separates from shuffled-control baselines.

### 4.3 Interpreting the 60-Fold Pass Rate Variation

Our preprint (P3) reports pass rates varying from 0.3% (EGFRvIII, HER2 Domain IV) to 19.8% (CEACAM5). In light of P1, this variation may partially reflect differential false-positive susceptibility rather than genuine design tractability:

- **CEACAM5 (19.8%)**: The concave A3B3 epitope may provide a "parking lot" geometry where many CDR conformations achieve low pAE/RMSD without necessarily achieving specific binding. This target's dominant pass rate could partly reflect geometric permissiveness rather than binding competence.
- **EGFRvIII (0.3%)**: The neoepitope junction's unusual geometry may simply be harder to dock plausibly (real or shuffled), resulting in low pass rates but potentially higher specificity among the rare passing designs.

This does not invalidate the rankings --- CEACAM5 candidates may still bind --- but it means pass rate is a compound metric reflecting both genuine design quality and geometric permissiveness.

---

## 5. Impact on DADB Benchmark

DADB's current structure component (25% of total score) uses pAE, RMSD, and pLDDT --- the same metric class that P1 demonstrates is non-discriminative for binding specificity. Specific changes needed:

1. **Add specificity discrimination sub-score**: For each candidate, generate a matched set of shuffled controls (same target, different nanobody sequences from ASD). The candidate's structural metrics must significantly exceed the shuffled-control distribution to earn specificity points.

2. **Reweight structure component**: The 25% allocation assumes structural metrics are partially informative. P1 suggests they are approximately random for discrimination (AP=0.187 vs baseline 0.011 for the best tool). Consider reducing structure weight or conditioning it on specificity discrimination.

3. **Add cross-seed consistency metric**: Run scoring with multiple seeds. Low variance across seeds indicates robust structural prediction, which P1 shows is weakly correlated with correct pairing but still provides marginal signal (Fig. 4E shows slightly higher RMSD for shuffled vs real complexes).

---

## 6. Actionable Recommendations

1. **Implement shuffled-control generation** in rfab-harness as a standard analysis step. For each campaign, generate N shuffled VHH-antigen pairs using ASD's nanobody pool, score them identically, and report the percentile rank of each design relative to shuffled controls.

2. **Add epitope recovery scoring** to `harness/analysis/`. Quantify the fraction of intended epitope residues contacted by each designed VHH in the RF2-predicted complex. Flag designs where the VHH docks to a non-epitope region despite low pAE.

3. **Integrate ASD lookup** into the target preparation step. For each target antigen, automatically query ASD (via UniProt ID or target name) to retrieve known binders as positive controls and unrelated nanobodies as negative control sources.

4. **Add cross-seed consensus** to the scoring pipeline. Run RF2 with 5 independent seeds per top-ranked design. Report score variance. Prioritize designs with consistent scores across seeds.

5. **Revise DADB structure component** to include a specificity discrimination sub-score based on shuffled-control comparison, reducing the weight of raw pAE/RMSD/pLDDT.

6. **Experimental validation priority**: Given P1's findings, prioritize experimental testing (yeast display, SPR) for the top 10 candidates rather than further computational optimization. P7's 1-2% experimental hit rate is the true benchmark; computational filtering beyond basic quality control may not improve this.

7. **Engage with Smorodyina et al.'s released dataset**: The ~1.8M AI-generated structures represent a community resource for developing specificity-aware scoring functions. Our 135 passing designs could be benchmarked against this dataset to assess whether they occupy the "real complex" or "shuffled complex" region of the structural landscape.

8. **Consider Rosetta ddG as a partially orthogonal signal**: Unlike pAE/ipTM (which derive from ML confidence), ddG captures energetic contributions. While not evaluated in P1, physics-based energy functions may encode specificity information that learned confidence scores miss. Test whether ddG discriminates shuffled from real complexes in P1's benchmark.

---

## References

- Smorodyina E, Ali M, Brumat K, et al. "Structural Plausibility Without Binding Specificity: Limits of AI-Based Antibody-Antigen Structure Prediction Confidence Scores." bioRxiv. 2026. doi:10.64898/2026.03.02.709004
- Czerwinski A, Dudzic P, Wojtowicz K, et al. "ASD: Antigen-Specific Antibody Database." bioRxiv. 2025. doi:10.1101/2025.11.24.690097
- Makani A. "Dogfooding rfab-harness: De Novo VHH Design for 10 Challenging Cancer Targets." Preprint, 2026.
- Watson JL, et al. "Atomically accurate de novo design of single-domain antibodies." Nature. 2025.
