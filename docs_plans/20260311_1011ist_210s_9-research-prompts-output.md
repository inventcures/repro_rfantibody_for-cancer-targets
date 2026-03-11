# 9-Step Literature Synthesis: 7-Work Corpus on AI Antibody Design Confidence and Data Infrastructure

> **Date**: 2026-03-11 | **Author**: Ashish Makani
> **Method**: `/synthesize-lit-review` framework (prompts #4-#12)

---

## Corpus

| ID | Work | Authors | Year |
|----|------|---------|------|
| P1 | Structural Plausibility Without Binding Specificity | Smorodyina, Greiff et al. | 2026 |
| P2 | ASD: Antigen-Specific Antibody Database | Czerwinski, Greiff, Krawczyk et al. | 2025 |
| P3 | Dogfooding rfab-harness on Cancer Targets (preprint) | Makani | 2026 |
| P4 | rfab-harness Platform (website/tool) | Makani | 2026 |
| P5 | FLAb Analyses (website) | Makani | 2026 |
| P6 | DADB Benchmark (website) | Makani | 2026 |
| P7 | 5-Platform Comparison (website) | Makani | 2026 |

---

## Prompt #4: Landscape Map

### Corpus Table

| ID | Author(s) | Year | Central Claim | Evidence Type |
|----|-----------|------|---------------|---------------|
| P1 | Smorodyina, Ali, Brumat, Salicari, Miklavc, Kappassov, Fu, Sormanni, de Marco, Greiff | 2026 | Confidence scores (ipTM) from AF3/Boltz-2/Chai-1 fail to discriminate cognate from shuffled VHH-antigen pairings (PR-AUC: AF3 AP=0.187, Boltz-2 AP=0.026) | Experimental benchmark: 106 complexes, 1.8M predictions |
| P2 | Czerwinski, Dudzic, Wojtowicz, Jaszczyszyn, Bielska, Wrobel, Demharter, Spreafico, Greiff, Krawczyk | 2025 | Unified 1.1M antibody-antigen interaction database enables model development through standardized formatting | Database: 15 sources, 9,575 unique antigens |
| P3 | Makani | 2026 | 10-target VHH campaign produces 135/4,085 passing designs (3.3%) with 60-fold pass rate variation (0.3%-19.8%) | Computational: RFAntibody pipeline on cancer targets |
| P4 | Makani | 2026 | Open-source campaign orchestration wrapping RFdiffusion-ProteinMPNN-RF2 with automated filtering/ranking | Software tool with 21 pre-built configs |
| P5 | Makani | 2026 | Training data representation biases structure prediction model performance on antibody targets | Computational analysis of FLAb benchmark |
| P6 | Makani | 2026 | Multi-metric gatekeeper benchmark (structure 25%, sequence, developability) for de novo antibody designs | Benchmark framework: DADB |
| P7 | Makani | 2026 | Cross-platform comparison of 5 de novo design tools reveals 1-2% experimental binding rates; only Latent-X2 reports immunogenicity | Comparative analysis |

### Thematic Clusters

**Cluster 1: Confidence Score Validity** (P1, P3, P5)
The central question: do computational confidence metrics reliably indicate binding? P1 answers definitively no for ipTM/pTM. P3 implicitly relies on the answer being yes (pAE < 10 filter). P5 shows training data bias contaminates even the metrics that do correlate with structural quality.

**Cluster 2: Data Infrastructure for Antibody AI** (P2, P5, P7)
The field lacks standardized, comprehensive antibody-antigen interaction datasets. P2 addresses this with ASD (1.1M interactions). P5 reveals that existing training data is biased. P7 shows experimental ground-truth across platforms is sparse (1-2% hit rates without systematic negative controls).

**Cluster 3: Pipeline Architecture and Filtering** (P3, P4, P6)
The generate-then-filter paradigm. P4 implements it (RFdiffusion-ProteinMPNN-RF2). P3 executes it on 10 cancer targets. P6 formalizes evaluation criteria. All three assume confidence metrics carry discriminative signal --- an assumption P1 undermines.

**Cluster 4: From Computation to Clinic** (P3, P6, P7)
Translational readiness. P3 produces candidates for experimental testing. P6 defines what "good enough" means. P7 establishes the experimental hit rate baseline (1-2%). None of the three addresses the gap between computational pass rates and experimental success rates with mechanistic rigor.

**Inter-work contradiction flags:** P3's 3.3% pass rate vs P7's 1-2% experimental hit rate; P3's pAE filter vs P1's pAE discrimination failure; P6's structure component vs P1's finding that structural metrics are non-discriminative.

---

## Prompt #5: Contradiction Finder

| Claim A (Source) | vs Claim B (Source) | Nature | Why They Conflict | Implications |
|------------------|---------------------|--------|-------------------|--------------|
| pAE < 10 filter selects for design quality (P3, `filter.py:40-81`) | pAE/ipTM fails to discriminate cognate from shuffled pairings; PR-AUC only 0.187 for best tool (P1, Fig. 2B) | Empirical | P3 uses pAE as a quality gate assuming lower pAE = better interface. P1 shows this metric assigns similar values to real and shuffled complexes. Different test: P3 tests self-consistency of de novo designs, P1 tests discrimination of known binders from decoys. | The 135 designs passing P3's filter may include structurally plausible non-binders. pAE threshold selects for geometric compatibility, not binding specificity. Resolving requires experimental testing of the 135 candidates. |
| 3.3% computational pass rate suggests appropriately stringent filtering (P3, preprint Discussion) | 1-2% experimental binding rate across platforms is the ground truth (P7) | Methodological | P3's 3.3% is computed from structural metrics. P7's 1-2% is from experimental binding assays. The computational rate exceeds the experimental rate, suggesting the filter is too permissive rather than too stringent --- consistent with P1's false-positive inflation finding. | Computational pass rate is not directly comparable to experimental binding rate, but the gap (3.3% vs 1-2%) is directionally consistent with P1's conclusion that confidence-based filtering inflates apparent success. |
| ASD provides 48% "fuzzy" (categorical h/m/l) affinity data suitable for model training (P2, Methods) | DADB benchmark requires quantitative scoring for each design property (P6) | Definitional | P2's "fuzzy" affinity is a categorical classification (high/medium/low) from the buzz dataset's binning. P6 requires numerical scores. Nearly half of ASD cannot directly feed DADB's scoring pipeline without imputation or re-categorization. | Using ASD for DADB validation requires either restricting to the 52% with numerical affinity or developing a mapping from categorical to quantitative scores. |
| Training data bias explains performance differences across targets (P5) | Discrimination failure persists equally for train-leakage and true-test complexes (P1, Fig. 2C) | Empirical | P5 finds that models perform better on targets in their training set. P1 finds that even train-leakage complexes (where both VHH and antigen were in training data) are not reliably discriminated from shuffled controls. Training data presence may improve structural accuracy but not binding discrimination. | Training bias affects structural quality (P5 is correct) but not specificity discrimination (P1 extends the finding). The two claims are compatible but address different failure modes: P5 addresses quality, P1 addresses specificity. Resolution: training data bias and discrimination failure are orthogonal problems. |
| Backbone topology is the primary determinant of design success in RFAntibody (P3, Discussion) | Seed selection determines qualitative docking mode; sampling refines within that mode (P1, Discussion) | Methodological | P3 attributes pass rate variation to backbone shape complementarity with the epitope. P1 attributes structural outcome to initial seed/path selection. Both describe the same phenomenon from different perspectives: the initial structural hypothesis (backbone or seed) determines outcome more than downstream refinement. | These claims are actually convergent. P3's backbone topology and P1's seed selection both describe the dominance of initial conditions over iterative refinement. Implication: invest in diverse initial conditions (more backbones, more seeds) rather than deeper sampling per condition. |

---

## Prompt #6: Citation Chain

### Concept 1: ipTM as a Discriminative Metric

**Intellectual lineage through the corpus:**
- **Origin**: ipTM (interface predicted Template Modeling score) was introduced as a confidence metric for multi-chain predictions in AlphaFold-Multimer (refs 35, 42 in P1). The field adopted it as the standard ranking metric for antibody-antigen complex predictions.
- **Adoption (P4, P3)**: rfab-harness uses the related pAE metric (which shares the same trunk-network pairwise representation as ipTM) as its primary filter. The composite score in `rank.py` weights pAE at 0.4 --- the highest weight among all metrics.
- **Challenge (P1)**: Smorodyina et al. directly test whether ipTM can discriminate cognate from shuffled pairings. Finding: it cannot (AP=0.187 for AF3, where random baseline is ~0.011). ipTM captures geometric plausibility, not biological specificity.
- **Consequence (P6, P7)**: DADB's structure component relies on the same metric class. P7's platform comparison lacks specificity controls entirely.
- **Concept drift**: ipTM was designed as a structural quality indicator (does the predicted structure match the true structure?). The field drifted to using it as a binding indicator (does this antibody bind this antigen?). P1 exposes this drift as invalid.

### Concept 2: Data Scarcity in Antibody-Antigen Interactions

**Intellectual lineage:**
- **Origin**: Long-standing problem in the field. SAbDab and PLAbDab attempted to address it but require significant preprocessing (P2, Table 1). Purpose-built affinity datasets revealed that data curation, not model architecture, is the bottleneck (refs 17-19 in P2).
- **Quantification (P2)**: ASD aggregates 1,097,946 interactions from 15 sources, demonstrating the scale achievable through systematic aggregation. But reveals new problems: HER2 dominance (524K/1.1M), heterogeneous affinity measures (48% fuzzy).
- **Interaction with bias (P5)**: FLAb analyses show that training data representation directly biases model performance. Models have seen more HER2 complexes (the dominant target in ASD) than rare cancer targets like GPC2 or EphA2.
- **Benchmarking consequence (P1, P7)**: P1's benchmark uses only 106 VHH-antigen complexes --- a tiny fraction of what ASD makes available --- because the shuffled-decoy design requires experimentally confirmed structures. P7's 1-2% hit rate is established without systematic negative controls from a database like ASD.
- **Gap**: No work in this corpus connects ASD's 1,258 VHH structures directly to P1's benchmarking methodology to scale up the discrimination test.

### Concept 3: De Novo Antibody Design Success Rates

**Intellectual lineage:**
- **Origin**: Watson et al. (Nature 2025, referenced in P3, P4) reported per-target experimental binding rates of 0-2% for RFAntibody designs.
- **Computational proxy (P3)**: 3.3% computational pass rate across 10 cancer targets, with 60-fold variation (0.3%-19.8%). Interpreted as reflecting genuine design difficulty variation.
- **Experimental ground truth (P7)**: 1-2% across 5 platforms. Establishes that the field-wide experimental success rate has not substantially improved despite different computational approaches.
- **Reinterpretation (P1)**: The gap between computational pass rates and experimental binding rates is at least partially explained by false-positive inflation from non-discriminative confidence scores. "Structurally plausible" does not equal "binding."
- **Benchmark integration (P6)**: DADB attempts to formalize what constitutes a successful design, but its structure component is now undermined by P1's findings.

---

## Prompt #7: Five Unanswered Research Questions

### Q1: Does Rosetta ddG discriminate shuffled from real VHH-antigen complexes?

**Urgency**: P1 tested ipTM/pTM/clash scores but not physics-based energy functions. P3's pipeline includes ddG (weight 0.3 in composite score). ddG captures van der Waals, electrostatic, and solvation contributions that ML confidence scores do not.
**Why gap exists**: P1 focused on ML structure prediction tools; Rosetta ddG is a separate physics-based scoring function that operates on predicted structures rather than predicting structures itself.
**Closest paper**: P1 (has the benchmark framework and structures)
**Methodology needed**: Run Rosetta's InterfaceAnalyzer on P1's 1.8M structures (both cognate and shuffled), compute ddG, and calculate PR-AUC for ddG-based discrimination.
**Answerability**: 1 (can be done immediately with existing data and tools)

### Q2: What is the experimental binding rate of rfab-harness designs passing pAE < 10 and RMSD < 2.0?

**Urgency**: P3 reports 135 passing designs. P1 suggests many may be false positives. P7 establishes 1-2% as the field-wide experimental rate. The actual binding rate of our specific 135 candidates is unknown.
**Why gap exists**: No experimental validation has been performed. The pipeline is purely computational.
**Closest paper**: P3 (has the designs), P7 (has the experimental baseline)
**Methodology needed**: Yeast surface display or SPR binding assays for the top 20-50 candidates. Compare observed binding rate to the 3.3% computational pass rate and P7's 1-2% baseline.
**Answerability**: 2 (requires wet lab work but is straightforward experimentally)

### Q3: Can cross-seed consensus scoring improve binding discrimination beyond single-seed confidence?

**Urgency**: P1 shows seed selection introduces DockQ variance of 0.04-0.05 (median) with 10-15% of systems strongly seed-sensitive. P1 also shows confidence scores are "locked in" to initial conditions. Cross-seed agreement might indicate robust binding vs. artifactual docking.
**Why gap exists**: P1 analyzed seeds for quality improvement but did not test whether cross-seed consistency discriminates real from shuffled.
**Closest paper**: P1 (has multi-seed data for all 106 systems x 3 tools)
**Methodology needed**: For each VHH-antigen pair in P1's benchmark, compute the variance of ipTM/DockQ across seeds. Test whether real complexes show lower variance (more consistent predictions) than shuffled complexes. Compute PR-AUC for variance-based discrimination.
**Answerability**: 1 (data already exists in P1's released dataset)

### Q4: How does ASD's antigen diversity correlate with prediction difficulty across cancer targets?

**Urgency**: P3 shows 60-fold pass rate variation but attributes it to epitope geometry. P2's ASD has variable representation across antigens. Targets with more ASD entries (more training exposure) might show higher pass rates due to model familiarity rather than genuine tractability.
**Why gap exists**: No work has linked ASD's antigen representation to computational design success rates.
**Closest paper**: P3 (has per-target pass rates), P2 (has per-antigen record counts)
**Methodology needed**: For each of P3's 10 targets, count matching ASD records. Correlate with pass rate. Confound with epitope geometry (convexity, area). Requires multivariate analysis.
**Answerability**: 2 (data exists but analysis requires careful confound control)

### Q5: Can a specificity-aware scoring function be trained on ASD + P1's shuffled benchmark?

**Urgency**: P1 demonstrates the problem (confidence scores are non-discriminative). P2 provides the data (1.1M interactions with affinity labels). Together they define the training signal: positive = ASD confirmed binders, negative = P1-style shuffled pairings.
**Why gap exists**: P1 focused on benchmarking existing tools, not building new scoring functions. P2 focused on data aggregation, not model development. The synthesis of both papers' contributions is unexplored.
**Closest paper**: P1 + P2 combined
**Methodology needed**: Train a binary classifier (binding vs shuffled) using structural features from predicted complexes (contact patterns, interface geometry, buried surface area) on ASD's structural subset + shuffled controls. Test generalization on held-out antigens.
**Answerability**: 3 (requires significant ML development and careful validation)

---

## Prompt #8: Methodology Comparison Matrix

| Aspect | P1 (Smorodyina) | P2 (ASD) | P3 (Dogfooding) | P4 (rfab-harness) | P5 (FLAb) | P6 (DADB) | P7 (5-Platform) |
|--------|-----------------|----------|-----------------|-------------------|-----------|-----------|-----------------|
| **Data source** | PDB (106 VHH-Ag complexes) | 15 sources (PDB, GenBank, patents, literature) | Modal cloud (10 cancer targets) | N/A (tool) | PDB + model predictions | Aggregated benchmarks | Published platform reports |
| **Sample size** | 106 real + 11,236 shuffled pairs; 1.8M predictions | 1,097,946 interactions | 4,085 designs, 135 passing | 21 pre-built configs | Variable per analysis | Variable per benchmark | 5 platforms |
| **Validation approach** | Shuffled-decoy discrimination (PR-AUC) | Overlap with SAbDab/PLAbDab; confidence labels | Computational scoring only (pAE, RMSD) | Unit tests + smoke campaigns | Cross-reference with experimental data | Multi-metric gating | Literature-reported experimental rates |
| **Statistical methods** | PR-AUC, Pearson/Spearman correlation, binomial enrichment test, Wilcoxon signed-rank | Deduplication, intraset aggregation | Descriptive statistics, distributions | N/A | Comparative statistics | Weighted scoring | Qualitative comparison |
| **Computational tools** | AF3, Boltz-1/2, Chai-1, TopModel (clash), DockQ | Claude Opus 4.1 (OCR), custom ETL pipeline | RFdiffusion, ProteinMPNN, RF2 | rfab-harness wrapping RFAntibody | Various structure prediction tools | Multiple tools per component | Various (RFAntibody, EvolvePro, etc.) |
| **Key assumptions** | Shuffled pairs are non-binders; off-target binding <1% | Source datasets are correctly annotated; OCR is accurate | pAE < 10 and RMSD < 2.0 indicate design quality | RFAntibody pipeline is correctly implemented | Training set dates are accurately tracked | Multi-metric gating captures design quality | Published experimental rates are comparable |
| **Reproducibility** | High (code + 1.8M structures released) | High (database publicly available at naturalantibody.com/agab/) | High (code + configs open-source) | High (open-source tool) | Medium (analyses partially documented) | Medium (framework defined, implementation varies) | Low (relies on published reports with different methodologies) |

### Dominant methodology: Computational structure prediction scoring
Five of seven works (P1, P3, P4, P5, P6) rely on or evaluate computational structure prediction. This creates a field bias toward optimizing for metrics that P1 shows are non-discriminative.

### Underused approaches:
- **Experimental binding assays**: Only P7 references experimental data, and indirectly through platform reports. No work in this corpus directly validates computational predictions with binding experiments.
- **Physics-based energy evaluation**: ddG/Rosetta energy functions are used in P3/P4 but not evaluated in P1's discrimination benchmark.
- **Molecular dynamics**: P1's Discussion recommends MD-informed filtering but no work implements it.

### Weakest methodological link:
P7's cross-platform comparison relies on published reports with different experimental protocols, different target panels, and different success criteria. The 1-2% hit rate is an approximate range, not a controlled measurement. This makes P7's ground-truth baseline less reliable than it appears.

---

## Prompt #9: 400-Word Field Synthesis

The field of computational antibody design has arrived at a paradox. Structure prediction tools can generate geometrically plausible antibody-antigen complexes at scale --- RFdiffusion produces thousands of backbone scaffolds per target, ProteinMPNN sequences them, and scoring tools like RF2, AlphaFold3, Boltz-2, and Chai-1 assign confidence scores that look reassuringly quantitative (P3, P4). The implicit contract is that higher confidence means higher likelihood of binding. Smorodyina et al. (P1) have broken this contract. Their controlled benchmark of 106 experimentally determined nanobody-antigen complexes against 11,236 shuffled decoys demonstrates that none of the three leading structure prediction tools can reliably distinguish real from fabricated pairings using their internal confidence scores. The best performer, AlphaFold3, achieves a PR-AUC of 0.187 against a random baseline of 0.011 --- above chance but far from useful for screening. Boltz-2 performs at 0.026, barely above noise.

This finding cascades through the entire computational antibody pipeline. The rfab-harness platform (P4), applied to 10 cancer targets (P3), filters 4,085 designs down to 135 using pAE and RMSD thresholds --- metrics from the same family that P1 demonstrates are non-discriminative. The 3.3% computational pass rate exceeds the 1-2% experimental binding rate observed across five design platforms (P7), consistent with confidence-based filtering inflating apparent success by passing structurally plausible non-binders.

Meanwhile, the data infrastructure has expanded. The ASD database (P2) aggregates 1,097,946 antibody-antigen interactions including 1,258 nanobody-specific structural records, providing the raw material for developing better discrimination metrics. But ASD has its own limitations: 48% of entries carry only categorical affinity and a single antigen (HER2) dominates nearly half the database. FLAb analyses (P5) confirm that training data bias systematically skews model performance, and the DADB benchmark (P6) uses structural metrics now shown to be non-discriminative for its structure component.

What the field collectively believes: generating structurally plausible complexes is a solved problem. What is genuinely contested: whether any existing computational metric can predict actual binding. What has been proven: structural plausibility and binding specificity are distinct properties, and current tools optimize the former without guaranteeing the latter. The single most important open question is whether physics-based energy functions, cross-seed consensus, or learned specificity classifiers --- trained on databases like ASD against shuffled controls --- can close the discrimination gap that confidence scores cannot.

---

## Prompt #10: Untested Assumptions

### P1 (Smorodyina et al.)
**Assumption**: Shuffled VHH-antigen pairs are non-binders.
**Relied on by**: P1's entire benchmark validity (shuffled pairs define the negative class).
**Consequence if wrong**: If even 1-5% of shuffled pairs bind (off-target effect, estimated <1% in P1's Discussion, ref 66), the negative class is contaminated and true PR-AUC is higher than reported. However, even correcting for this, the gap between AF3's AP=0.187 and useful discrimination (>0.5) remains enormous.
**How to test**: Experimentally test a random sample of 50 high-confidence shuffled complexes for binding (SPR or BLI).

### P2 (ASD)
**Assumption**: Antibody-antigen interactions labeled as "binding" in source datasets are correctly annotated.
**Relied on by**: P2's entire database integrity.
**Consequence if wrong**: Incorrectly annotated binders would propagate through any model trained on ASD. Patent-derived data (217K records, "medium" confidence) is particularly vulnerable --- the mapping assumes each patent family has only one heavy and one light chain.
**How to test**: Cross-validate a random subset of ASD entries against orthogonal binding databases (SAbDab structures, published SPR data).

### P3 (Dogfooding preprint)
**Assumption**: The NbBCII10 VHH framework is equally compatible with all 10 target epitopes.
**Relied on by**: P3's cross-target pass rate comparison.
**Consequence if wrong**: Some targets may fail not because the epitope is difficult but because the framework scaffold is geometrically incompatible. The 60-fold pass rate variation could partially reflect framework-epitope mismatch rather than intrinsic design difficulty.
**How to test**: Repeat campaigns with 2-3 alternative VHH frameworks (e.g., cAbBCII10, Nb.b201) and compare pass rate variation.

### P4 (rfab-harness)
**Assumption**: The RFAntibody pipeline (RFdiffusion-ProteinMPNN-RF2) is correctly implemented in the harness wrapper.
**Relied on by**: All results from P3.
**Consequence if wrong**: Implementation bugs (incorrect epitope definition, wrong truncation, misaligned frameworks) would invalidate all 135 passing candidates.
**How to test**: Reproduce Watson et al.'s published benchmark results (IL-7Ra, TrkA) through the harness and compare to their reported metrics.

### P5 (FLAb)
**Assumption**: Model training set cutoff dates are accurately documented.
**Relied on by**: P5's train/test split analyses, P1's leakage analysis.
**Consequence if wrong**: Misattributed train/test status would confound the analysis of whether training data bias explains performance differences.
**How to test**: Verify training set composition for AF3/Boltz-2/Chai-1 against their official documentation and model cards.

### P6 (DADB)
**Assumption**: The 25% weight for structure component reflects the relative importance of structural quality to overall design quality.
**Relied on by**: DADB's composite scoring.
**Consequence if wrong**: Given P1's findings, 25% may be too generous. A non-discriminative component receiving 25% weight dilutes the signal from potentially informative components (sequence quality, developability).
**How to test**: Correlate DADB composite scores with experimental binding data (when available) and test whether removing or downweighting the structure component improves predictive accuracy.

### P7 (5-Platform Comparison)
**Assumption**: Published experimental binding rates from different platforms are comparable.
**Relied on by**: P7's 1-2% hit rate baseline.
**Consequence if wrong**: Different platforms test different numbers of candidates, use different binding assays (yeast display, SPR, ELISA), and apply different success thresholds. The "1-2%" range may not reflect true methodological comparability.
**How to test**: Standardize: test the same set of designed sequences from multiple platforms against the same target using the same assay protocol.

---

## Prompt #11: Knowledge Map

```
Central Claim: AI structure prediction generates structurally plausible
but not necessarily binding-specific antibody-antigen complexes
  |
  |-- Pillar 1: Confidence scores fail to discriminate (P1)
  |     |-- PR-AUC: AF3=0.187, Chai-1=0.067, Boltz-2=0.026 (P1, Fig. 2B)
  |     |-- Cross-tool agreement is weak: r=0.13-0.18 (P1, Fig. 2A)
  |     |-- ipTM "locked-in" to initial conditions: delta-ipTM vs delta-DockQ r~0 (P1, Fig. 3C)
  |     |-- Epitope recall similar for real and shuffled complexes (P1, Fig. 4A-B)
  |
  |-- Pillar 2: Existing pipelines assume confidence = quality (P3, P4)
  |     |-- filter.py: pAE < 10.0, RMSD < 2.0 (P4)
  |     |-- rank.py: 0.4*pAE + 0.3*RMSD + 0.3*ddG (P4)
  |     |-- 135/4,085 pass (3.3%), 60-fold variation across targets (P3)
  |     |-- CEACAM5 dominance (19.8%) may reflect geometric permissiveness (P3 + P1)
  |
  |-- Pillar 3: Data infrastructure is growing but flawed (P2, P5)
  |     |-- ASD: 1,097,946 interactions, 9,575 antigens (P2)
  |     |-- HER2 dominance: 524K/1.1M records (P2)
  |     |-- 48% fuzzy affinity limits quantitative use (P2)
  |     |-- Training data bias contaminates predictions (P5)
  |     |-- VHH structural subset: 1,258 records enables nanobody-specific work (P2)
  |
  |-- Pillar 4: Experimental ground truth remains sparse (P7)
  |     |-- 1-2% binding rate across 5 platforms (P7)
  |     |-- Only Latent-X2 reports immunogenicity (P7)
  |     |-- No platform systematically generates negative controls (P7)
  |
  |-- Contested Zone: Pass rate variation --- intrinsic difficulty vs. false-positive susceptibility
  |     |-- P3: 60-fold variation reflects genuine epitope geometry differences
  |     |-- P1: Variation may partially reflect differential geometric permissiveness
  |     |-- Resolution requires: experimental validation of top candidates per target
  |
  |-- Contested Zone: Training data bias vs. fundamental discrimination failure
  |     |-- P5: Training bias explains some performance variation
  |     |-- P1: Discrimination failure persists across train/test splits
  |     |-- Resolution: Orthogonal problems; both contribute but to different aspects
  |
  |-- Frontier Question 1: Can physics-based ddG or cross-seed consensus discriminate?
  |-- Frontier Question 2: Can ASD + shuffled decoys train a specificity classifier?
  |-- Frontier Question 3: What is the actual experimental binding rate of rfab-harness designs?
```

### 3 Must-Read Works

1. **P1 (Smorodyina et al.)**: The foundational challenge. Every computational antibody design group using confidence-based filtering must reckon with this paper's evidence that their primary metric is non-discriminative.

2. **P2 (ASD Database)**: The data infrastructure. With 1.1M interactions and a VHH-specific subset of 1,258 structures, ASD is the largest publicly available resource for training and validating antibody-antigen interaction models.

3. **P3 (Dogfooding preprint)**: The practical case study. Shows exactly how confidence-based filtering plays out on real cancer targets, providing the concrete pipeline and candidate set that P1's findings reinterpret.

---

## Prompt #12: "So What" Test

### One sentence that is proven:

Confidence scores from current AI structure prediction tools (ipTM, pTM, pAE) measure structural self-consistency of predicted antibody-antigen complexes, not binding specificity, and cannot reliably distinguish designed binders from geometrically plausible non-binders.

### One honest admission:

We do not know how many of our 135 computationally passing VHH designs for 10 cancer targets actually bind their intended antigens --- and the best available evidence (P1's PR-AUC data, P7's 1-2% experimental rate) suggests the answer may be closer to 2-5 than 135.

### One real-world implication:

Any drug discovery team using AI-predicted confidence scores (pAE, ipTM, pLDDT) as the primary filter for advancing antibody candidates to experimental testing should immediately implement shuffled-control baselines: generate and score 100+ non-cognate pairings per target using the same pipeline, and only advance candidates that statistically outperform the shuffled distribution --- because without this control, they are selecting for structural plausibility, not binding.
