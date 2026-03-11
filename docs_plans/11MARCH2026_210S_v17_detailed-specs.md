# Detailed Specifications: Literature Synthesis --- Smorodyina/Greiff + ASD + Prior Work

> **Version**: 17 | **Date**: 2026-03-11 | **Status**: IMPLEMENTED
> **Session**: 210s (10:11 IST)

---

## 1. Overview

### Goal

Synthesize a 7-work corpus spanning two new Greiff-group papers (Smorodyina et al. 2026, ASD/Czerwinski et al. 2025) with 5 prior project analyses into a structured set of deliverables that (a) expose confidence score limitations in our rfab-harness pipeline, (b) identify actionable improvements, and (c) create a reusable skill for future literature syntheses.

### Key Insight Chain

1. Smorodyina et al. prove ipTM/pTM/pAE fail to discriminate real from shuffled VHH-antigen pairings (PR-AUC: AF3 AP=0.187)
2. Our rfab-harness uses pAE < 10 + RMSD < 2.0 as primary filters (`filter.py:40-81`)
3. These filters select for structural self-consistency, not binding specificity
4. 135 passing designs (3.3%) may include structurally plausible non-binders
5. ASD provides 1.1M interactions (including 1,258 VHH structures) to build shuffled controls
6. Proposed remedy: shuffled-control baselines, cross-seed consensus, ASD lookup, epitope recovery scoring

### Deliverable Summary

| # | Deliverable | Status |
|---|-------------|--------|
| 1 | Custom skill `/synthesize-lit-review` | DONE |
| 2 | Deep research results file (Part A) | DONE |
| 3 | 9 research prompts output (Part B) | DONE |
| 4 | This specifications document | DONE |

---

## 2. Corpus Definition

### The 7-Work Corpus

| ID | Title | Authors | Year | Type | Key Finding |
|----|-------|---------|------|------|-------------|
| P1 | Structural Plausibility Without Binding Specificity | Smorodyina, Ali, Brumat, Salicari, Miklavc, Kappassov, Fu, Sormanni, de Marco, Greiff | 2026 | bioRxiv preprint | ipTM fails discrimination: PR-AUC AF3=0.187, Boltz-2=0.026, Chai-1=0.067 |
| P2 | ASD: Antigen-Specific Antibody Database | Czerwinski, Dudzic, Wojtowicz, Jaszczyszyn, Bielska, Wrobel, Demharter, Spreafico, Greiff, Krawczyk | 2025 | bioRxiv preprint | 1,097,946 interactions, 9,575 antigens, 15 sources, 1,258 VHH structures |
| P3 | Dogfooding rfab-harness on Cancer Targets | Makani | 2026 | Preprint (v3) | 135/4,085 passing (3.3%), 60-fold pass rate variation across 10 targets |
| P4 | rfab-harness Platform | Makani | 2026 | Open-source tool | Campaign orchestration wrapping RFdiffusion-ProteinMPNN-RF2 |
| P5 | FLAb Analyses | Makani | 2026 | Website | Training data representation biases structure prediction performance |
| P6 | DADB Benchmark | Makani | 2026 | Website | Multi-metric gatekeeper: structure (25%), sequence, developability |
| P7 | 5-Platform Comparison | Makani | 2026 | Website | 1-2% experimental binding rate; only Latent-X2 has immunogenicity data |

### Source Papers (DOIs)

- P1: `10.64898/2026.03.02.709004` (bioRxiv, posted March 3, 2026)
- P2: `10.1101/2025.11.24.690097` (bioRxiv, posted November 26, 2025)

### Source PDFs (Local Paths)

- P1: `papers/11march2026_evo-smordoina/Structural Plausibility Without Binding Specificity [...].pdf`
- P2: `papers/11march2026_evo-smordoina/ASD: Antigen-Specific Antibody Database [...].pdf`

---

## 3. DR Results File Specification

**File**: `docs_plans/20260311_1011ist_210s_2papers-new-datassets_abag-db+developability_victor-grieff_eva-smordina_____dr-detailed-results.md`

### Section Outline

| Section | Description | Target Words |
|---------|-------------|--------------|
| 1. Executive Summary | Convergence of confidence failure + data infrastructure + prior work | ~400 |
| 2. Corpus Table | All 7 works: ID, authors, year, central claim, evidence type | ~200 |
| 3. Thematic Synthesis | 4 themes (A-D) with quantitative evidence | ~1,600 |
| 4. Impact on rfab-harness | Code-level consequences + proposed additions | ~500 |
| 5. Impact on DADB | Structure component revision | ~200 |
| 6. Actionable Recommendations | 8 numbered items | ~400 |
| **Total** | | **~3,300** |

### Theme Specifications

| Theme | Scope | Required Evidence |
|-------|-------|-------------------|
| A: Confidence Score Limitations | P1's PR-AUC results, cross-tool correlations, locked-in effect | AF3 AP=0.187, cross-tool r=0.13-0.18, delta-DockQ vs delta-ipTM r~0 |
| B: Data Scarcity & ASD | P2's scale, HER2 dominance, VHH subset, affinity heterogeneity | 1,097,946 total, 524K buzz, 48% fuzzy, 1,258 VHH structures |
| C: Developability Gaps | P5 training bias, P6 gatekeeper, P7 immunogenicity gap | FLAb bias finding, DADB 25% structure weight, Latent-X2 only |
| D: Benchmarking Standards | P1's shuffled-decoy design, ASD as validation oracle | 106 complexes, 1.8M structures released, 50 replicates |

### Code References Required

| Reference | Context |
|-----------|---------|
| `harness/analysis/filter.py:40-81` | pAE < 10.0, RMSD < 2.0, ddG < -20.0 thresholds |
| `harness/analysis/rank.py:14-45` | Composite: 0.4*pAE + 0.3*RMSD + 0.3*ddG |

---

## 4. Research Prompts Specification

**File**: `docs_plans/20260311_1011ist_210s_9-research-prompts-output.md`

### Input/Output Spec per Prompt

| Prompt # | Name | Input | Output Format | Key Requirement |
|----------|------|-------|---------------|-----------------|
| 4 | Landscape Map | 7 works | Table + 3-5 clusters + contradiction flags | All 7 works classified; inter-work contradictions identified |
| 5 | Contradiction Finder | Corpus analysis | Table: Claim A vs Claim B, Nature, Why, Implications | Min 4 contradictions with specific findings (numbers, figure refs) |
| 6 | Citation Chain | 3 recurring concepts | Lineage traces through all works per concept | 3 concepts: ipTM, data scarcity, de novo success rates |
| 7 | 5 Unanswered Questions | Corpus gaps | Each: question, 2+ citations, gap reason, closest paper, methodology, answerability rank | Ranked 1-5 by answerability |
| 8 | Methodology Matrix | All 7 works | Aspects x Works table + dominant/underused/weakest flags | 7+ aspects, all 7 works |
| 9 | 400-Word Synthesis | Full corpus | Continuous prose, no bullets/tables | All 7 works referenced, 380-420 words, no summaries |
| 10 | Untested Assumptions | Per-work analysis | 1-2 assumptions per work with consequences | All 7 works covered, testability described |
| 11 | Knowledge Map | Hierarchical synthesis | ASCII tree with pillars, contested zones, frontier questions, 3 must-reads | Central claim + 4+ pillars + 2+ contested zones + 3 must-reads |
| 12 | "So What" Test | Full corpus | 3 deliverables (proven, admission, implication) | Each exactly one statement |

### Cross-Reference Requirements

- Prompt #5 must cite specific findings from P1 (PR-AUC values, figure numbers) and P3 (filter thresholds, pass rates)
- Prompt #6 must trace each concept through all 7 works where relevant
- Prompt #9 must reference all 7 works by ID at least once
- Prompt #11 must include at least 2 contested zones

---

## 5. Custom Skill Specification

**File**: `/home/tp53/.claude/skills/synthesize-lit-review/SKILL.md`

### Format

```yaml
---
name: synthesize-lit-review
description: Run a 9-step structured literature synthesis on a corpus of papers/works. Usage: /synthesize-lit-review [optional corpus description or file path]
---
```

Followed by markdown body with:
- Workflow (4 steps: identify corpus, build table, execute 9 steps, output)
- The 9 Analysis Steps (#4-#12) with detailed instructions per step
- Tools section
- Rules section

### Invocation

```
/synthesize-lit-review
```

User provides corpus description or file path. Skill parses input, builds corpus table, and executes all 9 steps sequentially.

### Verification

The skill appears in the available skills list after creation. Confirmed: `synthesize-lit-review` appears in the skill listing.

---

## 6. Implementation Sequence

```
Step 1: Create skill file (no dependencies)
   |
   v
Step 2: Write DR results (Part A)
   |   Requires: Both paper PDFs read, filter.py + rank.py read,
   |             preprint tex read, specs format reference
   v
Step 3: Write 9 prompts output (Part B)
   |   Requires: Part A complete (shares corpus analysis)
   v
Step 4: Write specs file (this document)
        Requires: Steps 1-3 complete (references all deliverables)
```

All steps executed sequentially in a single session.

---

## 7. File Manifest

| # | File Path | Type | Status |
|---|-----------|------|--------|
| 1 | `~/.claude/skills/synthesize-lit-review/SKILL.md` | Custom skill | CREATED |
| 2 | `docs_plans/20260311_1011ist_210s_2papers-new-datassets_abag-db+developability_victor-grieff_eva-smordina_____dr-detailed-results.md` | Deep research (Part A) | CREATED |
| 3 | `docs_plans/20260311_1011ist_210s_9-research-prompts-output.md` | 9 prompts output (Part B) | CREATED |
| 4 | `docs_plans/11MARCH2026_210S_v17_detailed-specs.md` | This specs file | CREATED |

### Referenced Files (Not Modified)

| File | Role |
|------|------|
| `harness/analysis/filter.py` (lines 40-81) | pAE/RMSD/ddG filtering --- challenged by P1 |
| `harness/analysis/rank.py` (lines 14-45) | Composite scoring --- all metrics structural self-consistency |
| `~/.claude/skills/journal/SKILL.md` | SKILL.md format reference |
| `docs_plans/v1_challenging-but-validated+impactful-targets-in-cancer_detailed_specs.md` | Existing specs format reference |
| `out/v3_scientific_preprint.tex` | Preprint with 10-target campaign results |

---

## 8. Verification Checklist

### Skill Verification
- [x] `/synthesize-lit-review` appears in available skills list
- [x] YAML frontmatter has `name` and `description` fields
- [x] All 9 steps (#4-#12) defined with clear instructions
- [x] Tools section specifies Read, Write, Edit, Glob, Grep, Agent
- [x] Rules section includes no-fabrication and consistency requirements

### DR File (Part A) Verification
- [x] All 7 works appear in corpus table with correct attribution
- [x] No fabricated citations (all PR-AUC values, r-values, thresholds traced to paper content)
- [x] Specific code paths cited: `filter.py:40-81`, `rank.py:14-45`
- [x] 4 themes cover confidence scores, data infrastructure, developability, benchmarking
- [x] Actionable recommendations are numbered and specific

### 9 Prompts (Part B) Verification
- [x] All 9 prompts (#4-#12) executed with self-contained sections
- [x] All 7 works referenced where applicable
- [x] Prompt #5 contradictions cite specific findings: P1 PR-AUC=0.187, P3 pAE<10 filter, P7 1-2%
- [x] Prompt #6 traces 3 concepts (ipTM, data scarcity, success rates) through corpus
- [x] Prompt #9 synthesis is 380-420 words, references all 7 works
- [x] Prompt #11 knowledge map has central claim, 4 pillars, 2 contested zones, 3 must-reads
- [x] Prompt #12 produces exactly 3 deliverables

### Specs File Verification
- [x] References both papers with DOIs
- [x] References all 5 project websites (P3-P7)
- [x] All 4 deliverables listed with file paths and status
- [x] Implementation sequence matches execution order
- [x] File manifest matches actual created files

### Cross-Check: Contradictions
- [x] Prompt #5 Contradiction 1: P3 pAE filter vs P1 pAE discrimination failure --- P1 Fig. 2B cited
- [x] Prompt #5 Contradiction 2: P3 3.3% pass rate vs P7 1-2% experimental rate --- both values sourced
- [x] Prompt #5 Contradiction 3: P2 fuzzy affinities vs P6 quantitative scores --- 48% figure from P2
- [x] Prompt #5 Contradiction 4: P5 training bias vs P1 train/test independence --- P1 Fig. 2C cited
