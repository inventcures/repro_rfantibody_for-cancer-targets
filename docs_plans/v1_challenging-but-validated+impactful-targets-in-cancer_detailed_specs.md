# Challenging but Validated & Impactful Targets for De Novo Antibody Design in Cancer

> **Version**: 1.0 | **Date**: 2026-02-11 | **Author**: tp53 (Ashish)
>
> Deep research on cancer targets NOT in the original RFAntibody paper (Watson et al., Nature 2025),
> selected for high societal impact and unmet therapeutic need. Compiled from OpenTargets API,
> RCSB PDB, PubMed, and clinical trial registries (2023-2026).

---

## Table of Contents

1. [Executive Summary](#1-executive-summary)
2. [MPNST (Malignant Peripheral Nerve Sheath Tumor)](#2-mpnst)
3. [Pediatric Gliomas (DIPG/DMG)](#3-pediatric-gliomas)
4. [Neuroblastoma](#4-neuroblastoma)
5. [Glioblastoma (GBM)](#5-glioblastoma)
6. [Pancreatic Ductal Adenocarcinoma (PDAC)](#6-pancreatic-cancer)
7. [Cross-Indication Analysis](#7-cross-indication-analysis)
8. [Structural Readiness for RFAntibody Design](#8-structural-readiness)
9. [Recommended Campaign Configurations](#9-recommended-campaigns)
10. [Sources](#10-sources)

---

## 1. Executive Summary

### Why These Indications

| Indication | US Incidence | 5-Year Survival | Approved Ab Therapies | Unmet Need |
|---|---|---|---|---|
| **MPNST** | ~1,500/yr | ~35% (localized); ~15% (metastatic) | None | No systemic therapy beyond doxorubicin |
| **DIPG/DMG** | ~300/yr (pediatric) | <1% at 5 years | None | Median OS ~11 months; universally fatal |
| **Neuroblastoma** | ~800/yr (pediatric) | ~50% (high-risk) | Dinutuximab, Naxitamab (anti-GD2) | Relapsed/refractory remains lethal |
| **GBM** | ~13,000/yr | ~6.9% | Bevacizumab (no OS benefit) | Median OS ~15 months; no curative therapy |
| **PDAC** | ~64,000/yr | ~12% | None approved for PDAC specifically | 3rd leading cancer death; near-universal lethality |

### Key Cross-Indication Findings

**B7-H3 (CD276)** is the single most promising pan-cancer antibody target across these indications:
- Overexpressed in 4/5 indications (MPNST, pediatric gliomas, neuroblastoma, GBM)
- FDA Breakthrough Therapy Designation for DIPG (April 2025)
- >90% tumor expression with limited normal tissue expression
- Multiple modalities in clinical trials (CAR-T, ADC, radioimmunotherapy)

**GD2** is the most clinically validated target (3 FDA-approved antibodies for neuroblastoma) and shows dramatic responses in DIPG (Phase I, Nature 2024). However, it is a glycolipid, limiting protein-based de novo design approaches.

**Structural data is the bottleneck** for computational antibody design in these indications. Most high-priority targets lack publicly deposited antibody-antigen co-crystal structures:

| Structural Readiness | Targets |
|---|---|
| **Excellent** (multiple Ab-Ag complexes) | EGFR/EGFRvIII, HER2, PD-L1, VEGF-A |
| **Good** (1+ Ab-Ag complex) | MSLN, CD47, EphA2, GPC2 |
| **Limited** (target structure only) | B7-H3, IL-13Ra2, Nectin-4, TROP-2, ALK |
| **Poor** (no structure / glycolipid) | CLDN18.2, GD2, DLL3, MUC16, GPC1 |

### Top 5 Targets for RFAntibody De Novo Design (by structural readiness + clinical impact)

1. **GPC2** (neuroblastoma) -- 6WJL: D3 Fab complex at 3.3A; tumor-restricted expression
2. **MSLN** (pancreatic cancer) -- 4F3F: MORAb-009 Fab complex at 2.6A; 85-89% PDAC expression
3. **EGFR/EGFRvIII** (GBM, MPNST) -- 1YY9/8UKV/1I8I: multiple complexes at 1.8-2.8A
4. **EphA2** (GBM) -- 3SKJ: 1C1 Fab complex at 2.5A; tumor + vasculature expression
5. **CD47** (GBM) -- 5IWL: magrolimab complex; redirects TAMs (dominant GBM immune cell)

---

## 2. MPNST

### Disease Context

Malignant Peripheral Nerve Sheath Tumors are aggressive sarcomas arising from Schwann cells. ~50% occur in patients with Neurofibromatosis Type 1 (NF1). The core driver biology is loss of NF1 (constitutive RAS activation) and PRC2 components (SUZ12/EED). Most genetic drivers are intracellular, making surface target identification challenging. No systemic therapy has shown meaningful benefit beyond doxorubicin.

**OpenTargets Disease ID**: EFO_0000760 | **Total associated targets**: 1,104

### Target Ranking

| Rank | Target | UniProt | Evidence | Best PDB | Ab-Ag Complex? | Design Readiness |
|---|---|---|---|---|---|---|
| 1 | **B7-H3 (CD276)** | Q5ZPR3 | Strong preclinical (pan-sarcoma) | 5CMA (2.5A, Fab only) | No | Limited |
| 2 | **EGFR** | P00533 | Strong preclinical + correlative clinical | 1YY9 (2.6A) | Yes (cetuximab) | Excellent |
| 3 | **ErbB3/HER3** | P21860 | Strong preclinical (functional screens) | -- | No | Poor |
| 4 | **GD2** | glycolipid | Neural crest rationale; variable expression | 3VFG (1.65A, Fab only) | No | Limited |
| 5 | **PDGFRA** | P16234 | OT score 0.295; cooperates with NF1 loss | 5K5X (kinase only) | No | Poor |
| 6 | **HER2** | P04626 | Protein detected but genomically lost | 1N8Z (2.5A) | Yes (trastuzumab) | Excellent (but target invalid) |
| 7 | **MET** | P08581 | **NOT RECOMMENDED** -- literature explicitly excludes MPNST | 4K3J (2.8A) | Yes | N/A |

### Top Targets -- Detailed

#### B7-H3 (CD276) -- Rank #1

B7-H3/CD276 RNA overexpression observed in 58.3% (14/24) of sarcoma patients -- the highest frequency among immune checkpoint targets in sarcomas. An Fc-optimized anti-B7-H3 antibody induced NK cell reactivity against sarcoma cells in preclinical studies. >50% of preclinical studies have used B7-H3-targeting ADCs effective in sarcoma in vivo models.

**Known antibodies**: Enoblituzumab (Fc-engineered, Phase I), Omburtamab (131I-8H9, radioimmunotherapy), HS-20093 (ADC, ARTEMIS-002 Phase 2: ORR 25% in non-osteosarcoma sarcomas), DS-7300 (ADC, Phase 1/1b).

**PDB structures**: 5CMA (2.5A, ch8H9 Fab only), 4I0K (2.97A, murine B7-H3 ECD only). **No Ab-Ag complex structure available.**

**MPNST-specific rationale**: MPNSTs are immunologically "cold" -- B7-H3 ADC payload mechanism bypasses this. B7-H3 is immunosuppressive (inhibits T/NK cells), so targeting it could both kill tumor cells and relieve immune suppression.

**Key gap**: No MPNST-specific B7-H3 expression studies published yet.

#### EGFR -- Rank #2

EGFR amplification found in ~28% of MPNSTs. High EGFR expression correlates with poor disease-free and overall survival. EGFR knockdown or gefitinib inhibition causes dose-dependent inhibition of proliferation, migration, and invasion via PI3K/AKT and MAPK attenuation.

**Known antibodies**: Cetuximab (1YY9), Necitumumab (3B2U), Panitumumab (5SX4). All FDA-approved for other indications; none tested in MPNST trials.

**PDB structures**: 1YY9 (2.6A, cetuximab complex), 3B2U (2.58A, necitumumab complex), 3C09 (matuzumab complex), 4UV7 (GC1118 complex). **Excellent structural coverage.**

**Challenge**: NF1 loss drives constitutive RAS activation downstream of EGFR, creating a partial bypass. Combination with MEK inhibitor may be required.

#### ErbB3/HER3 -- Rank #3

The most robustly validated surface RTK in MPNST by functional screens. ErbB3 knockdown inhibits MPNST proliferation and survival. Pharmacologic and genome-scale shRNA screens consistently identify ErbB kinases as required for MPNST cell proliferation. NRG1-ErbB3 axis is a critical proliferative and survival signal specifically in Schwann cell-derived tumors.

**Known antibodies**: Patritumab (fully human), Seribantumab (Phase II in NRG1 fusion-positive cancers), Elgemtumab (locks inactive conformation).

**Key advantage**: ErbB3 is kinase-dead -- it MUST heterodimerize for signaling, making antibody-based blockade particularly well-suited vs. small molecules.

**Key gap**: No ErbB3 Ab-Ag complex structure in PDB.

### Strategic Recommendations for MPNST

1. **B7-H3 ADC**: Most immediately actionable. HS-20093 and DS-7300 already in sarcoma trials. Pursue MPNST-specific cohort.
2. **EGFR + MEK inhibitor**: Rational combination addressing both upstream signal and NF1-driven downstream bypass.
3. **Anti-ErbB3 for NRG1-high MPNST**: Screen MPNSTs for NRG1/ErbB3 expression; enroll in existing anti-HER3 basket trials (seribantumab).
4. **EGFR x B7-H3 bispecific**: Addresses multi-receptor biology while engaging immune mechanisms.

---

## 3. Pediatric Gliomas (DIPG/DMG)

### Disease Context

Diffuse Intrinsic Pontine Glioma (DIPG) / Diffuse Midline Glioma (DMG) is the deadliest pediatric brain cancer. ~80% harbor the H3K27M driver mutation. Median overall survival is ~11 months from diagnosis. No approved therapies exist. The brainstem location makes surgical resection impossible. The blood-brain barrier (BBB) blocks systemic antibody delivery (~0.1% penetration for IgG).

**OpenTargets Disease IDs**: EFO_1000026 (DIPG, 449 targets), EFO_0020983 (DMG, 162 targets), MONDO_1010030 (pediatric HGG, 58 targets)

### Target Ranking

| Rank | Target | Type | Evidence | Clinical Status | Design Readiness |
|---|---|---|---|---|---|
| 1 | **GD2** | Glycolipid | Phase I (Nature 2024): CR >30mo | RMAT designation | Limited (glycolipid) |
| 2 | **B7-H3 (CD276)** | Protein | Phase I complete; **FDA Breakthrough (2025)** | Phase 2 pivotal planned | Moderate (9LME: nanobody complex, 2.4A) |
| 3 | **H3K27M** | Intracellular | Phase I (vaccine): immunogenic | Multiple vaccine trials | N/A (intracellular) |
| 4 | **IL-13Ra2** | Protein | Phase I (adult GBM); pediatric ongoing | 65 patients treated (City of Hope) | Moderate (3LB6: ligand complex, 3.05A) |
| 5 | **PDGFRA** | Protein | Genetic (~15% HGG); avapritinib BBB-penetrant | Early clinical | Moderate (kinase structures only) |
| 6 | **EGFRvIII** | Protein | Phase I (adult GBM); heterogeneous in DIPG | CAR-T + EGFR806 trials | Good (1I8I: 1.8A; 8UKV: nanobody) |

### Top Targets -- Detailed

#### GD2 -- Rank #1 (Clinical Readiness)

**Phase I results (Stanford, Nature Nov 2024)**: 13 patients. IV + ICV GD2-CAR T cells. Four patients showed major volumetric tumor reductions (52%, 54%, 91%, 100%). One patient achieved **complete response ongoing >30 months**. Nine patients showed neurological benefit.

GD2 is **uniformly and highly expressed** on H3K27M-mutant DMG/DIPG cells. H3K27M mutation drives GD2 overexpression through upregulation of ganglioside biosynthesis genes.

**PDB**: 3VFG (1.65A, 3F8 Fab only). No GD2-antibody co-crystal exists (glycolipid).

**Challenge for de novo design**: GD2 is a ganglioside, not a protein. Cannot use standard protein-directed antibody discovery. Severe neuropathic pain from GD2 expression on peripheral nerves. Requires ICV delivery for sustained CNS responses.

#### B7-H3 (CD276) -- Rank #2 (Most Promising Protein Target)

**Phase I results (Seattle Children's, Nature Medicine Jan 2025)**: 21 DIPG patients. ICV B7-H3 CAR-T cells. Median survival from CAR-T: **10.7 months**; from diagnosis: **19.8 months**. Three patients alive at 44, 45, and 52 months from diagnosis. **FDA Breakthrough Therapy Designation** granted April 2025.

B7-H3 is highly and uniformly overexpressed on DIPG tumor cells AND tumor vasculature, while having limited expression on normal brain.

**PDB**: 5CMA (2.5A, ch8H9 Fab only), 4I0K (2.97A, murine B7-H3 ECD), **9LME (2.4A, B7-H3 + nanobody complex, deposited Jan 2025)**.

**Quad-targeting**: BrainChild-04 (NCT05768880) is the world's first **quad-targeting CAR** trial: B7-H3 + HER2 + EGFR806 + IL-13Ra2 zetakine, open for DMG/DIPG patients ages 1-26.

#### IL-13Ra2 -- Rank #4

Overexpressed in 75-80% of pediatric HGG. A GPI-anchored decoy receptor with minimal normal brain expression. Phase I (City of Hope): 65 patients, 50% achieved stable disease or better with locoregional delivery.

**PDB**: 3LB6 (3.05A, IL-13 ligand complex). **No antibody complex structure available.**

**Key for de novo design**: The IL-13 binding interface (3LB6) defines the targetable surface. Antibodies must not cross-react with IL13RA1/IL4RA signaling complex.

### Cross-Cutting Challenge: Blood-Brain Barrier

All successful clinical approaches in DIPG required locoregional delivery:
- **ICV**: B7-H3 CAR-T, GD2 CAR-T (after initial IV)
- **CED**: 8H9 radioimmunotherapy, D2C7-IT immunotoxin
- **Combined IV + intracranial**: GD2 CAR-T (Stanford protocol)

Systemic antibodies achieve <0.1% brain penetration. Smaller formats (VHH ~15 kDa) may have improved BBB penetration vs. full IgG (~150 kDa), making RFAntibody-designed VHH nanobodies potentially advantageous for CNS tumors.

---

## 4. Neuroblastoma

### Disease Context

Neuroblastoma is the most common extracranial solid tumor of childhood, arising from neural crest cells. High-risk neuroblastoma (MYCN-amplified, >18 months at diagnosis, metastatic) has ~50% survival despite intensive multimodal therapy. Anti-GD2 immunotherapy (dinutuximab) is standard of care but relapsed/refractory disease remains lethal.

### Target Ranking

| Rank | Target | UniProt | Evidence | Best PDB | Ab-Ag Complex? | Design Readiness |
|---|---|---|---|---|---|---|
| 1 | **GPC2** | Q8N158 | Phase I CAR-T (NCT05650749) | 6WJL (3.3A) | **Yes (D3 Fab)** | **HIGH** |
| 2 | **B7-H3** | Q5ZPR3 | Phase 1/1b (ADC, CAR-T) | 5CMA / 9LME (2.4A) | Nanobody (9LME) | Moderate |
| 3 | **DLK1** | P80370 | Phase 1 FIH (ADCT-701, CBA-1205) | 9D20 (2.67A, ACVR2B complex) | No | Moderate |
| 4 | **ALK** | Q9UM73 | OT score 0.81; lorlatinib Phase I | 7LRZ (1.91A, GRD) | No | Moderate |
| 5 | **L1CAM** | P32004 | Phase 1 CAR-T; Ab612 entering FIH | 8AFO (1.99A, FNIII domains) | No | Moderate |
| 6 | **PHOX2B** | Q99453 | OT score 0.65; intracellular TF | None (pMHC) | No | Poor |
| 7 | **NCAM1** | P13591 | Phase 2 (lorvotuzumab): limited activity | 5AEA | No | Deprioritized |
| 8 | **GD2** | glycolipid | **FDA-approved** (dinutuximab, naxitamab) | 3VFG (1.65A) | No | Established; glycolipid |

### Top Targets -- Detailed

#### GPC2 (Glypican-2) -- Rank #1 for De Novo Design

**The strongest candidate for computational antibody design among all neuroblastoma targets.** GPC2 is a heparan sulfate proteoglycan highly expressed on most neuroblastomas but essentially absent from normal childhood tissues -- the best differential expression profile among novel targets.

**Structural data**: 6WJL (3.3A, GPC2 core protein + D3 Fab complex). D3 recognizes a discontinuous, conformation-dependent epitope spanning helices 2, 12, and 13 of the GPC2 core protein. The D3 Fab buries 1,069 A^2 with multiple salt bridges and hydrophobic contacts. This is a tumor-specific epitope not shared with other glypicans.

**Clinical status**: First-in-human CAR-T trial ongoing (NCT05650749). D3-GPC2-PBD ADC effective in preclinical models. Combined GPC2 + B7-H3 dual-targeting CAR-T being developed.

**RFAntibody design opportunity**: The 6WJL complex structure provides the epitope definition needed for de novo VHH/scFv design targeting the D3 epitope or exploring alternative epitopes on the GPC2 surface.

#### DLK1 (Delta-Like Canonical Notch Ligand 1) -- Rank #3

Identified through unbiased proteogenomic surfaceome analysis at CHOP (Cancer Cell, Oct 2024). DLK1 is robustly expressed on the neuroblastoma cell surface and largely restricted to adrenal medulla and pituitary gland in normal tissues.

**Clinical status**: Two agents in Phase 1: ADCT-701 (ADC, NCT06041516) and CBA-1205 (humanized anti-DLK1 mAb with enhanced ADCC).

**PDB**: 9D20 (2.67A, DLK1 + ACVR2B complex -- not an antibody complex). No antibody-DLK1 co-crystal available.

#### ALK -- Rank #4

Strongest genetic evidence (OT score 0.81). Mutated/amplified in ~14% of high-risk neuroblastoma. Lorlatinib (small molecule) showed 30% response rate (single agent) and 63% with chemo combo (Phase I, Nature Medicine 2023). Now in Phase 3 COG trial.

**PDB**: 7LRZ (1.91A, ALK glycine-rich domain of ECD). No antibody-ALK ECD complex exists. The ECD is enormous (~1000 aa) with MAM, LDLa, and glycine-rich domains.

**Antibody opportunity**: ALK inhibitor resistance mutations may not affect antibody epitopes, providing a differentiated therapeutic approach. ALK.CAR-T + ALK inhibitor synergy demonstrated preclinically.

---

## 5. Glioblastoma (GBM)

### Disease Context

GBM is the most aggressive primary brain tumor in adults. Median overall survival ~15 months with standard Stupp protocol (surgery + TMZ + RT). Bevacizumab (anti-VEGF-A) is FDA-approved for recurrent GBM but provides NO overall survival benefit. No antibody therapy has demonstrated OS improvement in GBM. The immunosuppressive tumor microenvironment (TAMs constitute up to 50% of tumor mass), BBB, and extreme intratumoral heterogeneity are the core challenges.

**OpenTargets Disease ID**: EFO_0000519 | **Total associated targets**: 9,906

### Target Ranking

| Rank | Target | UniProt | Evidence | Best PDB | Ab-Ag Complex? | Design Readiness |
|---|---|---|---|---|---|---|
| 1 | **EGFR/EGFRvIII** | P00533 | OT #1 (0.683); ~57% amplified; Phase I-III | 1YY9/8UKV/1I8I (1.8-2.8A) | Yes (multiple) | **Excellent** |
| 2 | **IL-13Ra2** | Q14627 | Phase I CAR-T (dramatic CR published) | 3LB6 (1.85A, ligand only) | No | Moderate |
| 3 | **B7-H3** | Q5ZPR3 | Phase I (CAR-T, radioimmunotherapy) | 5CMA/9LME | Nanobody (9LME) | Moderate |
| 4 | **HER2** | P04626 | Phase I CAR-T; ~80% GBM expression | 1N8Z (2.5A) | Yes (trastuzumab) | **Excellent** |
| 5 | **CD47** | Q08722 | Preclinical; Phase I/II (other tumors) | 5IWL | Yes (magrolimab) | **Good** |
| 6 | **EphA2** | P29317 | Phase I ADC (halted); preclinical CAR-T | 3SKJ (2.5A) | Yes (1C1 mAb) | **Good** |
| 7 | **VEGF-A** | P15692 | OT #14 (0.550); **FDA-approved** (bevacizumab) | 1BJ1 | Yes | Excellent (but no OS benefit) |
| 8 | **PD-L1** | Q9NZQ7 | Phase II/III (all negative in GBM) | 5X8L (2.9A) | Yes (atezolizumab) | Excellent (but target failed) |
| 9 | **DLL3** | Q9NYJ7 | Preclinical GBM; FDA-approved for SCLC | None (homology from DLL1) | No | Poor |
| 10 | **GD2** | glycolipid | Phase I CAR-T (DMG, Nature 2022) | 3VFG (Fab only) | No | Limited |
| 11 | **TROP-2** | P09758 | Phase 0 GBM: PFS 2 months only | 7PEE (ECD only) | No | Moderate |

### Top Targets -- Detailed

#### EGFR/EGFRvIII -- Rank #1

EGFR amplification in ~57% of primary GBM. EGFRvIII (exons 2-7 deletion) in ~25-30% of newly diagnosed GBM, completely absent from normal tissues.

**Most promising clinical approach**: Bivalent CART-EGFR-IL13Ra2 (NCT05168423). ASCO 2025 update: 7 of 10 patients surviving >8 months. This dual-targeting approach addresses the heterogeneity problem.

**PDB structures**:
- 1YY9 (2.8A): EGFR ECD + cetuximab Fab. Domain III epitope.
- 8UKV: Nanobody 34E5 + EGFRvIII ECD. EGFRvIII-specific.
- 8UKX: EGFRvIII ECD alone at pH 7.0. Compact conformation.
- 1I8I (1.8A): dsFv MR1 + EGFRvIII junction neoepitope. **Highest resolution EGFRvIII-specific structure.**
- mAb 806 epitope (residues 287-302): cryptic epitope exposed only on EGFRvIII and amplified/overexpressed EGFR, NOT on normal EGFR. Provides tumor selectivity.

**RFAntibody design opportunity**: Design VHH nanobodies targeting the EGFRvIII-specific conformation (using 8UKV/8UKX as templates) for improved BBB penetration over full IgG.

#### CD47 -- Rank #5 (Mechanistically Distinct)

CD47 overexpression on GBM cells inhibits phagocytosis by tumor-associated macrophages (TAMs), which constitute up to 50% of the GBM tumor mass. Anti-CD47 redirects the **dominant immune cell population** against the tumor, rather than relying on T cells.

**PDB**: 5IWL (CD47 + magrolimab diabody complex), 2JJS (CD47 + SIRPa).

Magrolimab binds N-terminal pyroglutamate of CD47 and BC/FG loops, competitively blocking SIRPa. Development halted by Gilead for hematologic malignancies but VT-1021 (cyclic peptide) in Phase III for recurrent GBM (NCT03970447).

**Why this is important**: GBM is "cold" with minimal T-cell infiltration. T-cell-dependent mechanisms (BiTEs, checkpoint inhibitors, CAR-T) face an inherently hostile microenvironment. Macrophage-redirecting anti-CD47 is mechanistically orthogonal.

#### EphA2 -- Rank #6

Overexpressed in GBM, including tumor neovasculature. The 1C1 antibody (MEDI-547) binds the ligand-binding domain through HCDR3-mediated ligand mimicry, inducing receptor agonism, internalization, and degradation. Phase I halted due to hemorrhagic toxicity (EphA2 on normal vasculature).

**PDB**: 3SKJ (2.5A, EphA2 + 1C1 Fab), 2X10 (2.4A, full ECD), 2X11 (ECD + ephrinA5).

**RFAntibody opportunity**: Design antibodies targeting non-overlapping epitopes (away from ligand-binding domain) that may avoid the vascular toxicity of 1C1.

---

## 6. Pancreatic Cancer (PDAC)

### Disease Context

Pancreatic ductal adenocarcinoma is the 3rd leading cause of cancer death. Five-year survival ~12%. The hallmark challenges for antibody therapy are:

1. **Dense desmoplastic stroma**: Up to 80% of tumor mass is extracellular matrix generated by cancer-associated fibroblasts (CAFs). IgG (~150 kDa) has severely limited diffusion.
2. **Poor vascularization**: Hypovascular tumors with elevated interstitial fluid pressure (2-3x normal).
3. **Immunosuppressive TME**: CXCL12 blocks CD8+ T-cell trafficking. High TGF-beta and IL-10. Checkpoint inhibitors have failed as single agents.
4. **Antigen shedding**: MSLN and MUC16/CA-125 are shed, creating antigen sinks.

These challenges make PDAC uniquely resistant to antibody-based therapies. **ADCs outperform naked antibodies**, and **smaller formats (VHH, scFv) could improve penetration.**

### Target Ranking

| Rank | Target | UniProt | Expression in PDAC | Clinical Status | Design Readiness |
|---|---|---|---|---|---|
| 1 | **CLDN18.2** | P56856 | ~28-33% | Phase II (zolbetuximab); ADC IBI343 FDA Fast Track | **Poor** (no structure) |
| 2 | **CEACAM5/6** | P06731 | Broadly expressed | EBC-129 ADC: 20% ORR, FDA Fast Track | Moderate (8BW0: 3.11A) |
| 3 | **MSLN** | Q13421 | 85-89% | Phase I/Ib (anetumab ravtansine) | **Good** (4F3F: 2.6A) |
| 4 | **EGFR** | P00533 | ~5-10% eligible (KRAS-WT) | Phase III (nimotuzumab: OS 10.9 vs 8.5 mo) | Excellent |
| 5 | **Nectin-4** | Q96NY8 | ~71% positivity | Preclinical (enfortumab vedotin) | Moderate (4FRW + 9MW2821) |
| 6 | **HER2** | P04626 | 2-3% (IHC 3+) | T-DXd: **4% ORR** in PDAC (failed) | Excellent (but target failed) |
| 7 | **TROP-2** | P09758 | Broadly expressed | Phase I/II basket (limited PDAC activity) | Moderate (7PEE) |
| 8 | **GPC1** | P35052 | Overexpressed | Phase Ib planned (177Lu-Miltuximab) | Poor |
| 9 | **MUC16** | Q8WXI7 | Co-expressed with MSLN | Phase I discontinued | Poor |
| 10 | **DLL3** | Q9NYJ7 | PNEC only | FDA-approved for SCLC; not PDAC | Poor |

### Top Targets -- Detailed

#### CLDN18.2 -- Rank #1 (Most Clinically Advanced New Target)

Claudin-18.2 is a tight junction protein with gastric lineage-restricted expression. ~28-33% of PDAC tumors express it. Zolbetuximab (anti-CLDN18.2 mAb) in Phase II randomized trial (NCT03816163, n=396). IBI343 (anti-CLDN18.2 ADC) received FDA Fast Track designation for advanced PDAC.

**Critical structural gap**: No CLDN18.2 crystal structure exists in PDB. Claudin-4 structures (7KP4, 8U4V) serve as homology templates. Zolbetuximab epitope (ECL1 loop) described but not structurally resolved.

**RFAntibody challenge**: Would require homology modeling from claudin-4, introducing uncertainty. Epitope resides in the short extracellular loops of a 4-pass transmembrane protein.

#### MSLN (Mesothelin) -- Rank #3 (Best for Computational Design)

Expressed in 85-89% of PDAC. The highest expression rate among all PDAC targets. Multiple Phase I/Ib trials completed (anetumab ravtansine ADC).

**PDB structures**:
- 8CX3: Full-length mesothelin structure
- 7UED: Full-length MSLN + MORAb-009 (amatuximab) Fab complex
- **4F3F (2.6A)**: MSLN fragment + MORAb-009 Fab. **The gold-standard structure for PDAC antibody design.**
- 7U8C: MSLN C-terminal + 15B6 Fab (different epitope)

**RFAntibody design opportunity**: Multiple epitopes defined structurally. Design VHH nanobodies (~15 kDa) for improved stromal penetration, targeting either the N-terminal or C-terminal epitope.

#### CEACAM5/6 -- Rank #2 (Best Clinical ORR)

EBC-129 (anti-CEACAM5/6 ADC) showed **20% ORR and 71.4% DCR** in heavily pretreated 2L+ PDAC at ASCO 2025. FDA Fast Track. Phase II launching 2026.

**PDB**: 8BW0 (3.11A, CEACAM5 A3-B3 domains + tusamitamab Fab, cryo-EM). Epitope-paratope interaction paper published (Nature Communications 2024).

#### Unique Target: GPC1 (Glypican-1)

GPC1 is uniquely interesting because it is expressed on **both tumor cells AND cancer-associated fibroblasts (CAFs)**, enabling dual targeting of tumor and stroma. Phase Ib therapeutic dose escalation with 177Lu-DOTA-Miltuximab (radioimmunotherapy) planned. The radioimmunotherapy approach may be less dependent on deep tissue penetration.

---

## 7. Cross-Indication Analysis

### Targets Appearing Across Multiple Indications

| Target | MPNST | DIPG/DMG | Neuroblastoma | GBM | PDAC | Total |
|---|---|---|---|---|---|---|
| **B7-H3 (CD276)** | #1 | #2 | #2 | #3 | -- | **4/5** |
| **GD2** | #4 | #1 | Established | #10 | -- | **4/5** |
| **EGFR/EGFRvIII** | #2 | #6 | -- | #1 | #4 | **4/5** |
| **HER2** | #6 | -- | -- | #4 | #6 | **3/5** |
| **PDGFRA** | #5 | #5 | -- | -- | -- | **2/5** |
| **IL-13Ra2** | -- | #4 | -- | #2 | -- | **2/5** |
| **DLL3** | -- | -- | -- | #9 | #10 | **2/5** |
| **TROP-2** | -- | -- | -- | #11 | #7 | **2/5** |

### Shared Structural Data Opportunities

B7-H3 structural data benefits 4 indications simultaneously. The recently deposited 9LME (B7-H3 + nanobody, 2.4A, Jan 2025) is the first publicly available B7-H3 antibody complex and enables de novo design for MPNST, DIPG, neuroblastoma, and GBM.

EGFR structural data (1YY9, 8UKV, 1I8I) is directly applicable to MPNST (wild-type EGFR) and GBM (EGFRvIII-specific), with different epitope strategies for each indication.

### Indication-Specific Unique Targets

| Indication | Unique High-Value Targets |
|---|---|
| MPNST | ErbB3/HER3 (functionally validated, kinase-dead) |
| DIPG/DMG | H3K27M neoantigen (100% tumor-specific, intracellular) |
| Neuroblastoma | GPC2 (best differential expression + structural data), DLK1 (newest) |
| GBM | CD47 (redirects TAMs, mechanistically orthogonal to T-cell approaches) |
| PDAC | MSLN (85-89% expression + best structural data), CLDN18.2 (most clinically advanced), GPC1 (targets tumor + stroma) |

---

## 8. Structural Readiness for RFAntibody Design

The RFAntibody pipeline requires:
1. Target PDB structure (ideally extracellular domain)
2. Defined epitope residues (ideally from Ab-Ag co-crystal)
3. Hotspot selection for CDR loop placement

### Tier 1: Excellent (Multiple Ab-Ag Complexes, <3A Resolution)

| Target | Best Complex PDB | Resolution | Epitopes Defined | Indication Priority |
|---|---|---|---|---|
| **EGFR** | 1YY9 (cetuximab) | 2.6A | Domain III (cetuximab, necitumumab) | GBM, MPNST |
| **EGFRvIII** | 1I8I (MR1 dsFv) | 1.8A | Junction neoepitope (LEEKKGNYVVTDH) | GBM |
| **EGFRvIII** | 8UKV (nanobody 34E5) | -- | EGFRvIII-specific conformation | GBM |
| **HER2** | 1N8Z (trastuzumab) | 2.5A | Domain IV (trastuzumab), Domain II (pertuzumab) | GBM |
| **PD-L1** | 5X8L (atezolizumab) | 2.9A | IgV domain front face | GBM (failed) |
| **VEGF-A** | 1BJ1 (bevacizumab) | -- | Receptor binding site | GBM (no OS benefit) |

### Tier 2: Good (1 Ab-Ag Complex Available)

| Target | Best Complex PDB | Resolution | Notes | Indication Priority |
|---|---|---|---|---|
| **GPC2** | 6WJL (D3 Fab) | 3.3A | Discontinuous epitope on helices 2/12/13 | Neuroblastoma |
| **MSLN** | 4F3F (MORAb-009 Fab) | 2.6A | N-terminal epitope; also 7U8C (C-term) | PDAC |
| **CD47** | 5IWL (magrolimab) | -- | Competitive with SIRPa | GBM |
| **EphA2** | 3SKJ (1C1 Fab) | 2.5A | Ligand mimicry by HCDR3 | GBM |
| **B7-H3** | 9LME (nanobody) | 2.4A | First Ab complex (Jan 2025) | MPNST, DIPG, NB, GBM |
| **CEACAM5** | 8BW0 (tusamitamab) | 3.11A | A3-B3 domain epitope | PDAC |

### Tier 3: Limited (Target Structure Only, No Ab Complex)

| Target | Best Target PDB | Resolution | Notes | Indication Priority |
|---|---|---|---|---|
| **IL-13Ra2** | 3LB6 (+ IL-13 ligand) | 3.05A | Ligand defines binding surface | GBM, DIPG |
| **Nectin-4** | 4FRW (D1-D2 ECD) | -- | Cryo-EM complex (9MW2821, 3.26A) | PDAC |
| **TROP-2** | 7PEE (ECD) | -- | Sacituzumab epitope mapped but no complex | PDAC |
| **ALK** | 7LRZ (GRD domain) | 1.91A | Enormous ECD; no Ab complex | Neuroblastoma |
| **DLK1** | 9D20 (+ ACVR2B) | 2.67A | Protein topology defined; no Ab complex | Neuroblastoma |
| **L1CAM** | 8AFO (FNIII domains) | 1.99A | Multiple Ig+FNIII domains; no Ab complex | Neuroblastoma |

### Tier 4: Insufficient (No Usable Structure)

| Target | Issue | Indication Priority |
|---|---|---|
| **CLDN18.2** | No structure; homology from claudin-4 only | PDAC |
| **GD2** | Glycolipid, not protein | DIPG, NB, MPNST, GBM |
| **ErbB3/HER3** | No ECD structure or Ab complex | MPNST |
| **DLL3** | No structure; homology from DLL1 | GBM, PDAC |
| **PDGFRA** | No ECD structure (kinase domain only) | MPNST, DIPG |
| **GPC1** | No Ab complex; minimal structural data | PDAC |
| **MUC16** | ~22,000 aa mucin; minimal structural data | PDAC |

---

## 9. Recommended Campaign Configurations

For each recommendation: PDB ID for target preparation, suggested epitope region, antibody format (VHH preferred for BBB penetration and stromal access), and key design considerations.

### Campaign 1: GPC2 for Neuroblastoma

| Parameter | Value |
|---|---|
| **Target PDB** | 6WJL (GPC2 + D3 Fab complex) |
| **Epitope** | Helices 2, 12, 13 of GPC2 core protein (D3 epitope) |
| **Format** | VHH (nanobody) |
| **Rationale** | Best differential expression + structural data among novel neuroblastoma targets |
| **Design count** | 5,000-10,000 |
| **Validation** | Compare to D3 binding; assess cross-reactivity with other glypicans |

### Campaign 2: MSLN for Pancreatic Cancer

| Parameter | Value |
|---|---|
| **Target PDB** | 4F3F (MSLN + MORAb-009 Fab, 2.6A) + 7U8C (C-terminal epitope) |
| **Epitope** | N-terminal region (MORAb-009 site) OR C-terminal region (15B6 site) |
| **Format** | VHH (for stromal penetration; 15 kDa vs 150 kDa IgG) |
| **Rationale** | 85-89% PDAC expression; best structural data among PDAC targets; ADC-compatible |
| **Design count** | 5,000-10,000 per epitope |
| **Validation** | SPR binding; competition with amatuximab/15B6; internalization assay for ADC |

### Campaign 3: EGFRvIII for GBM

| Parameter | Value |
|---|---|
| **Target PDB** | 8UKV (EGFRvIII ECD + nanobody 34E5) or 8UKX (EGFRvIII ECD alone) |
| **Epitope** | EGFRvIII-specific conformation (tumor-selective) |
| **Format** | VHH (for BBB penetration) |
| **Rationale** | Completely tumor-specific; absent from normal tissue; VHH may cross BBB |
| **Design count** | 5,000-10,000 |
| **Validation** | Verify binding to EGFRvIII but NOT wild-type EGFR; BBB transcytosis assay |

### Campaign 4: B7-H3 for MPNST / DIPG / Neuroblastoma

| Parameter | Value |
|---|---|
| **Target PDB** | 9LME (B7-H3 + nanobody, 2.4A, Jan 2025) |
| **Epitope** | Nanobody epitope from 9LME or alternative face of IgV domain |
| **Format** | VHH (matches 9LME template; BBB access for DIPG) |
| **Rationale** | Pan-cancer target across 4/5 indications; FDA Breakthrough for DIPG |
| **Design count** | 5,000-10,000 |
| **Validation** | Confirm binding to 4Ig human isoform (not just 2Ig); immune checkpoint blockade assay |

### Campaign 5: CD47 for GBM

| Parameter | Value |
|---|---|
| **Target PDB** | 5IWL (CD47 + magrolimab diabody) |
| **Epitope** | SIRPa-competitive face (N-terminal pyroglutamate, BC/FG loops) |
| **Format** | VHH (BBB penetration) or scFv |
| **Rationale** | Mechanistically orthogonal to T-cell approaches; redirects TAMs (50% of GBM mass) |
| **Design count** | 5,000 |
| **Validation** | SIRPa competition assay; phagocytosis assay with GBM cell lines; RBC binding safety |

### Campaign 6: CEACAM5 for Pancreatic Cancer

| Parameter | Value |
|---|---|
| **Target PDB** | 8BW0 (CEACAM5 A3-B3 + tusamitamab Fab, 3.11A cryo-EM) |
| **Epitope** | A3-B3 domain interface (tusamitamab site) |
| **Format** | VHH (stromal penetration) |
| **Rationale** | Best clinical ORR in PDAC (EBC-129: 20%); broadly expressed |
| **Design count** | 5,000 |
| **Validation** | Competition with tusamitamab; ADC conjugation compatibility |

### Campaign 7: EphA2 for GBM

| Parameter | Value |
|---|---|
| **Target PDB** | 3SKJ (EphA2 + 1C1 Fab, 2.5A) + 2X10 (full ECD, 2.4A) |
| **Epitope** | Alternative face (AVOID ligand-binding domain to reduce vascular toxicity) |
| **Format** | VHH |
| **Rationale** | Dual tumor + vasculature expression; excellent structural data; MEDI-547 toxicity was epitope-specific |
| **Design count** | 5,000 |
| **Validation** | Verify NO agonist activity (unlike 1C1); normal vasculature binding safety |

---

## 10. Sources

### MPNST
- [OpenTargets Platform - MPNST (EFO_0000760)](https://platform.opentargets.org/disease/EFO_0000760)
- [PDGFRA, PDGFRB, EGFR activation in MPNST (PMC2802393)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2802393/)
- [PDGFRA cooperates with NF1 loss in MPNST (Nature/Oncogene)](https://www.nature.com/articles/onc2016269)
- [EGFR expression and significance in MPNST (J Neuro-Oncol 2009)](https://link.springer.com/article/10.1007/s11060-009-9862-z)
- [EGFR and ErbB2 in MPNST (PubMed 18650488)](https://pubmed.ncbi.nlm.nih.gov/18650488/)
- [Pan-ErbB inhibitor in MPNST (PMC3923431)](https://pmc.ncbi.nlm.nih.gov/articles/PMC3923431/)
- [ErbB3 inhibition in MPNST (PMC10477957, 2023)](https://pmc.ncbi.nlm.nih.gov/articles/PMC10477957/)
- [NRG1/ErbB signaling in neoplastic Schwann cells (PubMed 15897877)](https://pubmed.ncbi.nlm.nih.gov/15897877/)
- [B7-H3 Fc-optimized antibody in sarcoma (PubMed 36275693)](https://pubmed.ncbi.nlm.nih.gov/36275693/)
- [B7-H3 Inhibitors Clinical Trials (PMC10846634)](https://pmc.ncbi.nlm.nih.gov/articles/PMC10846634/)
- [B7-H3 ADC Phase 1/1b (Nature Medicine 2025)](https://www.nature.com/articles/s41591-025-03600-2)
- [GD2 in solid tumors (PMC7358363)](https://pmc.ncbi.nlm.nih.gov/articles/PMC7358363/)
- [GD2 14G2a structural basis (PMC4597138)](https://pmc.ncbi.nlm.nih.gov/articles/PMC4597138/)
- [MET unlikely to benefit MPNST (PMC4386816)](https://pmc.ncbi.nlm.nih.gov/articles/PMC4386816/)
- [MPNST pathogenesis (PMC9954030)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9954030/)

### Pediatric Gliomas (DIPG/DMG)
- [GD2 CAR-T Phase I (Nature Nov 2024)](https://www.nature.com/articles/s41586-024-08171-9)
- [GD2 CAR T Cells Show Promise (NCI 2025)](https://www.cancer.gov/news-events/cancer-currents-blog/2025/car-t-cell-therapy-gd2-diffuse-midline-gliomas)
- [B7-H3 CAR-T Phase 1 (Nature Medicine Jan 2025)](https://www.nature.com/articles/s41591-024-03451-3)
- [FDA Breakthrough for BCB-276 (April 2025)](https://www.businesswire.com/news/home/20250422418430/en/)
- [B7-H3 CAR-T First-in-Human (Cancer Discovery 2023)](https://aacrjournals.org/cancerdiscovery/article/13/1/114/712703/)
- [IL-13Ra2 CAR-T Phase 1 (Nature Medicine 2024)](https://www.nature.com/articles/s41591-024-02875-1)
- [PDGFRA Avapritinib (Cancer Cell 2025)](https://www.cell.com/cancer-cell/fulltext/S1535-6108(25)00070-4)
- [H3K27M Vaccine Adult DMG (Nature Medicine 2023)](https://www.nature.com/articles/s41591-023-02555-6)
- [H3K27M Neoepitope responses (Science Advances 2024)](https://www.science.org/doi/10.1126/sciadv.adi9091)
- [GD2.CART with IL-7R (JCO 2024)](https://ascopubs.org/doi/10.1200/JCO.23.02019)
- [EGFR-EDVsMIT Pediatric Phase 1 (Targeted Oncology 2024)](https://link.springer.com/article/10.1007/s11523-024-01051-2)
- [IL13RA2 GEMM for Pediatric HGG (Acta Neuropathologica 2025)](https://link.springer.com/article/10.1186/s40478-025-01991-4)
- [B7-H3 CAR-T + ONC206 (bioRxiv 2025)](https://www.biorxiv.org/content/10.1101/2025.09.01.673023v1.full)

### Neuroblastoma
- [Naxitamab Phase 2 (Nature Communications 2025)](https://www.nature.com/articles/s41467-025-56619-x)
- [GD2-CART01 Phase 1/2 (PubMed 40841488)](https://pubmed.ncbi.nlm.nih.gov/40841488/)
- [Emerging targeted therapies for neuroblastoma (Frontiers Oncology 2025)](https://www.frontiersin.org/journals/oncology/articles/10.3389/fonc.2025.1553511/full)
- [Targets and Antibody Formats for Neuroblastoma (JCO, PMC7255979)](https://pmc.ncbi.nlm.nih.gov/articles/PMC7255979/)
- [DLK1 proteogenomic surfaceome (Cancer Cell Oct 2024)](https://www.cell.com/cancer-cell/fulltext/S1535-6108(24)00366-0)
- [CBA-1205 Phase 1 FIH (PubMed 39832211)](https://pubmed.ncbi.nlm.nih.gov/39832211/)
- [GPC2 ADC reprograms immune milieu (PMC9723962)](https://pmc.ncbi.nlm.nih.gov/articles/PMC9723962/)
- [D3-GPC2 ADC conformational epitope (Cell Reports Medicine)](https://www.cell.com/cell-reports-medicine/fulltext/S2666-3791(21)00193-2)
- [D3-GPC2 CAR T preclinical 2025 (PubMed 41026583)](https://pubmed.ncbi.nlm.nih.gov/41026583/)
- [Lorlatinib Phase 1 (Nature Medicine 2023)](https://www.nature.com/articles/s41591-023-02297-5)
- [ALK.CAR-T + ALK inhibitor synergy](https://onlinelibrary.wiley.com/doi/pdfdirect/10.1002/ctm2.1732)
- [CE7-CAR T preclinical (Clin Cancer Res)](https://aacrjournals.org/clincancerres/article/23/2/466/275307/)
- [Ab612 nonclinical safety 2025](https://link.springer.com/article/10.1007/s00262-025-04142-9)
- [ADVL1522 Phase 2 lorvotuzumab (Cancer)](https://acsjournals.onlinelibrary.wiley.com/doi/abs/10.1002/cncr.33195)
- [2035 Clinical Research Vision for Neuroblastoma (Ped Blood Cancer 2025)](https://onlinelibrary.wiley.com/doi/10.1002/pbc.31660)

### GBM
- [OpenTargets Platform - GBM (EFO_0000519)](https://platform.opentargets.org/disease/EFO_0000519)
- [Targeted Glioma Therapy (PMC 2024)](https://pmc.ncbi.nlm.nih.gov/articles/PMC10820492/)
- [CAR-T for adult HGG (Nature 2024)](https://www.nature.com/articles/s41698-024-00753-0)
- [Immunotherapy for GBM (CMI 2024)](https://www.nature.com/articles/s41423-024-01226-x)
- [Overcoming GBM immunotherapy resistance (Frontiers 2025)](https://www.frontiersin.org/journals/pharmacology/articles/10.3389/fphar.2025.1584688/full)
- [GBM targeted therapies (EHO 2024)](https://ehoonline.biomedcentral.com/articles/10.1186/s40164-024-00512-8)
- [CAR-T advances in GBM (PMC 2025)](https://pmc.ncbi.nlm.nih.gov/articles/PMC12401244/)
- [Bispecific Antibodies in GBM (MedSci 2025)](https://www.medsci.org/v22p4250.pdf)
- [GD2-CAR T for DMG (Nature 2022)](https://www.nature.com/articles/s41586-022-04489-4)
- [DLL3 in IDH-mutant glioma (PMC 2020)](https://pmc.ncbi.nlm.nih.gov/articles/PMC7365589/)
- [Sacituzumab govitecan Phase 0 GBM (Nature Comms 2024)](https://www.nature.com/articles/s41467-024-50558-9)
- [BBB drug delivery for GBM (PMC 2024)](https://pmc.ncbi.nlm.nih.gov/articles/PMC11016125/)
- [CD47-SIRPa therapeutics (PMC 2022)](https://pmc.ncbi.nlm.nih.gov/articles/PMC9216458/)
- [Bevacizumab for newly diagnosed GBM (NEJM)](https://www.nejm.org/doi/full/10.1056/NEJMoa1308573)
- [Glioblastoma microenvironment (Genes & Dev 2024)](https://genesdev.cshlp.org/content/38/9-10/360.full)

### Pancreatic Cancer
- [Anti-CLDN18.2 ADC FDA Fast Track (CancerNetwork)](https://www.cancernetwork.com/view/anti-cldn18-2-adc-receives-fda-fast-track-designation-in-advanced-pdac)
- [EBC-129 in PDAC (Targeted Oncology)](https://www.targetedonc.com/view/early-results-show-promise-for-ebc-129-a-novel-adc-in-pancreatic-cancer)
- [EBC-129 FDA Fast Track (OncLive)](https://www.onclive.com/view/fda-grants-fast-track-status-to-ebc-129-for-pancreatic-ductal-adenocarcinoma)
- [CEACAM6-targeted ADC in PDAC (Nature Communications)](https://www.nature.com/articles/s41467-024-46167-1)
- [Zolbetuximab Phase II PDAC (Annals of Oncology)](https://www.annalsofoncology.org/article/S0923-7534(24)03114-4/fulltext)
- [CLDN18.2 ASCO 2024 update (J Hematol Oncol)](https://link.springer.com/article/10.1186/s13045-024-01595-w)
- [Anetumab ravtansine Phase I (JCO)](https://ascopubs.org/doi/10.1200/JCO.19.02085)
- [DESTINY-PanTumor02 results (JCO)](https://ascopubs.org/doi/10.1200/JCO.23.02005)
- [Nimotuzumab Phase III KRAS-WT PDAC (JCO)](https://ascopubs.org/doi/10.1200/JCO.22.02630)
- [MSLN overexpression in PDAC (PMC7053310)](https://pmc.ncbi.nlm.nih.gov/articles/PMC7053310/)
- [MSLN-MUC16 interaction (PMC3660778)](https://pmc.ncbi.nlm.nih.gov/articles/PMC3660778/)
- [Nectin-4 in PDAC (PMC11875723)](https://pmc.ncbi.nlm.nih.gov/articles/PMC11875723/)
- [GPC1 dual targeting (AACR)](https://aacrjournals.org/mct/article/20/12/2495/675151/)
- [PDAC TME barrier (Frontiers Immunology 2024)](https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2024.1287459/full)
- [CEACAM5 epitope-paratope insights (Nature Communications 2024)](https://www.nature.com/articles/s41467-024-53746-9)
- [Novel TROP2 amanitin ADC (AACR MCT)](https://aacrjournals.org/mct/article/24/4/485/754278/)
- [Zenocutuzumab approval for NRG1 PDAC (GI Oncology Now)](https://www.gioncologynow.com/post/advances-move-precision-oncology-forward-in-pancreatic-cancer)

### PDB Structures
- [1YY9: EGFR + cetuximab](https://www.rcsb.org/structure/1YY9)
- [8UKV: EGFRvIII + nanobody 34E5](https://www.rcsb.org/structure/8UKV)
- [8UKX: EGFRvIII ECD](https://www.rcsb.org/structure/8UKX)
- [1I8I: EGFRvIII peptide + MR1 dsFv](https://www.rcsb.org/structure/1I8I)
- [3B2U: EGFR + necitumumab](https://www.rcsb.org/structure/3B2U)
- [1N8Z: HER2 + trastuzumab](https://www.rcsb.org/structure/1N8Z)
- [1S78: HER2 + pertuzumab](https://www.rcsb.org/structure/1S78)
- [6OGE: HER2 + trastuzumab + pertuzumab](https://www.rcsb.org/structure/6OGE)
- [5CMA: B7-H3 ch8H9 Fab](https://www.rcsb.org/structure/5CMA)
- [4I0K: Murine B7-H3 ECD](https://www.rcsb.org/structure/4I0K)
- [9LME: B7-H3 + nanobody](https://www.rcsb.org/structure/9LME)
- [6WJL: GPC2 + D3 Fab](https://www.rcsb.org/structure/6WJL)
- [3LB6: IL-13Ra2 + IL-13](https://www.rcsb.org/structure/3LB6)
- [3VFG: Anti-GD2 3F8 Fab](https://www.rcsb.org/structure/3VFG)
- [4TUJ: Anti-GD2 14G2a Fab](https://www.rcsb.org/structure/4TUJ)
- [4F3F: MSLN + MORAb-009 Fab](https://www.rcsb.org/structure/4F3F)
- [7UED: Full-length MSLN + amatuximab Fab](https://www.rcsb.org/structure/7UED)
- [7U8C: MSLN C-term + 15B6 Fab](https://www.rcsb.org/structure/7U8C)
- [8CX3: Full-length mesothelin](https://www.rcsb.org/structure/8CX3)
- [8BW0: CEACAM5 + tusamitamab Fab](https://www.rcsb.org/structure/8BW0)
- [5IWL: CD47 + magrolimab](https://www.rcsb.org/structure/5IWL)
- [3SKJ: EphA2 + 1C1 Fab](https://www.rcsb.org/structure/3SKJ)
- [2X10: EphA2 full ECD](https://www.rcsb.org/structure/2X10)
- [5X8L: PD-L1 + atezolizumab](https://www.rcsb.org/structure/5X8L)
- [1BJ1: VEGF-A + bevacizumab](https://www.rcsb.org/structure/1BJ1)
- [7LRZ: ALK GRD domain](https://www.rcsb.org/structure/7LRZ)
- [9D20: DLK1 + ACVR2B](https://www.rcsb.org/structure/9D20)
- [8AFO: L1CAM FNIII domains](https://www.rcsb.org/structure/8AFO)
- [4K3J: MET + onartuzumab](https://www.rcsb.org/structure/4K3J)
- [5K5X: PDGFRA kinase domain](https://www.rcsb.org/structure/5K5X)
- [4FRW: Nectin-4 D1-D2](https://www.rcsb.org/structure/4FRW)
- [7PEE: TROP-2 ECD](https://www.rcsb.org/structure/7PEE)

### RFAntibody Method
- [Watson et al., De novo design of high-affinity antibodies (Nature 2025)](https://www.nature.com/articles/s41586-025-09721-5)
