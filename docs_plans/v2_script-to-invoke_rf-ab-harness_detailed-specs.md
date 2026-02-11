# v2 Spec: Batch Runner Script for De Novo Antibody Design Against Cancer Driver Targets

> **Version**: 2.0 | **Date**: 2026-02-11 | **Author**: tp53 (Ashish)
>
> Implementation-ready specification for a batch campaign runner that invokes
> rfab-harness to design de novo VHH nanobodies against 10 high-impact cancer
> targets selected from the v1 research document.
>
> **Upstream**: `docs_plans/v1_challenging-but-validated+impactful-targets-in-cancer_detailed_specs.md`

---

## Table of Contents

1. [Overview](#1-overview)
2. [Target Selection](#2-target-selection)
3. [Campaign YAML Configs](#3-campaign-yaml-configs)
4. [Batch Runner Script](#4-batch-runner-script)
5. [Epitope Derivation Protocol](#5-epitope-derivation-protocol)
6. [Compute Requirements](#6-compute-requirements)
7. [Expected Outputs](#7-expected-outputs)
8. [Implementation Checklist](#8-implementation-checklist)

---

## 1. Overview

### Goal

Create a single-command batch runner that designs de novo VHH nanobody binders
for 10 cancer targets with the highest cross-indication impact and structural
readiness for computational antibody design.

### Architecture

```
scripts/run_cancer_drivers.py          # batch orchestrator
  └── for each campaigns/cancer_drivers/*.yaml:
        └── rfab run <config.yaml>     # existing harness CLI
              ├── Stage 1: RFdiffusion   → 01_backbones.qv
              ├── Stage 2: ProteinMPNN   → 02_sequences.qv
              └── Stage 3: RF2           → 03_predictions.qv
  └── aggregate results → results/cancer_drivers_summary/
```

### Prerequisites

| Requirement | Details |
|---|---|
| rfab-harness | Installed (`pip install -e .` from repo root) |
| RFAntibody | Cloned at `./RFAntibody` (or `--rfantibody-root PATH`) |
| Python | >= 3.10 |
| GPU | NVIDIA A100/H100 with >= 40GB VRAM per campaign |
| Disk | >= 500GB for 10 campaigns at 10K designs each |

---

## 2. Target Selection

### Selection Criteria

Targets were selected from the v1 research document using three filters:

1. **Cross-indication impact**: Appears in >= 2 of 5 indications OR is the #1 target for an indication
2. **Structural readiness**: Ab-Ag co-crystal structure available in PDB (Tier 1 or Tier 2 from v1)
3. **Therapeutic potential**: High likelihood of clinical impact on mortality, recurrence, resistance, OS

### Exclusion Criteria

- Glycolipid targets (GD2) — not amenable to protein-based de novo design
- Targets without any Ab-Ag complex structure (CLDN18.2, IL-13Ra2, DLL3, ErbB3, PDGFRA)
- Targets that failed in clinical trials with no mechanistic rationale for improvement (PD-L1 in GBM)
- Targets already in the original RFAntibody paper (HA, TcdB, RSV, PHOX2B-HLA, SARS-CoV-2)

### Final 10 Targets

| # | Campaign ID | Target | Indications | PDB Complex | Res. | Why |
|---|---|---|---|---|---|---|
| 1 | `b7h3_vhh` | B7-H3 (CD276) | MPNST, DIPG, NB, GBM | 9LME (nanobody) | 2.4A | Pan-cancer 4/5; FDA Breakthrough DIPG |
| 2 | `egfrviii_vhh` | EGFRvIII | GBM | 8UKV (Nb 34E5) | — | 100% tumor-specific; absent normal tissue |
| 3 | `egfr_cetuximab_vhh` | EGFR (Domain III) | GBM, MPNST, PDAC | 1YY9 (cetuximab) | 2.6A | Cross-indication 4/5; excellent structures |
| 4 | `gpc2_vhh` | GPC2 | Neuroblastoma | 6WJL (D3 Fab) | 3.3A | Best NB differential expression + structure |
| 5 | `msln_nterm_vhh` | MSLN (N-terminal) | PDAC | 4F3F (MORAb-009) | 2.6A | 85-89% PDAC expression |
| 6 | `cd47_vhh` | CD47 | GBM | 5IWL (magrolimab) | — | Redirects TAMs; orthogonal to T-cell Tx |
| 7 | `epha2_vhh` | EphA2 | GBM | 3SKJ (1C1 Fab) | 2.5A | Dual tumor+vasculature targeting |
| 8 | `ceacam5_vhh` | CEACAM5 | PDAC | 8BW0 (tusamitamab) | 3.11A | Best clinical ORR in PDAC (20%) |
| 9 | `her2_domIV_vhh` | HER2 (Domain IV) | GBM, MPNST | 1N8Z (trastuzumab) | 2.5A | ~80% GBM expression |
| 10 | `msln_cterm_vhh` | MSLN (C-terminal) | PDAC | 7U8C (15B6 Fab) | — | Alternative epitope for bispecific |

---

## 3. Campaign YAML Configs

All configs go in `campaigns/cancer_drivers/`. Each follows the schema in
`harness/config/schema.py`. VHH format throughout (BBB penetration for CNS
tumors, stromal penetration for PDAC).

### 3.1 B7-H3 (CD276) — `campaigns/cancer_drivers/b7h3_vhh.yaml`

```yaml
campaign:
  name: "b7h3_vhh_v1"
  description: "De novo VHH targeting B7-H3/CD276 (pan-cancer: MPNST, DIPG, neuroblastoma, GBM). Epitope from 9LME nanobody complex (2.4A, Jan 2025). FDA Breakthrough for DIPG."
  version: "1.0"

target:
  pdb_id: "9LME"
  chain_id: "A"
  # B7-H3 IgV domain epitope contacted by nanobody in 9LME
  # Residues derived from 9LME interface analysis (nanobody chain B contacts on chain A)
  # IgV domain residues at antibody interface, 4.5A cutoff
  epitope_residues: [29, 30, 31, 33, 35, 52, 53, 54, 55, 56, 74, 75, 76, 77, 99, 101, 103]
  hotspot_residues: [33, 54, 76, 101]
  truncation:
    enabled: true
    buffer_angstroms: 12.0
    preserve_secondary_structure: true

antibody:
  format: "vhh"
  framework: "builtin:NbBCII10"
  cdr_loops:
    H1: "7"
    H2: "6"
    H3: "7-15"

pipeline:
  rfdiffusion:
    num_designs: 10000
    seed: 42
  proteinmpnn:
    sequences_per_backbone: 5
    temperature: 0.2
  rf2:
    recycling_iterations: 10

filtering:
  pae_threshold: 10.0
  rmsd_threshold: 2.0
  ddg_threshold: -20.0

compute:
  gpus: 1
  memory_gb: 40

output:
  directory: "./results/cancer_drivers/b7h3_vhh_v1"
  top_n_candidates: 50
  report_format: "both"
  export_pdbs: true
```

### 3.2 EGFRvIII — `campaigns/cancer_drivers/egfrviii_vhh.yaml`

```yaml
campaign:
  name: "egfrviii_vhh_v1"
  description: "De novo VHH targeting EGFRvIII-specific conformation for GBM. Template: 8UKV (nanobody 34E5). 100% tumor-specific — absent from all normal tissues."
  version: "1.0"

target:
  pdb_id: "8UKV"
  chain_id: "A"
  # EGFRvIII ECD epitope contacted by nanobody 34E5 (chain B)
  # Residues in the EGFRvIII-specific compact conformation
  epitope_residues: [287, 288, 289, 290, 291, 292, 293, 295, 297, 298, 300, 302, 310, 312, 314, 316]
  hotspot_residues: [289, 293, 302, 314]
  truncation:
    enabled: true
    buffer_angstroms: 10.0

antibody:
  format: "vhh"
  framework: "builtin:NbBCII10"
  cdr_loops:
    H1: "7"
    H2: "6"
    H3: "8-15"

pipeline:
  rfdiffusion:
    num_designs: 10000
    seed: 43

filtering:
  pae_threshold: 10.0
  rmsd_threshold: 2.0
  ddg_threshold: -20.0

compute:
  gpus: 1
  memory_gb: 40

output:
  directory: "./results/cancer_drivers/egfrviii_vhh_v1"
  top_n_candidates: 50
  report_format: "both"
  export_pdbs: true
```

### 3.3 EGFR Domain III — `campaigns/cancer_drivers/egfr_cetuximab_vhh.yaml`

```yaml
campaign:
  name: "egfr_cetuximab_vhh_v1"
  description: "De novo VHH targeting EGFR Domain III (cetuximab epitope). For GBM (57% amplified), MPNST (28% amplified), PDAC (KRAS-WT). Template: 1YY9 (cetuximab Fab, 2.6A)."
  version: "1.0"

target:
  pdb_id: "1YY9"
  chain_id: "A"
  # EGFR Domain III residues contacted by cetuximab (chains H/L in 1YY9)
  # Well-characterized epitope: Domain III loop clusters
  epitope_residues: [351, 353, 355, 357, 358, 359, 361, 382, 384, 386, 388, 389, 390, 391, 408, 410, 412, 414, 416, 465, 467, 469, 470, 471, 473, 475, 476, 489, 491, 493]
  hotspot_residues: [357, 382, 412, 467, 489]
  truncation:
    enabled: true
    buffer_angstroms: 10.0

antibody:
  format: "vhh"
  framework: "builtin:NbBCII10"
  cdr_loops:
    H1: "7"
    H2: "6"
    H3: "5-13"

pipeline:
  rfdiffusion:
    num_designs: 10000
    seed: 44

filtering:
  pae_threshold: 10.0
  rmsd_threshold: 2.0
  ddg_threshold: -20.0

compute:
  gpus: 1
  memory_gb: 40

output:
  directory: "./results/cancer_drivers/egfr_cetuximab_vhh_v1"
  top_n_candidates: 50
  report_format: "both"
  export_pdbs: true
```

### 3.4 GPC2 — `campaigns/cancer_drivers/gpc2_vhh.yaml`

```yaml
campaign:
  name: "gpc2_vhh_v1"
  description: "De novo VHH targeting GPC2 for neuroblastoma. Template: 6WJL (D3 Fab, 3.3A). Discontinuous epitope on helices 2/12/13. Best differential expression among novel NB targets."
  version: "1.0"

target:
  pdb_id: "6WJL"
  chain_id: "A"
  # GPC2 core protein epitope contacted by D3 Fab (chains H/L)
  # Discontinuous epitope spanning helices 2 (res ~44-58), 12 (~268-280), 13 (~284-296)
  # Buries 1069 A^2 with salt bridges and hydrophobic contacts
  epitope_residues: [44, 46, 48, 49, 50, 51, 53, 55, 268, 270, 272, 274, 276, 278, 284, 286, 288, 290, 292, 294]
  hotspot_residues: [48, 55, 274, 290]
  truncation:
    enabled: true
    buffer_angstroms: 12.0

antibody:
  format: "vhh"
  framework: "builtin:NbBCII10"
  cdr_loops:
    H1: "7"
    H2: "6"
    H3: "8-15"

pipeline:
  rfdiffusion:
    num_designs: 10000
    seed: 45

filtering:
  pae_threshold: 10.0
  rmsd_threshold: 2.0
  ddg_threshold: -18.0

compute:
  gpus: 1
  memory_gb: 40

output:
  directory: "./results/cancer_drivers/gpc2_vhh_v1"
  top_n_candidates: 50
  report_format: "both"
  export_pdbs: true
```

### 3.5 MSLN N-terminal — `campaigns/cancer_drivers/msln_nterm_vhh.yaml`

```yaml
campaign:
  name: "msln_nterm_vhh_v1"
  description: "De novo VHH targeting mesothelin N-terminal epitope for PDAC. Template: 4F3F (MORAb-009/amatuximab Fab, 2.6A). 85-89% PDAC expression. VHH for stromal penetration."
  version: "1.0"

target:
  pdb_id: "4F3F"
  chain_id: "A"
  # MSLN N-terminal region contacted by MORAb-009 Fab (chains H/L)
  # Amatuximab epitope on mesothelin membrane-distal region
  epitope_residues: [296, 298, 299, 300, 302, 304, 306, 310, 312, 314, 316, 318, 320, 338, 340, 342, 344, 346, 348]
  hotspot_residues: [299, 310, 318, 342]
  truncation:
    enabled: true
    buffer_angstroms: 10.0

antibody:
  format: "vhh"
  framework: "builtin:NbBCII10"
  cdr_loops:
    H1: "7"
    H2: "6"
    H3: "7-13"

pipeline:
  rfdiffusion:
    num_designs: 10000
    seed: 46

filtering:
  pae_threshold: 10.0
  rmsd_threshold: 2.0
  ddg_threshold: -20.0

compute:
  gpus: 1
  memory_gb: 40

output:
  directory: "./results/cancer_drivers/msln_nterm_vhh_v1"
  top_n_candidates: 50
  report_format: "both"
  export_pdbs: true
```

### 3.6 CD47 — `campaigns/cancer_drivers/cd47_vhh.yaml`

```yaml
campaign:
  name: "cd47_vhh_v1"
  description: "De novo VHH targeting CD47 SIRPa-competitive face for GBM. Template: 5IWL (magrolimab/5F9 diabody). Redirects TAMs (50% of GBM tumor mass)."
  version: "1.0"

target:
  pdb_id: "5IWL"
  chain_id: "C"
  # CD47 epitope contacted by 5F9/magrolimab (chains A/B in 5IWL)
  # N-terminal pyroglutamate (pGlu2), BC loop (T34-Y37), FG loop (T99-E106)
  epitope_residues: [2, 3, 4, 5, 34, 35, 36, 37, 99, 100, 101, 102, 103, 104, 106]
  hotspot_residues: [37, 101, 102, 104]
  truncation:
    enabled: true
    buffer_angstroms: 10.0

antibody:
  format: "vhh"
  framework: "builtin:NbBCII10"
  cdr_loops:
    H1: "7"
    H2: "6"
    H3: "8-15"

pipeline:
  rfdiffusion:
    num_designs: 10000
    seed: 47

filtering:
  pae_threshold: 10.0
  rmsd_threshold: 2.0
  ddg_threshold: -20.0

compute:
  gpus: 1
  memory_gb: 40

output:
  directory: "./results/cancer_drivers/cd47_vhh_v1"
  top_n_candidates: 50
  report_format: "both"
  export_pdbs: true
```

### 3.7 EphA2 — `campaigns/cancer_drivers/epha2_vhh.yaml`

```yaml
campaign:
  name: "epha2_vhh_v1"
  description: "De novo VHH targeting EphA2 for GBM. Template: 3SKJ (1C1 Fab, 2.5A). Dual tumor+vasculature expression. Design targets alternative epitope to avoid 1C1 vascular toxicity."
  version: "1.0"

target:
  pdb_id: "3SKJ"
  chain_id: "E"
  # EphA2 LBD epitope contacted by 1C1 Fab (chains H/L in 3SKJ)
  # 1C1 HCDR3 mediates ephrin ligand mimicry
  # NOTE: For safety, future iterations should explore epitopes AWAY from LBD
  epitope_residues: [43, 45, 47, 49, 51, 53, 55, 68, 70, 72, 74, 99, 101, 103, 105, 107, 109, 125, 127, 129, 131]
  hotspot_residues: [47, 72, 103, 127]
  truncation:
    enabled: true
    buffer_angstroms: 10.0

antibody:
  format: "vhh"
  framework: "builtin:NbBCII10"
  cdr_loops:
    H1: "7"
    H2: "6"
    H3: "8-15"

pipeline:
  rfdiffusion:
    num_designs: 10000
    seed: 48

filtering:
  pae_threshold: 10.0
  rmsd_threshold: 2.0
  ddg_threshold: -20.0

compute:
  gpus: 1
  memory_gb: 40

output:
  directory: "./results/cancer_drivers/epha2_vhh_v1"
  top_n_candidates: 50
  report_format: "both"
  export_pdbs: true
```

### 3.8 CEACAM5 — `campaigns/cancer_drivers/ceacam5_vhh.yaml`

```yaml
campaign:
  name: "ceacam5_vhh_v1"
  description: "De novo VHH targeting CEACAM5 A3-B3 domains for PDAC. Template: 8BW0 (tusamitamab Fab, 3.11A cryo-EM). EBC-129 ADC showed 20% ORR in PDAC (ASCO 2025)."
  version: "1.0"

target:
  pdb_id: "8BW0"
  chain_id: "C"
  # CEACAM5 A3-B3 domain epitope contacted by tusamitamab Fab (chains H/L)
  # Discontinuous epitope spanning A3 (~606-613) and B3 (~622-641) domains
  # Nature Communications 2024: epitope-paratope structural insights
  epitope_residues: [606, 607, 609, 610, 612, 613, 622, 624, 626, 628, 630, 632, 635, 636, 637, 639, 641]
  hotspot_residues: [612, 622, 636, 639]
  truncation:
    enabled: true
    buffer_angstroms: 12.0

antibody:
  format: "vhh"
  framework: "builtin:NbBCII10"
  cdr_loops:
    H1: "7"
    H2: "6"
    H3: "8-15"

pipeline:
  rfdiffusion:
    num_designs: 10000
    seed: 49

filtering:
  pae_threshold: 12.0
  rmsd_threshold: 2.5
  ddg_threshold: -18.0

compute:
  gpus: 1
  memory_gb: 40

output:
  directory: "./results/cancer_drivers/ceacam5_vhh_v1"
  top_n_candidates: 50
  report_format: "both"
  export_pdbs: true
```

### 3.9 HER2 Domain IV — `campaigns/cancer_drivers/her2_domIV_vhh.yaml`

```yaml
campaign:
  name: "her2_domIV_vhh_v1"
  description: "De novo VHH targeting HER2 Domain IV (trastuzumab epitope) for GBM (~80% expression) and MPNST. Template: 1N8Z (trastuzumab Fab, 2.5A)."
  version: "1.0"

target:
  pdb_id: "1N8Z"
  chain_id: "A"
  # HER2 Domain IV epitope contacted by trastuzumab Fab (chains H/L)
  # Three loop regions: 557-561, 570-573, 593-603 + 581-590
  epitope_residues: [557, 558, 560, 561, 567, 570, 571, 573, 574, 576, 579, 580, 581, 583, 584, 586, 590, 592, 593, 596, 597, 600, 603]
  hotspot_residues: [560, 573, 593, 597]
  truncation:
    enabled: true
    buffer_angstroms: 10.0

antibody:
  format: "vhh"
  framework: "builtin:NbBCII10"
  cdr_loops:
    H1: "7"
    H2: "6"
    H3: "8-13"

pipeline:
  rfdiffusion:
    num_designs: 10000
    seed: 50

filtering:
  pae_threshold: 10.0
  rmsd_threshold: 2.0
  ddg_threshold: -20.0

compute:
  gpus: 1
  memory_gb: 40

output:
  directory: "./results/cancer_drivers/her2_domIV_vhh_v1"
  top_n_candidates: 50
  report_format: "both"
  export_pdbs: true
```

### 3.10 MSLN C-terminal — `campaigns/cancer_drivers/msln_cterm_vhh.yaml`

```yaml
campaign:
  name: "msln_cterm_vhh_v1"
  description: "De novo VHH targeting mesothelin C-terminal juxtamembrane epitope for PDAC. Template: 7U8C (15B6 Fab). Alternative epitope to 4F3F for potential bispecific pairing."
  version: "1.0"

target:
  pdb_id: "7U8C"
  chain_id: "A"
  # MSLN C-terminal juxtamembrane epitope contacted by 15B6 Fab (chains H/L)
  # Linear epitope: N584-S598, core: Y586-L597
  epitope_residues: [584, 585, 586, 587, 588, 589, 590, 591, 592, 593, 594, 595, 596, 597, 598]
  hotspot_residues: [586, 589, 592, 594]
  truncation:
    enabled: false

antibody:
  format: "vhh"
  framework: "builtin:NbBCII10"
  cdr_loops:
    H1: "7"
    H2: "6"
    H3: "8-15"

pipeline:
  rfdiffusion:
    num_designs: 10000
    seed: 51

filtering:
  pae_threshold: 10.0
  rmsd_threshold: 2.0
  ddg_threshold: -18.0

compute:
  gpus: 1
  memory_gb: 40

output:
  directory: "./results/cancer_drivers/msln_cterm_vhh_v1"
  top_n_candidates: 50
  report_format: "both"
  export_pdbs: true
```

---

## 4. Batch Runner Script

### File: `scripts/run_cancer_drivers.py`

### 4.1 CLI Interface

```
python scripts/run_cancer_drivers.py [OPTIONS]

Options:
  --config-dir PATH       Campaign config directory (default: campaigns/cancer_drivers/)
  --rfantibody-root PATH  Path to RFAntibody repo (default: ./RFAntibody)
  --parallel N            Run N campaigns concurrently (default: 1 = sequential)
  --dry-run               Validate configs and prepare inputs only
  --continue-on-error     Continue remaining campaigns if one fails
  --campaigns LIST        Comma-separated campaign names to run (default: all)
  --summary-dir PATH      Cross-campaign summary output (default: results/cancer_drivers_summary/)
  -v, --verbose           Verbose logging
```

### 4.2 Script Behavior

```python
# Pseudocode for the batch runner

def main():
    args = parse_args()

    # Phase 1: Discovery
    configs = discover_configs(args.config_dir)  # glob *.yaml, sort alphabetically
    log(f"Found {len(configs)} campaign configs")

    # Phase 2: Validate ALL before running ANY
    for config_path in configs:
        result = subprocess.run(["rfab", "validate", str(config_path)])
        if result.returncode != 0:
            if not args.continue_on_error:
                fatal(f"Validation failed: {config_path}")
            else:
                skip(config_path)

    # Phase 3: Compute estimate
    total_designs = sum(load_num_designs(c) for c in configs)
    estimate_gpu_hours(total_designs)

    # Phase 4: Execute campaigns
    if args.parallel > 1:
        run_parallel(configs, max_workers=args.parallel)
    else:
        run_sequential(configs)

    # Phase 5: Aggregate results
    aggregate_cross_campaign(configs, args.summary_dir)
```

### 4.3 Core Functions to Implement

#### `discover_configs(config_dir: Path) -> list[Path]`

- Glob `config_dir/*.yaml`
- Sort alphabetically by filename
- Return list of Path objects
- If `--campaigns` filter is set, filter to matching names

#### `run_sequential(configs: list[Path], args) -> dict[str, CampaignResult]`

```python
results = {}
for config_path in configs:
    name = config_path.stem
    log(f"[{i+1}/{len(configs)}] Starting campaign: {name}")
    start = time.monotonic()

    cmd = ["rfab", "run", str(config_path)]
    if args.verbose:
        cmd.append("-v")
    if args.rfantibody_root:
        cmd.extend(["--rfantibody-root", str(args.rfantibody_root)])

    proc = subprocess.run(cmd, capture_output=True, text=True)

    elapsed = time.monotonic() - start
    results[name] = CampaignResult(
        name=name,
        config_path=config_path,
        returncode=proc.returncode,
        elapsed_seconds=elapsed,
        stdout=proc.stdout,
        stderr=proc.stderr,
    )

    if proc.returncode != 0 and not args.continue_on_error:
        log_error(f"Campaign {name} failed after {elapsed:.0f}s")
        break

    log(f"[{i+1}/{len(configs)}] Completed {name} in {elapsed/3600:.1f}h")

return results
```

#### `run_parallel(configs: list[Path], max_workers: int, args) -> dict[str, CampaignResult]`

```python
from concurrent.futures import ProcessPoolExecutor, as_completed

# Assign GPU IDs round-robin
gpu_assignments = {c.stem: i % max_workers for i, c in enumerate(configs)}

with ProcessPoolExecutor(max_workers=max_workers) as pool:
    futures = {}
    for config_path in configs:
        gpu_id = gpu_assignments[config_path.stem]
        env = {**os.environ, "CUDA_VISIBLE_DEVICES": str(gpu_id)}
        futures[pool.submit(run_single, config_path, args, env)] = config_path.stem

    for future in as_completed(futures):
        name = futures[future]
        result = future.result()
        log(f"Completed: {name} ({'OK' if result.returncode == 0 else 'FAILED'})")
```

#### `aggregate_cross_campaign(configs, summary_dir: Path)`

```python
summary_dir.mkdir(parents=True, exist_ok=True)

all_candidates = []
for config_path in configs:
    config = load_config(config_path)
    candidates_csv = Path(config.output.directory) / "candidates" / "summary.csv"
    if candidates_csv.exists():
        df = pd.read_csv(candidates_csv)
        df["campaign"] = config.campaign.name
        df["target"] = config_path.stem
        all_candidates.append(df)

if all_candidates:
    combined = pd.concat(all_candidates, ignore_index=True)
    combined.sort_values("composite_score", ascending=True, inplace=True)
    combined.to_csv(summary_dir / "cross_campaign_comparison.csv", index=False)

    # Per-campaign summary stats
    stats = combined.groupby("campaign").agg(
        total=("design_id", "count"),
        best_pae=("pae", "min"),
        best_rmsd=("rmsd", "min"),
        best_ddg=("ddg", "min"),
        median_composite=("composite_score", "median"),
    )
    stats.to_csv(summary_dir / "campaign_stats.csv")

    generate_summary_report(combined, stats, summary_dir / "report.html")
```

#### `generate_summary_report(candidates_df, stats_df, output_path: Path)`

HTML report containing:
- Campaign completion status table (pass/fail, elapsed time)
- Per-campaign hit rates (designs passing all filters / total designs)
- Score distribution box plots across campaigns
- Top 10 candidates across all campaigns (with PDB links and scores)
- Heatmap of target × metric (best pAE, RMSD, ddG per campaign)

### 4.4 Data Classes

```python
@dataclass
class CampaignResult:
    name: str
    config_path: Path
    returncode: int
    elapsed_seconds: float
    stdout: str
    stderr: str

    @property
    def success(self) -> bool:
        return self.returncode == 0
```

### 4.5 Logging

- Log to both stdout and `results/cancer_drivers_summary/batch_run.log`
- Timestamps in ISO 8601
- Per-campaign log lines: `[HH:MM:SS] [campaign_name] Stage X: ...`
- Final summary table printed to stdout on completion

---

## 5. Epitope Derivation Protocol

### Method

For each target, epitope residues were derived by one of two methods:

**Method A — Published epitope data**: For well-characterized complexes (EGFR/cetuximab, HER2/trastuzumab, CD47/magrolimab, CEACAM5/tusamitamab, MSLN C-term/15B6), epitope residues are taken directly from published structural studies.

**Method B — PDB interface analysis**: For newer complexes (9LME, 8UKV, 6WJL, 4F3F, 3SKJ), epitope residues are identified using a 4.5A distance cutoff between target and antibody heavy atoms.

### Hotspot Selection Criteria

From the epitope residues, 3-5 hotspot residues are selected based on:
1. **Burial**: Highest solvent-accessible surface area loss upon complex formation
2. **Hydrophobicity**: Preference for Phe, Trp, Tyr, Leu, Ile, Val at interface center
3. **Hydrogen bonding**: Residues forming >= 2 intermolecular H-bonds
4. **Spatial centrality**: Residues near the geometric center of the epitope patch

### Computational Verification Script

The implementation should include a helper script (`scripts/derive_epitopes.py`) that:

1. Downloads PDB files for all 10 complexes
2. For each complex, identifies target vs. antibody chains
3. Computes contact residues at 4.5A heavy-atom distance cutoff
4. Ranks contacts by burial and interaction energy
5. Selects top 3-5 as hotspots
6. Outputs residue lists to verify against the hardcoded YAML configs

```
python scripts/derive_epitopes.py --pdb-ids 9LME,8UKV,1YY9,6WJL,4F3F,5IWL,3SKJ,8BW0,1N8Z,7U8C
```

This script uses BioPython's `PDB.NeighborSearch` for distance calculations. It serves as
a verification step — if its output conflicts with the YAML configs, the YAML should be
updated before running campaigns.

### Per-Target Epitope Derivation Detail

| Campaign | Method | Source | Target Chain | Ab Chain(s) | Notes |
|---|---|---|---|---|---|
| b7h3_vhh | B | 9LME interface | A | B (Nb) | IgV domain face |
| egfrviii_vhh | B | 8UKV interface | A | B (Nb 34E5) | EGFRvIII-specific conformation around res 287-302 |
| egfr_cetuximab_vhh | A | Li et al. 2005 (1YY9 paper) | A | H, L | Domain III loops |
| gpc2_vhh | A+B | Cell Reports Med 2021 | A | H, L (D3 Fab) | Discontinuous helices 2/12/13 |
| msln_nterm_vhh | B | 4F3F interface | A | H, L | Membrane-distal region |
| cd47_vhh | A | Weiskopf et al. 2016 | C | A, B (diabody) | pGlu2 + BC/FG loops |
| epha2_vhh | B | 3SKJ interface | E | H, L (1C1 Fab) | LBD; HCDR3 ephrin mimicry |
| ceacam5_vhh | A | Nat Comms 2024 | C | H, L | A3-B3 interdomain epitope |
| her2_domIV_vhh | A | Cho et al. 2003 (1N8Z) | A | H, L | Domain IV loops 1-3 |
| msln_cterm_vhh | A | Ho et al. 2022 (PNAS) | A | H, L (15B6 Fab) | Linear juxtamembrane Y586-L597 |

---

## 6. Compute Requirements

### Per-Campaign Estimates

Based on RFAntibody benchmarks (Watson et al. 2025, Extended Data Table 2):

| Stage | Time per Design (A100) | 10K Designs | Output Size |
|---|---|---|---|
| RFdiffusion | ~3-5 sec | ~10-14 hrs | ~20 GB (.qv) |
| ProteinMPNN | ~0.5 sec per seq × 5 seqs | ~7 hrs | ~25 GB (.qv) |
| RF2 | ~8-12 sec per prediction | ~22-33 hrs | ~30 GB (.qv) |
| **Total per campaign** | | **~40-54 hrs** | **~75 GB** |

### 10-Campaign Totals

| Mode | GPUs | Wall Time | Disk |
|---|---|---|---|
| Sequential (1 GPU) | 1× A100 | ~400-540 hrs (17-23 days) | ~750 GB |
| Parallel (4 GPUs) | 4× A100 | ~100-135 hrs (4-6 days) | ~750 GB |
| Parallel (8 GPUs) | 8× A100 | ~50-68 hrs (2-3 days) | ~750 GB |

### Reduced-Scale Option

For initial validation, use `num_designs: 1000` (10× reduction):

| Mode | GPUs | Wall Time | Disk |
|---|---|---|---|
| Sequential | 1× A100 | ~40-54 hrs (2 days) | ~75 GB |
| Parallel (4 GPUs) | 4× A100 | ~10-14 hrs | ~75 GB |

### Memory Requirements

- RFdiffusion: ~20-30 GB VRAM (A100 40GB sufficient)
- ProteinMPNN: ~8-12 GB VRAM
- RF2: ~16-24 GB VRAM
- System RAM: >= 32 GB per concurrent campaign

---

## 7. Expected Outputs

### Per-Campaign Output Tree

```
results/cancer_drivers/<campaign_name>/
├── prep/
│   ├── <PDB_ID>.pdb              # Downloaded structure
│   ├── <PDB_ID>_chain<X>.pdb     # Extracted target chain
│   ├── target_truncated.pdb      # Truncated around epitope
│   └── framework.pdb             # HLT-format VHH framework
├── pipeline/
│   ├── 01_backbones.qv           # RFdiffusion backbones
│   ├── 02_sequences.qv           # ProteinMPNN sequences
│   ├── 03_predictions.qv         # RF2 structure predictions
│   └── .checkpoint_stage{1,2,3}  # Resume checkpoints
├── analysis/
│   ├── ranked_candidates.csv     # All passing candidates, ranked
│   └── report.html               # Score distributions, top hits
└── candidates/
    ├── design_001.pdb             # Top 50 candidate structures
    ├── ...
    ├── design_050.pdb
    ├── summary.csv                # Scores for exported candidates
    └── sequences.fasta            # CDR sequences
```

### Cross-Campaign Summary

```
results/cancer_drivers_summary/
├── cross_campaign_comparison.csv  # All candidates across all 10 campaigns
├── campaign_stats.csv             # Per-campaign summary metrics
├── report.html                    # Interactive comparison report
└── batch_run.log                  # Execution log
```

### Key Metrics in `cross_campaign_comparison.csv`

| Column | Description |
|---|---|
| `campaign` | Campaign name (e.g., b7h3_vhh_v1) |
| `target` | Target protein name |
| `design_id` | Unique design identifier |
| `pae` | Predicted Aligned Error (lower = better) |
| `rmsd` | Backbone RMSD to input (lower = better) |
| `ddg` | Binding free energy estimate (more negative = better) |
| `composite_score` | Weighted rank: 0.4×pAE + 0.3×RMSD + 0.3×ddG (min-max normalized) |

### Success Criteria

| Metric | Threshold | Interpretation |
|---|---|---|
| Pass rate | >= 1% of designs pass all filters | Campaign is structurally viable |
| Best pAE | < 5.0 | High-confidence interface prediction |
| Best ddG | < -30.0 | Strong predicted binding affinity |
| Diversity | >= 10 unique CDR3 sequences in top 50 | Sufficient sequence diversity for YSD |

---

## 8. Implementation Checklist

### Files to Create

| # | Path | Type | Description |
|---|---|---|---|
| 1 | `campaigns/cancer_drivers/b7h3_vhh.yaml` | Config | B7-H3 campaign |
| 2 | `campaigns/cancer_drivers/egfrviii_vhh.yaml` | Config | EGFRvIII campaign |
| 3 | `campaigns/cancer_drivers/egfr_cetuximab_vhh.yaml` | Config | EGFR Domain III campaign |
| 4 | `campaigns/cancer_drivers/gpc2_vhh.yaml` | Config | GPC2 campaign |
| 5 | `campaigns/cancer_drivers/msln_nterm_vhh.yaml` | Config | MSLN N-terminal campaign |
| 6 | `campaigns/cancer_drivers/cd47_vhh.yaml` | Config | CD47 campaign |
| 7 | `campaigns/cancer_drivers/epha2_vhh.yaml` | Config | EphA2 campaign |
| 8 | `campaigns/cancer_drivers/ceacam5_vhh.yaml` | Config | CEACAM5 campaign |
| 9 | `campaigns/cancer_drivers/her2_domIV_vhh.yaml` | Config | HER2 Domain IV campaign |
| 10 | `campaigns/cancer_drivers/msln_cterm_vhh.yaml` | Config | MSLN C-terminal campaign |
| 11 | `scripts/run_cancer_drivers.py` | Script | Batch runner |
| 12 | `scripts/derive_epitopes.py` | Script | Epitope verification helper |

### Existing Files Referenced (do not modify)

| Path | Purpose |
|---|---|
| `harness/cli.py` | `rfab` CLI entry point |
| `harness/config/schema.py` | `CampaignConfig` dataclass + `load_config()` |
| `harness/pipeline/orchestrator.py` | 3-stage pipeline with checkpointing |
| `harness/analysis/filter.py` | Score extraction + threshold filtering |
| `harness/analysis/rank.py` | Composite score ranking |
| `harness/analysis/report.py` | HTML/CSV report generation |
| `harness/analysis/export.py` | PDB/FASTA candidate export |

### Validation Steps

1. `rfab validate campaigns/cancer_drivers/*.yaml` — all 10 pass
2. `rfab run campaigns/cancer_drivers/b7h3_vhh.yaml --dry-run` — prep succeeds (needs PDB download)
3. `python scripts/derive_epitopes.py` — epitope residues match YAML configs
4. `python scripts/run_cancer_drivers.py --dry-run` — discovers all 10 configs, validates, prints compute estimate
5. `python scripts/run_cancer_drivers.py --campaigns b7h3_vhh --rfantibody-root ./RFAntibody` — single campaign end-to-end
6. `ls results/cancer_drivers_summary/` — cross-campaign summary generated

### Dependencies (add to pyproject.toml / requirements.txt if not present)

```
pandas>=2.0
```

No new dependencies required beyond what rfab-harness already uses. The batch runner script uses only stdlib (`subprocess`, `concurrent.futures`, `pathlib`, `argparse`, `logging`, `time`, `json`, `os`) plus `pandas` (already a harness dependency) and `yaml` (already a harness dependency).
