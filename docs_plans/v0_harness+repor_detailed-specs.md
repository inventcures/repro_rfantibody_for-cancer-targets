# RFAntibody Harness + Reproduction: Detailed Computational Specifications

**Version**: 0.1 (Draft)
**Date**: 2026-02-11
**Source Paper**: Bennett et al., "Atomically accurate de novo design of antibodies with RFdiffusion", Nature (2025)
**Upstream Repo**: https://github.com/RosettaCommons/RFantibody

---

## Table of Contents

1. [Executive Summary](#1-executive-summary)
2. [RFAntibody Pipeline Reference](#2-rfantibody-pipeline-reference)
3. [Harness Architecture](#3-harness-architecture)
4. [Campaign Configuration Schema](#4-campaign-configuration-schema)
5. [Target Preparation Module](#5-target-preparation-module)
6. [Pipeline Orchestration Engine](#6-pipeline-orchestration-engine)
7. [Results Analysis Module](#7-results-analysis-module)
8. [Paper Reproduction Campaigns](#8-paper-reproduction-campaigns)
9. [New Cancer Target Campaigns](#9-new-cancer-target-campaigns)
10. [Rare Disease Target Campaigns](#10-rare-disease-target-campaigns)
11. [Experimental Validation Pipeline](#11-experimental-validation-pipeline)
12. [Verification & Validation](#12-verification--validation)
13. [Project Directory Structure](#13-project-directory-structure)
14. [Implementation Roadmap](#14-implementation-roadmap)

---

## 1. Executive Summary

### Problem

RFAntibody is a powerful 3-stage pipeline for de novo antibody design, but it currently requires manual orchestration of individual scripts, manual input preparation, and manual result analysis. There is no integrated way to define a "design campaign" (target + epitope + parameters) and have the system produce ranked antibody candidates end-to-end.

### Solution

Build a Python harness that wraps the existing RFAntibody pipeline and provides:

1. **Campaign-as-config**: Define a complete design campaign in a single YAML file
2. **Automated target prep**: Fetch PDBs, truncate targets, convert frameworks, validate inputs
3. **Pipeline orchestration**: Run all 3 stages with Quiver I/O, GPU parallelization, checkpointing
4. **Results analysis**: Filter, rank, and report on candidate designs
5. **Pre-built campaigns**: Reproduce paper results + test new cancer/rare disease targets
6. **Experimental planning**: Bridge to wet-lab validation workflows

### Design Principles

- **Wrap, don't reimplement**: The harness calls existing RFAntibody scripts as subprocesses
- **Minimal dependencies**: Python stdlib + PyYAML + pandas + BioPython (already in RFAntibody env)
- **Clear separation**: `harness/` directory sits alongside the cloned RFAntibody repo
- **Reproducible**: Every campaign is fully defined by its YAML config + git SHA

---

## 2. RFAntibody Pipeline Reference

### 2.1 Three-Stage Pipeline

```
Target PDB + Framework HLT + Hotspots + CDR Ranges
                    |
                    v
    ┌───────────────────────────────┐
    │  Stage 1: RFdiffusion         │
    │  Backbone Generation          │
    │  SE3-equivariant diffusion    │
    │  Output: N backbone designs   │
    └───────────────┬───────────────┘
                    |
                    v
    ┌───────────────────────────────┐
    │  Stage 2: ProteinMPNN         │
    │  Sequence Design              │
    │  CDR-specific masking         │
    │  Output: N sequenced designs  │
    └───────────────┬───────────────┘
                    |
                    v
    ┌───────────────────────────────┐
    │  Stage 3: RF2                 │
    │  Structure Prediction         │
    │  Filtering (pAE/RMSD/ddG)    │
    │  Output: Ranked candidates    │
    └───────────────────────────────┘
```

### 2.2 Upstream CLI Reference

**Stage 1 - RFdiffusion:**
```bash
rfdiffusion \
  -t target.pdb \
  -f framework.pdb \
  -q designs.qv \
  -n 10000 \
  -l "H1:7,H2:6,H3:5-13" \
  -h "T305,T456"
```

Hydra-style alternative:
```bash
python src/rfantibody/scripts/rfdiffusion_inference.py \
  --config-name antibody \
  antibody.target_pdb=/path/to/target.pdb \
  antibody.framework_pdb=/path/to/framework.pdb \
  inference.ckpt_override_path=/path/to/weights/RFdiffusion_Ab.pt \
  'ppi.hotspot_res=[T305,T456]' \
  'antibody.design_loops=[L1:8-13,L2:7,L3:9-11,H1:7,H2:6,H3:5-13]' \
  inference.num_designs=10000 \
  inference.quiver=/path/to/output.qv
```

**Stage 2 - ProteinMPNN:**
```bash
proteinmpnn \
  -q backbones.qv \
  --output-quiver sequences.qv \
  -n 5 \
  -t 0.2
```

Script-style:
```bash
python scripts/proteinmpnn_interface_design.py \
  -inquiver /path/to/backbones.qv \
  -outquiver /path/to/sequences.qv
```

**Stage 3 - RF2:**
```bash
rf2 \
  -q sequences.qv \
  --output-quiver predictions.qv \
  -r 10
```

Script-style:
```bash
python scripts/rf2_predict.py \
  input.quiver=/path/to/sequences.qv \
  output.quiver=/path/to/predictions.qv
```

### 2.3 Key File Formats

**HLT Format**: PDB files with chains labeled H (heavy), L (light), T (target). Chain order: H → L → T. CDR annotations via PDB REMARK lines:
```
REMARK PDBinfo-LABEL: 32 H1
REMARK PDBinfo-LABEL: 52 H2
REMARK PDBinfo-LABEL: 78 H3
```

**Quiver Format**: Single-file container for thousands of PDB designs. Text format with `QV_TAG` delimiters and `QV_SCORE` metadata. Eliminates filesystem strain from thousands of individual PDB files.

Quiver tools: `qvls`, `qvextract`, `qvextractspecific`, `qvscorefile`, `qvsplit`, `qvslice`, `qvrename`, `qvfrompdbs`

### 2.4 Filtering Thresholds (from paper + README)

| Metric | Threshold | Source |
|--------|-----------|--------|
| pAE (predicted aligned error) | < 10 | Paper + README |
| RMSD (design vs predicted) | < 2 Å | Paper + README |
| ddG (binding energy) | < -20 REU | README (optional) |
| iPTM (AF3, retrospective) | > 0.6 | Paper |
| pBind (RF2) | > 0.99 | Paper |

**Caveat from README**: "The lack of an effective filter is the main limitation of the RFantibody pipeline at the moment."

### 2.5 Model Weights

Downloaded via `include/download_weights.sh`:
- `weights/RFdiffusion_Ab.pt` - Antibody-finetuned RFdiffusion
- `weights/` - ProteinMPNN weights (base, not finetuned)
- `weights/` - Antibody-finetuned RF2

---

## 3. Harness Architecture

### 3.1 High-Level Architecture

```
harness/
├── cli.py                    # CLI entry point: `rfab campaign run config.yaml`
├── config/
│   ├── schema.py             # Campaign config validation (Pydantic or dataclass)
│   └── defaults.py           # Default parameter values
├── target_prep/
│   ├── fetch_pdb.py          # Download PDB from RCSB
│   ├── truncate.py           # Truncate target around epitope
│   ├── convert_framework.py  # Chothia → HLT conversion
│   └── validate.py           # Input validation
├── pipeline/
│   ├── orchestrator.py       # Main pipeline runner
│   ├── rfdiffusion.py        # Stage 1 wrapper
│   ├── proteinmpnn.py        # Stage 2 wrapper
│   ├── rf2.py                # Stage 3 wrapper
│   └── parallel.py           # Multi-GPU parallelization via qvsplit
├── analysis/
│   ├── filter.py             # Apply filtering thresholds
│   ├── rank.py               # Rank candidates by composite score
│   ├── report.py             # Generate HTML/CSV reports
│   └── export.py             # Export top candidates (PDB files)
├── campaigns/                # Pre-built campaign YAML configs
│   ├── paper_repro/          # Paper reproduction campaigns
│   │   ├── ha_vhh.yaml
│   │   ├── tcdb_vhh.yaml
│   │   ├── tcdb_scfv.yaml
│   │   ├── rsv_vhh.yaml
│   │   ├── phox2b_scfv.yaml
│   │   └── sars_cov2_vhh.yaml
│   ├── cancer/               # Cancer target campaigns
│   │   ├── pdl1_vhh.yaml
│   │   ├── her2_vhh.yaml
│   │   ├── egfr_vhh.yaml
│   │   ├── ctla4_vhh.yaml
│   │   ├── trop2_vhh.yaml
│   │   └── ...
│   └── rare_disease/         # Rare disease campaigns
│       ├── gne_stabilizer_vhh.yaml
│       ├── complement_c5_vhh.yaml
│       └── il6r_vhh.yaml
└── experimental/             # Experimental validation planning
    ├── yeast_display.py      # Yeast display screening plan generator
    ├── spr_protocol.py       # SPR binding protocol generator
    └── maturation.py         # OrthoRep maturation planning
```

### 3.2 Key Interfaces

```python
# Core campaign runner
class CampaignRunner:
    def __init__(self, config: CampaignConfig, rfantibody_root: Path)
    def prepare_inputs(self) -> PreparedInputs
    def run_pipeline(self, inputs: PreparedInputs) -> PipelineResults
    def analyze_results(self, results: PipelineResults) -> AnalysisReport
    def run(self) -> AnalysisReport  # end-to-end

# Config loader
class CampaignConfig:
    target: TargetConfig          # PDB, epitope, hotspots
    antibody: AntibodyConfig      # Format (VHH/scFv), framework, CDR ranges
    pipeline: PipelineConfig      # Num designs, sequences per backbone, recycling
    filtering: FilteringConfig    # pAE, RMSD, ddG thresholds
    compute: ComputeConfig        # GPU count, parallelization
    output: OutputConfig          # Output directory, report format

# Stage wrappers (each calls upstream script via subprocess)
class RFdiffusionRunner:
    def run(self, target_pdb: Path, framework_pdb: Path, output_qv: Path,
            hotspots: list[str], design_loops: list[str], num_designs: int) -> Path

class ProteinMPNNRunner:
    def run(self, input_qv: Path, output_qv: Path,
            seqs_per_backbone: int, temperature: float) -> Path

class RF2Runner:
    def run(self, input_qv: Path, output_qv: Path,
            recycling_iterations: int) -> Path
```

### 3.3 Data Flow

```
Campaign YAML
    │
    ├──> TargetPrep
    │       ├── fetch_pdb(pdb_id) → raw.pdb
    │       ├── truncate(raw.pdb, epitope_residues, buffer=10Å) → truncated.pdb
    │       └── validate(truncated.pdb, hotspots)
    │
    ├──> FrameworkPrep
    │       ├── select_framework(antibody_format) → chothia.pdb
    │       └── convert_to_hlt(chothia.pdb) → framework.hlt
    │
    ├──> Stage 1: RFdiffusion
    │       └── run(truncated.pdb, framework.hlt, hotspots, cdr_loops, n=10000) → backbones.qv
    │
    ├──> Stage 2: ProteinMPNN
    │       └── run(backbones.qv, seqs_per_backbone=5, temp=0.2) → sequences.qv
    │
    ├──> Stage 3: RF2
    │       └── run(sequences.qv, recycling=10) → predictions.qv
    │
    └──> Analysis
            ├── extract_scores(predictions.qv) → scores.tsv
            ├── filter(scores.tsv, pAE<10, RMSD<2, ddG<-20) → filtered.tsv
            ├── rank(filtered.tsv) → ranked.tsv
            ├── export_top(predictions.qv, ranked.tsv, n=50) → top_candidates/
            └── generate_report(ranked.tsv) → report.html
```

---

## 4. Campaign Configuration Schema

### 4.1 Full YAML Schema

```yaml
# campaign.yaml
campaign:
  name: "pdl1_vhh_v1"
  description: "De novo VHH design targeting PD-L1 immune checkpoint"
  version: "1.0"

target:
  pdb_id: "5N2C"                    # RCSB PDB ID (fetched automatically)
  # pdb_file: "/path/to/custom.pdb" # OR provide a local file
  chain_id: "A"                     # Target chain in the PDB
  epitope_residues: [54, 56, 58, 60, 62, 66, 68, 113, 115, 117]  # Residue numbers
  hotspot_residues: [56, 60, 115]   # Subset of epitope - hydrophobic anchors
  truncation:
    enabled: true
    buffer_angstroms: 10.0          # Keep residues within 10Å of epitope
    preserve_secondary_structure: true

antibody:
  format: "vhh"                     # "vhh" (nanobody) or "scfv"
  framework: "builtin:NbBCII10"     # Built-in framework or path to PDB
  # framework: "/path/to/custom_framework.pdb"  # Custom Chothia-annotated framework
  cdr_loops:
    H1: "7"                         # Fixed 7 residues
    H2: "6"                         # Fixed 6 residues
    H3: "5-13"                      # Variable 5-13 residues
    # For scFv, also specify:
    # L1: "8-13"
    # L2: "7"
    # L3: "9-11"

pipeline:
  rfdiffusion:
    num_designs: 10000
    weights: "default"              # "default" uses bundled weights
    seed: null                      # null = random, or integer for reproducibility
  proteinmpnn:
    sequences_per_backbone: 5
    temperature: 0.2
  rf2:
    recycling_iterations: 10
    hotspot_percentage: 0.1

filtering:
  pae_threshold: 10.0
  rmsd_threshold: 2.0
  ddg_threshold: -20.0              # Optional, set null to skip
  # Future: af3_iptm_threshold: 0.6

compute:
  gpus: 1                           # Number of GPUs for parallelization
  # If gpus > 1, pipeline splits Quiver files via qvsplit
  memory_gb: 16
  container: "docker"               # "docker", "apptainer", or "local"

output:
  directory: "./results/pdl1_vhh_v1"
  top_n_candidates: 50
  report_format: "html"             # "html", "csv", or "both"
  export_pdbs: true
  keep_intermediates: false         # Keep intermediate Quiver files

experimental:
  enabled: false                    # Set true to generate experimental protocols
  screening_method: "yeast_display"
  target_concentrations: [1400, 78]  # nM, for flow cytometry
  spr_concentrations: [0.1, 1, 10, 100, 1000]  # nM
```

### 4.2 Built-in Frameworks

| ID | Name | Format | Source |
|----|------|--------|--------|
| `NbBCII10` | h-NbBCII10 | VHH (nanobody) | `scripts/examples/example_inputs/h-NbBCII10.pdb` |
| `hu4D5-8` | hu-4D5-8 Fv | scFv | `scripts/examples/example_inputs/hu-4D5-8_Fv.pdb` |

### 4.3 Config Validation Rules

- `target.pdb_id` XOR `target.pdb_file` must be set (not both)
- `target.hotspot_residues` must be subset of `target.epitope_residues`
- `target.hotspot_residues` must contain >= 3 hydrophobic residues
- `antibody.format` = "vhh" → only H1/H2/H3 CDR loops allowed
- `antibody.format` = "scfv" → must specify both H and L CDR loops
- `pipeline.rfdiffusion.num_designs` >= 20 (minimum for pilot)
- CDR loop ranges must be within biologically plausible bounds

---

## 5. Target Preparation Module

### 5.1 PDB Fetching

```python
def fetch_pdb(pdb_id: str, output_dir: Path) -> Path:
    """Download PDB structure from RCSB.

    Uses RCSB REST API: https://data.rcsb.org/rest/v1/core/entry/{pdb_id}
    Downloads mmCIF or PDB format.
    Extracts specified chain(s).
    """
```

### 5.2 Target Truncation

The paper recommends truncating large targets to reduce O(N^2) computational cost. Keep residues within ~10 Å of the epitope, preserving secondary structure elements.

```python
def truncate_target(
    pdb_path: Path,
    epitope_residues: list[int],
    chain_id: str,
    buffer_angstroms: float = 10.0,
    preserve_secondary_structure: bool = True,
    output_path: Path = None
) -> Path:
    """Truncate target to residues near epitope.

    1. Load PDB via BioPython
    2. Identify residues within buffer of any epitope residue (Cα distance)
    3. If preserve_secondary_structure: extend selection to complete SS elements
    4. Write truncated PDB with chain renamed to 'T'
    """
```

### 5.3 Hotspot Validation

From the paper: hotspot residues are defined as antigen residues with Cβ distance < 8 Å to nearest CDR residues. For initial design, the user specifies them. The harness validates they are surface-exposed and contain hydrophobic residues.

```python
def validate_hotspots(
    target_pdb: Path,
    hotspot_residues: list[int],
    chain_id: str
) -> ValidationResult:
    """Validate hotspot selection.

    Checks:
    - All residue numbers exist in the structure
    - Residues are surface-exposed (SASA > threshold)
    - At least 3 hydrophobic residues (A, V, L, I, M, F, W, P, Y)
    - Residues form a contiguous patch (not scattered)
    - Warns if heavily glycosylated region detected
    """
```

### 5.4 Framework Conversion

Wraps existing `scripts/util/chothia_to_HLT.py`:

```python
def prepare_framework(
    framework_source: str,  # "builtin:NbBCII10" or "/path/to/chothia.pdb"
    antibody_format: str,   # "vhh" or "scfv"
    output_path: Path
) -> Path:
    """Prepare antibody framework in HLT format.

    1. If builtin: copy from examples directory
    2. If custom: run chothia_to_HLT.py conversion
    3. Validate output has correct chain labels (H, L, T)
    4. Validate REMARK CDR annotations present
    """
```

### 5.5 Input Assembly

```python
def assemble_inputs(
    target_pdb: Path,       # Truncated target
    framework_hlt: Path,    # HLT-formatted framework
    hotspot_residues: list[int],
    chain_id: str
) -> PreparedInputs:
    """Assemble validated inputs for pipeline.

    1. Merge target into framework HLT as 'T' chain
    2. Generate hotspot string: "T305,T456,T489"
    3. Validate merged structure integrity
    4. Return PreparedInputs dataclass
    """
```

---

## 6. Pipeline Orchestration Engine

### 6.1 Sequential Stage Execution

```python
class PipelineOrchestrator:
    def __init__(self, config: CampaignConfig, rfantibody_root: Path):
        self.config = config
        self.root = rfantibody_root
        self.work_dir = Path(config.output.directory) / "pipeline"

    def run(self, inputs: PreparedInputs) -> PipelineResults:
        """Run all 3 stages sequentially with checkpoint support."""

        # Stage 1: RFdiffusion
        backbones_qv = self.work_dir / "01_backbones.qv"
        if not self._checkpoint_exists("stage1"):
            self._run_rfdiffusion(inputs, backbones_qv)
            self._save_checkpoint("stage1")

        # Stage 2: ProteinMPNN
        sequences_qv = self.work_dir / "02_sequences.qv"
        if not self._checkpoint_exists("stage2"):
            self._run_proteinmpnn(backbones_qv, sequences_qv)
            self._save_checkpoint("stage2")

        # Stage 3: RF2
        predictions_qv = self.work_dir / "03_predictions.qv"
        if not self._checkpoint_exists("stage3"):
            self._run_rf2(sequences_qv, predictions_qv)
            self._save_checkpoint("stage3")

        return PipelineResults(
            backbones=backbones_qv,
            sequences=sequences_qv,
            predictions=predictions_qv
        )
```

### 6.2 Stage Wrappers

Each stage wrapper constructs the subprocess command and executes it:

```python
class RFdiffusionRunner:
    def run(self, inputs: PreparedInputs, output_qv: Path) -> Path:
        """Run RFdiffusion backbone generation.

        Constructs command:
          rfdiffusion -t {target} -f {framework} -q {output}
                      -n {num_designs} -l "{cdr_loops}" -h "{hotspots}"

        Or Hydra-style if using scripts directly.
        """
        cmd = [
            "rfdiffusion",
            "-t", str(inputs.target_pdb),
            "-f", str(inputs.framework_hlt),
            "-q", str(output_qv),
            "-n", str(self.config.pipeline.rfdiffusion.num_designs),
            "-l", self._format_cdr_loops(),
            "-h", self._format_hotspots(inputs.hotspot_string),
        ]
        self._execute(cmd)
        return output_qv
```

### 6.3 Multi-GPU Parallelization

For campaigns > 1 GPU, use `qvsplit` to distribute work:

```python
class ParallelPipelineRunner:
    def run_stage_parallel(self, input_qv: Path, stage_fn, num_gpus: int) -> Path:
        """Split Quiver file across GPUs, run in parallel, merge results.

        1. qvsplit input.qv into N chunks (one per GPU)
        2. Launch N subprocesses with CUDA_VISIBLE_DEVICES set
        3. Wait for all to complete
        4. Merge output Quiver files
        """
        chunk_size = self._count_designs(input_qv) // num_gpus
        chunks = self._split_quiver(input_qv, chunk_size)

        processes = []
        for gpu_id, chunk in enumerate(chunks):
            env = {"CUDA_VISIBLE_DEVICES": str(gpu_id)}
            output = chunk.with_suffix(".out.qv")
            p = stage_fn(chunk, output, env=env)
            processes.append(p)

        for p in processes:
            p.wait()

        return self._merge_quivers([p.output for p in processes])
```

### 6.4 Checkpoint/Resume

Simple file-based checkpointing:

```python
def _checkpoint_exists(self, stage: str) -> bool:
    return (self.work_dir / f".checkpoint_{stage}").exists()

def _save_checkpoint(self, stage: str):
    (self.work_dir / f".checkpoint_{stage}").write_text(
        f"completed:{datetime.now().isoformat()}"
    )
```

Resume simply skips stages that have checkpoint files.

### 6.5 Error Handling

- Each subprocess call captures stdout/stderr
- Non-zero exit codes raise `StageFailedError` with full log
- GPU OOM errors detected via `pynvml` or subprocess stderr parsing
- On OOM: suggest reducing batch size or truncating target further
- All logs written to `{output_dir}/logs/{stage}_{timestamp}.log`

---

## 7. Results Analysis Module

### 7.1 Score Extraction

```python
def extract_scores(predictions_qv: Path) -> pd.DataFrame:
    """Extract all QV_SCORE values from predictions Quiver file.

    Uses: qvscorefile predictions.qv > scores.tsv

    Returns DataFrame with columns:
    - tag: design identifier
    - pae: predicted aligned error
    - rmsd: design vs prediction RMSD
    - ddg: Rosetta binding energy (if available)
    """
```

### 7.2 Filtering

```python
def apply_filters(
    scores: pd.DataFrame,
    pae_threshold: float = 10.0,
    rmsd_threshold: float = 2.0,
    ddg_threshold: float | None = -20.0
) -> pd.DataFrame:
    """Apply filtering thresholds to score table.

    Designs must pass ALL enabled filters.
    Returns filtered DataFrame with added 'pass' column.
    """
```

### 7.3 Ranking

```python
def rank_candidates(filtered: pd.DataFrame) -> pd.DataFrame:
    """Rank filtered candidates by composite score.

    Composite = weighted combination:
      0.4 * normalize(pAE, lower=better)
    + 0.3 * normalize(RMSD, lower=better)
    + 0.3 * normalize(ddG, lower=better)

    Returns DataFrame sorted by composite_score ascending.
    """
```

### 7.4 Report Generation

```python
def generate_report(
    ranked: pd.DataFrame,
    config: CampaignConfig,
    output_dir: Path
) -> Path:
    """Generate campaign summary report.

    Contents:
    - Campaign parameters summary
    - Design count at each stage
    - Filter pass rates (what % passed pAE, RMSD, ddG)
    - Score distributions (histograms)
    - Top 50 candidates table
    - PyMOL visualization commands for top 10
    - Recommendations for experimental testing
    """
```

### 7.5 Candidate Export

```python
def export_candidates(
    predictions_qv: Path,
    ranked: pd.DataFrame,
    top_n: int,
    output_dir: Path
) -> Path:
    """Export top N candidates as individual PDB files.

    Uses: qvextractspecific predictions.qv tag1 tag2 ... > candidate.pdb
    Also generates:
    - FASTA file with CDR sequences
    - PyMOL session file (.pse) for visual inspection
    - Summary CSV with scores
    """
```

---

## 8. Paper Reproduction Campaigns

### 8.1 Influenza Hemagglutinin (HA) - VHH

**Paper reference**: VHH_flu_01, cryo-EM 3.0 Å, RMSD 1.45 Å, ~66% HA particles bound

```yaml
# campaigns/paper_repro/ha_vhh.yaml
campaign:
  name: "ha_vhh_paper_repro"
  description: "Reproduce HA stem VHH designs from Bennett et al. 2025"

target:
  pdb_id: "4HMG"                      # HA from A/USA/Iowa/1943 H1N1
  chain_id: "A"
  epitope_residues: [18, 19, 20, 21, 38, 39, 40, 41, 42, 45, 46, 49, 52, 53, 56]  # HA stem conserved region
  hotspot_residues: [19, 40, 42, 45, 53]  # Hydrophobic stem contacts
  truncation:
    enabled: true
    buffer_angstroms: 10.0

antibody:
  format: "vhh"
  framework: "builtin:NbBCII10"       # NBBcH10F-GLA7 (paper framework)
  cdr_loops:
    H1: "7"
    H2: "6"
    H3: "5-13"

pipeline:
  rfdiffusion:
    num_designs: 9000                  # Paper: ~9000 VHH designs per target
  proteinmpnn:
    sequences_per_backbone: 5
    temperature: 0.2
  rf2:
    recycling_iterations: 10

filtering:
  pae_threshold: 10.0
  rmsd_threshold: 2.0

# Expected results for validation:
# - Multiple binders at 78 nM - 1.4 μM
# - Top design RMSD ~1.45 Å vs cryo-EM
# - ~66% HA particles bound in cryo-EM
```

### 8.2 C. difficile TcdB (Frizzled site) - VHH

**Paper reference**: VHH_TcdB_H2, 262 nM initial, EC50 neutralization 460 nM

```yaml
# campaigns/paper_repro/tcdb_vhh.yaml
campaign:
  name: "tcdb_vhh_paper_repro"
  description: "Reproduce TcdB Frizzled-site VHH from Bennett et al. 2025"

target:
  pdb_id: "6C0B"                      # TcdB DRBD domain (Frizzled-binding region)
  chain_id: "A"
  epitope_residues: [1285, 1286, 1287, 1288, 1340, 1341, 1342, 1343, 1344, 1400, 1401, 1402]  # Frizzled-binding interface
  hotspot_residues: [1286, 1340, 1342, 1401]
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
    num_designs: 9000
  proteinmpnn:
    sequences_per_backbone: 5
    temperature: 0.2
  rf2:
    recycling_iterations: 10

# Expected results:
# - Binders at 5.5 μM (RBD) to 262 nM (H2)
# - Neutralization EC50 = 460 nM (TcdB_H2)
# - Cryo-EM resolution 3.0-4.6 Å
```

### 8.3 C. difficile TcdB - scFv (Combinatorial)

**Paper reference**: scFv6, 72 nM, cryo-EM 3.6 Å, CDR RMSD 0.2-1.1 Å

```yaml
# campaigns/paper_repro/tcdb_scfv.yaml
campaign:
  name: "tcdb_scfv_paper_repro"
  description: "Reproduce TcdB scFv combinatorial designs from Bennett et al. 2025"

target:
  pdb_id: "6C0B"
  chain_id: "A"
  epitope_residues: [1285, 1286, 1287, 1288, 1340, 1341, 1342, 1343, 1344, 1400, 1401, 1402]
  hotspot_residues: [1286, 1340, 1342, 1401]
  truncation:
    enabled: true
    buffer_angstroms: 10.0

antibody:
  format: "scfv"
  framework: "builtin:hu4D5-8"
  cdr_loops:
    H1: "7"
    H2: "6"
    H3: "5-13"
    L1: "8-13"
    L2: "7"
    L3: "9-11"

pipeline:
  rfdiffusion:
    num_designs: 9000                  # Paper used ~10M combinatorial but from multiple runs
  proteinmpnn:
    sequences_per_backbone: 5
    temperature: 0.2
  rf2:
    recycling_iterations: 10

# Expected results:
# - 6 scFv binders identified
# - Best: scFv6, Kd = 72 nM
# - CDR RMSD 0.2-1.1 Å per loop
# - Overall RMSD 0.9 Å (full scFv)
```

### 8.4 RSV Site III - VHH

**Paper reference**: Multiple VHH binders at 78 nM and 1.4 μM

```yaml
# campaigns/paper_repro/rsv_vhh.yaml
campaign:
  name: "rsv_vhh_paper_repro"
  description: "Reproduce RSV site III VHH from Bennett et al. 2025"

target:
  pdb_id: "4JHW"                      # RSV G protein
  chain_id: "A"
  epitope_residues: []                # NOTE: needs literature lookup for RSV site III residues
  hotspot_residues: []                # TODO: identify from known neutralizing antibody epitopes
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
    num_designs: 9000

# Expected results:
# - Multiple binders at 78 nM - 1.4 μM
# - Competition with known site III antibodies confirmed
```

### 8.5 PHOX2B-HLA-C*07:02 - scFv (Neuroblastoma)

**Paper reference**: 6 binders, best scFv6 Kd = 400 nM (monomer), 68 nM (tetramer)

```yaml
# campaigns/paper_repro/phox2b_scfv.yaml
campaign:
  name: "phox2b_scfv_paper_repro"
  description: "Reproduce PHOX2B-HLA scFv for neuroblastoma from Bennett et al. 2025"

target:
  pdb_id: "TBD"                       # PHOX2B peptide-HLA-C*07:02 complex
  chain_id: "A"
  epitope_residues: []                # R6 residue critical
  hotspot_residues: []                # R6 + surrounding MHC residues
  truncation:
    enabled: false                    # pMHC is already small

antibody:
  format: "scfv"
  framework: "builtin:hu4D5-8"
  cdr_loops:
    H1: "7"
    H2: "6"
    H3: "5-13"
    L1: "8-13"
    L2: "7"
    L3: "9-11"

pipeline:
  rfdiffusion:
    num_designs: 9000

# Expected results:
# - 6 binders, best Kd = 400 nM monomer
# - Tetramer Kd = 68 nM
# - Specific for PHOX2B, no cross-reactivity with CSPG4
# CHALLENGE: Peptide-MHC = flat surface, low stability at R6
```

### 8.6 SARS-CoV-2 RBD - VHH

**Paper reference**: Design partial failure - correct fold/epitope but incorrect binding orientation

```yaml
# campaigns/paper_repro/sars_cov2_vhh.yaml
campaign:
  name: "sars_cov2_vhh_paper_repro"
  description: "Reproduce SARS-CoV-2 RBD VHH (known challenge) from Bennett et al. 2025"

target:
  pdb_id: "6M0J"                      # SARS-CoV-2 RBD
  chain_id: "E"
  epitope_residues: [417, 446, 449, 453, 455, 456, 475, 486, 487, 489, 493, 496, 498, 500, 501, 502, 505]
  hotspot_residues: [455, 486, 489, 493, 501]
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
    num_designs: 9000

# Expected results:
# - Modest affinity binders
# - KNOWN CHALLENGE: Framework-mediated binding instead of CDR-mediated
# - Use as negative control / benchmark for harness
```

---

## 9. New Cancer Target Campaigns

### 9.1 Priority Tier 1: Excellent Structural Data

#### 9.1.1 PD-L1 (CD274) - Immune Checkpoint

**Rationale**: Blocks PD-1/PD-L1 interaction, reactivates exhausted T cells. Approved: atezolizumab, durvalumab, avelumab.

```yaml
# campaigns/cancer/pdl1_vhh.yaml
campaign:
  name: "pdl1_vhh_v1"
  description: "De novo VHH targeting PD-L1 BC/DE loop region"

target:
  pdb_id: "5N2C"                      # PD-1/PD-L1 complex (2.3 Å)
  # Alternative: "6C0A" (atezolizumab Fab/PD-L1 complex)
  chain_id: "B"                       # PD-L1 chain
  epitope_residues: [54, 56, 58, 60, 62, 66, 68, 113, 115, 117, 121, 123, 124]
  hotspot_residues: [56, 115, 121]    # Y56, M115, I121 - hydrophobic core
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

# CHALLENGE: High surface flexibility, epitope accessibility varies with glycosylation
# STRATEGY: Target BC/DE loop interface (atezolizumab-like epitope)
```

#### 9.1.2 HER2 (ERBB2) - Breast/Gastric Cancer

**Rationale**: Growth factor receptor overexpressed in ~20% breast cancers. Approved: trastuzumab, pertuzumab.

```yaml
# campaigns/cancer/her2_vhh.yaml
campaign:
  name: "her2_vhh_v1"
  description: "De novo VHH targeting HER2 Domain IV (trastuzumab epitope region)"

target:
  pdb_id: "3N85"                      # HER2 ECD / trastuzumab Fab complex (1.55 Å)
  chain_id: "A"                       # HER2 chain
  epitope_residues: [557, 558, 560, 561, 567, 570, 571, 573, 574, 576, 579, 580, 583, 584, 586, 592, 596, 603]
  hotspot_residues: [560, 573, 579, 592]  # Hydrophobic contacts in Domain IV
  truncation:
    enabled: true
    buffer_angstroms: 10.0

antibody:
  format: "vhh"
  framework: "builtin:NbBCII10"
  cdr_loops:
    H1: "7"
    H2: "6"
    H3: "8-13"                        # Longer H3 for deep pocket engagement

pipeline:
  rfdiffusion:
    num_designs: 10000

# CHALLENGE: Strong steric constraints, must achieve picomolar affinity
# STRATEGY: Target Domain IV (trastuzumab epitope) - well-characterized, deep pocket
```

#### 9.1.3 EGFR - NSCLC / Colorectal Cancer

**Rationale**: Overexpressed in NSCLC and colorectal cancers. Approved: cetuximab, panitumumab.

```yaml
# campaigns/cancer/egfr_vhh.yaml
campaign:
  name: "egfr_vhh_v1"
  description: "De novo VHH targeting EGFR Domain III (cetuximab epitope)"

target:
  pdb_id: "5IUS"                      # EGFR ECD / panitumumab Fab complex
  # Alternative: "1YY9" (EGFR ECD / cetuximab Fab)
  chain_id: "A"
  epitope_residues: [384, 388, 390, 391, 465, 467, 470, 475, 476, 489]
  hotspot_residues: [467, 470, 489]   # L467, V470 hydrophobic + surrounding
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

# CHALLENGE: Multiple conformational states (active/inactive), DE loop flexibility
```

#### 9.1.4 CTLA-4 - Immune Checkpoint

**Rationale**: Blocks negative co-stimulatory signal. Approved: ipilimumab, tremelimumab.

```yaml
# campaigns/cancer/ctla4_vhh.yaml
campaign:
  name: "ctla4_vhh_v1"
  description: "De novo VHH targeting CTLA-4 FG loop region"

target:
  pdb_id: "4XIS"                      # CTLA-4 / ipilimumab Fab complex (2.6 Å)
  chain_id: "C"                       # CTLA-4 chain
  epitope_residues: [45, 47, 49, 50, 59, 62, 64, 68, 72, 74, 90, 122]
  hotspot_residues: [50, 59, 68, 72]  # Y50, W59, V68, F72 - aromatic/hydrophobic
  truncation:
    enabled: false                    # CTLA-4 ECD is already small (~15 kDa)

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

# CHALLENGE: CTLA-4 exists as homodimer; may need to bind dimer interface
```

#### 9.1.5 TROP-2 - Solid Tumors

**Rationale**: Overexpressed in breast, lung, colorectal cancers. Approved: sacituzumab govitecan (ADC).

```yaml
# campaigns/cancer/trop2_vhh.yaml
campaign:
  name: "trop2_vhh_v1"
  description: "De novo VHH targeting TROP-2 ectodomain"

target:
  pdb_id: "6YD6"                      # TROP-2 ectodomain (2.0 Å)
  # Alternative: "7JFF" (TROP-2 / sacituzumab Fab complex)
  chain_id: "A"
  epitope_residues: [62, 63, 68, 70, 71, 72, 73, 75, 76, 77, 79]
  hotspot_residues: [62, 68, 79]      # L62, I68, L79 - hydrophobic core
  truncation:
    enabled: false                    # TROP-2 ECD is manageable size

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

# CHALLENGE: Must achieve picomolar affinity for ADC applications
# OPPORTUNITY: VHH format enables novel epitopes inaccessible to conventional mAbs
```

### 9.2 Priority Tier 2: Moderate Structural Data

#### 9.2.1 CD20 - B-cell Lymphoma

```yaml
# campaigns/cancer/cd20_vhh.yaml
target:
  pdb_id: "5Y5M"                      # CD20/rituximab complex
  chain_id: "A"
  epitope_residues: [160, 161, 162, 167, 168, 169, 170, 172, 174, 175, 176, 177]
  hotspot_residues: [167, 168, 172]   # L167, I168, L172 - hydrophobic
  # CHALLENGE: Small accessible epitope (transmembrane protein)
  # CHALLENGE: Requires CDC for efficacy - VHH may lack Fc effector functions
```

#### 9.2.2 GPC3 (Glypican-3) - Hepatocellular Carcinoma

```yaml
# campaigns/cancer/gpc3_vhh.yaml
target:
  pdb_id: "6V3G"                      # GPC3 with antibody Fab
  chain_id: "A"
  epitope_residues: [45, 47, 52, 54, 65, 68, 71, 76, 79]  # Predicted from homology
  hotspot_residues: [45, 52, 71, 79]
  # CHALLENGE: Extensive heparan sulfate modifications mask epitopes
  # CHALLENGE: Limited high-resolution structural data
```

### 9.3 Priority Tier 3: Limited Structural Data (Exploratory)

#### 9.3.1 TIGIT - Emerging Immune Checkpoint

```yaml
# campaigns/cancer/tigit_vhh.yaml
target:
  pdb_id: "6VIC"                      # TIGIT ectodomain
  chain_id: "A"
  epitope_residues: [49, 56, 85, 93, 97, 116, 119]  # CC'/FG loops
  hotspot_residues: [49, 85, 93]      # Y49, W85, I93
  # CHALLENGE: Incomplete structural data, epitope overlap with PVRL2
```

#### 9.3.2 Claudin-18.2 - Gastric/Pancreatic Cancer

```yaml
# campaigns/cancer/cldn18_vhh.yaml
target:
  # NO DIRECT CRYSTAL STRUCTURE - requires homology model
  pdb_file: "/path/to/cldn18_homology_model.pdb"  # Built from claudin-15 (6OV7)
  chain_id: "A"
  epitope_residues: [108, 109, 110, 111, 113, 121, 124, 130, 137]  # Predicted
  hotspot_residues: [108, 111, 130]
  # CRITICAL WARNING: Homology model uncertainty is HIGH
  # RECOMMENDED: Wait for experimental structure or use AF2/AF3 prediction
```

#### 9.3.3 CD19 - B-cell Malignancies

```yaml
# campaigns/cancer/cd19_scfv.yaml
target:
  pdb_id: "6WJ7"                      # CD19 / inebilizumab Fab
  chain_id: "A"
  # NOTE: scFv format critical for CD19 (CAR-T applications)
antibody:
  format: "scfv"
  # CHALLENGE: Limited structural data, scFv stability requirements for CAR-T
```

### 9.4 Cancer Target Summary

| Target | PDB | Format | Designs | Priority | Key Challenge |
|--------|-----|--------|---------|----------|---------------|
| PD-L1 | 5N2C | VHH | 10,000 | Tier 1 | Flexibility, glycosylation |
| HER2 | 3N85 | VHH | 10,000 | Tier 1 | Affinity bar (picomolar) |
| EGFR | 5IUS | VHH | 10,000 | Tier 1 | Conformational states |
| CTLA-4 | 4XIS | VHH | 10,000 | Tier 1 | Homodimer |
| TROP-2 | 6YD6 | VHH | 10,000 | Tier 1 | ADC-grade affinity needed |
| CD20 | 5Y5M | VHH | 10,000 | Tier 2 | Small epitope, needs CDC |
| GPC3 | 6V3G | VHH | 10,000 | Tier 2 | Glycan masking |
| TIGIT | 6VIC | VHH | 10,000 | Tier 3 | Limited structures |
| CLDN18.2 | homology | VHH | 10,000 | Tier 3 | No crystal structure |
| CD19 | 6WJ7 | scFv | 10,000 | Tier 3 | CAR-T constraints |

---

## 10. Rare Disease Target Campaigns

### 10.1 GNE Myopathy - Enzyme Stabilization (Unconventional)

**Disease**: GNE myopathy (hereditary inclusion body myopathy, HIBM). Autosomal recessive, caused by loss-of-function mutations in GNE (UDP-N-acetylglucosamine 2-epimerase / N-acetylmannosamine kinase). Impairs sialic acid biosynthesis, leading to progressive muscle weakness.

**Therapeutic Hypothesis**: Design a VHH that binds the GNE domain-domain interface and stabilizes the active conformation of the mutant enzyme, partially restoring catalytic activity. This is an unconventional antibody application - stabilization rather than blockade.

**Caveats**:
- No precedent for antibody-mediated enzyme stabilization in clinical practice
- Intracellular target - antibody delivery to cytoplasm is a major challenge
- May require intrabody format or cell-penetrating peptide fusion
- Alternative: target extracellular sialylation-related pathways

```yaml
# campaigns/rare_disease/gne_stabilizer_vhh.yaml
campaign:
  name: "gne_stabilizer_vhh_v1"
  description: "VHH stabilizer for GNE mutant enzyme (GNE myopathy) - EXPLORATORY"

target:
  pdb_id: "1FOZ"                      # GNE full-length structure (2.8 Å)
  # Alternative: "3A3F" (catalytic domain), "2QM8" (with NAG substrate)
  chain_id: "A"
  epitope_residues: [10, 14, 41, 92, 94, 95, 98, 137, 140, 143]  # Domain interface
  hotspot_residues: [95, 98, 140, 143]  # L95, V98, L140, I143 - hydrophobic core
  truncation:
    enabled: false                    # GNE is not oversized

antibody:
  format: "vhh"
  framework: "builtin:NbBCII10"
  cdr_loops:
    H1: "7"
    H2: "6"
    H3: "8-16"                        # Extended H3 to span domain interface

pipeline:
  rfdiffusion:
    num_designs: 10000
  proteinmpnn:
    sequences_per_backbone: 10        # More sequence diversity for novel application
    temperature: 0.3

# CRITICAL NOTES:
# - Must NOT occlude NAG binding pocket (substrate accessibility)
# - Must accommodate known GNE mutations: L127P, D207N, M712V
# - Enzyme must remain catalytically active when bound
# - Design must be validated by enzyme activity assay, not just binding
```

### 10.2 PCSK9 - Familial Hypercholesterolemia (Conventional Alternative)

**Disease**: Familial hypercholesterolemia (FH). PCSK9 gain-of-function mutations cause extremely high LDL cholesterol. Also applicable to broader cardiovascular disease.

**Therapeutic Hypothesis**: Design VHH that blocks PCSK9-LDLR interaction, preventing LDL receptor degradation. Well-validated target with approved therapeutics (evolocumab, alirocumab).

```yaml
# campaigns/rare_disease/pcsk9_vhh.yaml
campaign:
  name: "pcsk9_vhh_v1"
  description: "De novo VHH targeting PCSK9 LDLR-binding site (FH)"

target:
  pdb_id: "3BPS"                      # PCSK9 / evolocumab Fab complex
  # Alternative: "3H42" (PCSK9 / LDLR-EGFa complex)
  chain_id: "A"                       # PCSK9 chain
  epitope_residues: [153, 155, 157, 159, 194, 197, 206, 207, 208, 238, 343, 366, 367, 369, 374, 375, 378, 380]
  hotspot_residues: [155, 207, 369, 374]  # Hydrophobic core at LDLR interface
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

# Well-validated target - good benchmark for rare disease campaign
# Approved therapeutics provide affinity/specificity benchmarks
# Expected Kd < 100 nM for therapeutic relevance
```

### 10.3 Complement C5 - PNH / aHUS

**Disease**: Paroxysmal nocturnal hemoglobinuria (PNH), atypical hemolytic uremic syndrome (aHUS). Terminal complement pathway hyperactivation.

**Therapeutic Hypothesis**: Block C5 cleavage to prevent C5a generation and MAC formation. Approved: eculizumab.

```yaml
# campaigns/rare_disease/complement_c5_vhh.yaml
campaign:
  name: "c5_vhh_v1"
  description: "De novo VHH targeting complement C5 cleavage site"

target:
  pdb_id: "3CU7"                      # C5 structure
  chain_id: "A"
  epitope_residues: [91, 94, 95, 97, 99, 102, 150, 156, 160, 163, 164, 168, 170]
  hotspot_residues: [94, 97, 160, 163]  # L94, V97, I160, L163
  truncation:
    enabled: true                     # C5 is large (~190 kDa)
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

# CHALLENGE: Very large target, conformational dynamics
# VHH advantage: smaller size may access cleavage-site pocket better than mAb
```

### 10.4 IL-6R - Autoimmune / Inflammatory

**Disease**: Rare autoimmune conditions (multicentric Castleman disease, systemic sclerosis). Also RA.

```yaml
# campaigns/rare_disease/il6r_vhh.yaml
campaign:
  name: "il6r_vhh_v1"
  description: "De novo VHH targeting IL-6R (autoimmune)"

target:
  pdb_id: "6S36"                      # IL-6R / tocilizumab Fab complex (cryo-EM 3.2 Å)
  # Alternative: "1ILR" (IL-6/IL-6R complex)
  chain_id: "A"                       # IL-6R chain
  epitope_residues: [96, 98, 99, 101, 103, 105, 174, 176, 177, 179, 180, 183, 185]
  hotspot_residues: [99, 103, 177, 180]  # L99, V103, L177, I180
  truncation:
    enabled: false

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

# CHALLENGE: Must block IL-6 binding without affecting trans-signaling selectivity
# Well-validated target (tocilizumab, sarilumab) provides benchmarks
```

---

## 11. Experimental Validation Pipeline

### 11.1 Overview

The computational pipeline produces ranked antibody candidates. This section specifies the experimental validation workflow to bridge from in silico designs to validated binders.

```
Computational Pipeline Output (Top 50-100 candidates)
                    |
                    v
    ┌───────────────────────────────┐
    │  Gene Synthesis & Cloning     │
    │  CDR sequences → DNA          │
    │  Clone into display vector    │
    └───────────────┬───────────────┘
                    |
                    v
    ┌───────────────────────────────┐
    │  Yeast Surface Display (YSD)  │
    │  High-throughput screening    │
    │  Flow cytometry sorting       │
    └───────────────┬───────────────┘
                    |
                    v
    ┌───────────────────────────────┐
    │  SPR Binding Characterization │
    │  Kd, kon, koff measurement    │
    │  Dose-response curves         │
    └───────────────┬───────────────┘
                    |
                    v
    ┌───────────────────────────────┐
    │  Functional Assays            │
    │  Competition, neutralization  │
    │  Cell-based binding           │
    └───────────────┬───────────────┘
                    |
                    v (optional)
    ┌───────────────────────────────┐
    │  OrthoRep Affinity Maturation │
    │  In vivo hypermutation        │
    │  2-3 log affinity improvement │
    └───────────────┬───────────────┘
                    |
                    v (for top hits)
    ┌───────────────────────────────┐
    │  Structural Validation        │
    │  Cryo-EM or X-ray             │
    │  Confirm binding mode         │
    └───────────────────────────────┘
```

### 11.2 Gene Synthesis & Cloning

**Input**: Top N candidate CDR sequences from computational pipeline

**Protocol**:
1. Extract CDR sequences from top candidates (FASTA export from harness)
2. Codon-optimize for yeast (S. cerevisiae) expression
3. Synthesize DNA oligos (Twist Bioscience, IDT, or similar)
4. Clone into yeast surface display vector:
   - VHH: pYDS2.0 or equivalent (Aga2p fusion)
   - scFv: VH-linker-VL format in same vector
5. Transform into EBY100 yeast strain
6. Verify by Sanger sequencing

**Scale**: 50-200 candidates per campaign (limited by synthesis cost)

**Harness integration**: `harness/experimental/synthesis_order.py` generates:
- FASTA file of CDR sequences
- Codon-optimized DNA sequences
- Synthesis order spreadsheet (compatible with vendor upload)
- Cloning strategy document

### 11.3 Yeast Surface Display (YSD) Screening

**Protocol** (from paper methods):

1. **Display induction**: Grow transformed yeast in galactose media (48h, 20°C)
2. **Target labeling**: Fluorescently label target protein (biotinylated + streptavidin-PE, or direct conjugation)
3. **Binding assay**:
   - Incubate yeast with labeled target at multiple concentrations
   - Paper used: 1.4 μM and 78 nM
   - Include anti-HA/anti-myc tag antibody for display detection
4. **Flow cytometry sorting**:
   - Gate: display-positive (tag+) AND binding-positive (target+)
   - Sort top 1-5% of binding population
   - 2-3 rounds of enrichment sorting
5. **Sequencing**: Individual clone sequencing of sorted populations

**Key parameters**:
| Parameter | Value | Notes |
|-----------|-------|-------|
| Target concentrations | 1.4 μM, 78 nM | Coarse then fine screening |
| Display detection | Anti-HA tag (FITC) | Reports surface display level |
| Binding detection | Target-PE or Target-biotin + SA-PE | Reports binding |
| Sort gates | Top 1-5% dual-positive | Balance specificity vs recovery |
| Sort rounds | 2-3 | Enrichment for binders |
| Clones sequenced | 48-96 per sort | Identify unique sequences |

**Harness integration**: `harness/experimental/yeast_display.py` generates:
- Experimental protocol document (step-by-step)
- Flow cytometry gating strategy
- Reagent list with quantities
- Expected timeline (2-3 weeks per round)

### 11.4 SPR Binding Characterization

**Protocol**:

1. **Immobilization**: Amine-couple target protein to CM5 sensor chip
2. **Analyte preparation**: Serial dilutions of purified antibody (0.1 - 1000 nM)
3. **Kinetics measurement**: Single-cycle or multi-cycle kinetics
4. **Data fitting**: 1:1 Langmuir binding model
5. **Output**: Kd, kon, koff for each candidate

**Key parameters**:
| Parameter | Value |
|-----------|-------|
| Analyte concentrations | 0.1, 1, 10, 100, 1000 nM |
| Association time | 120-300 seconds |
| Dissociation time | 300-600 seconds |
| Regeneration | 10 mM glycine pH 2.0 |
| Temperature | 25°C |

**Harness integration**: `harness/experimental/spr_protocol.py` generates:
- SPR protocol with target-specific parameters
- Concentration series calculator
- Expected Kd ranges based on computational predictions
- Data analysis template

### 11.5 Functional Assays

**Competition assay** (confirm epitope targeting):
- Pre-incubate target with known antibody (e.g., trastuzumab for HER2)
- Add designed antibody
- Measure binding displacement
- Confirms designed antibody targets intended epitope

**Neutralization assay** (for relevant targets):
- TcdB: Toxin neutralization on cells (measure cell viability)
- Viral targets (HA, RSV): Pseudovirus neutralization
- Immune checkpoints (PD-L1, CTLA-4): Reporter cell assays

**Cell-based binding** (for membrane targets):
- Flow cytometry on target-expressing cell lines
- Dose-dependent binding (MFI quantification)
- Specificity panel (target-negative cells as control)

### 11.6 OrthoRep Affinity Maturation

For binders with modest initial affinity (> 100 nM), the paper describes OrthoRep for in vivo continuous evolution:

**Protocol** (from paper):
1. Clone designed antibody into OrthoRep vector (error-prone DNA polymerase)
2. Express on yeast surface with mutagenic replication
3. Serial passaging with increasing stringency (decreasing target concentration)
4. FACS sorting at each passage for high-affinity variants
5. Multiple rounds (~weeks)
6. Isolate and sequence improved variants

**Expected improvement**: ~2 orders of magnitude (100x)
- Example from paper: VHH_TcdB_H2: 262 nM → ~10-20 nM

**Key consideration**: Binding mode must be preserved (cryo-EM validates this)

### 11.7 Structural Validation (Cryo-EM)

For top validated binders, cryo-EM confirms:
- Binding mode matches computational prediction (RMSD < 2 Å)
- CDR loop conformations are accurate
- Epitope contacts match intended design
- No unexpected binding modes (framework-mediated, etc.)

**Resolution targets**: < 4 Å for overall structure, < 3 Å for interface

---

## 12. Verification & Validation

### 12.1 Harness Smoke Test

**Minimal end-to-end test** using the upstream example:

```bash
# Run with example inputs (RSV site III target, NbBCII10 framework)
rfab campaign run campaigns/smoke_test.yaml

# smoke_test.yaml uses:
# - Target: scripts/examples/example_inputs/rsv_site3.pdb
# - Framework: scripts/examples/example_inputs/h-NbBCII10.pdb
# - Only 20 designs (fast)
```

**Validation criteria**:
- All 3 pipeline stages complete without error
- Output Quiver file contains 20 designs
- Score extraction produces valid TSV
- At least 1 design passes pAE < 10 filter (not guaranteed but expected)
- Report generated successfully
- Runtime < 30 min on single GPU

### 12.2 Paper Reproduction Validation

For each paper target, compare harness output metrics against published values:

| Target | Metric | Paper Value | Our Threshold |
|--------|--------|-------------|---------------|
| HA VHH | Hit rate | 2-3 binders / 9000 | >= 1 design with pAE < 5 |
| HA VHH | CDR RMSD | 0.91-1.45 Å | Mean RMSD < 3 Å |
| TcdB VHH | Hit rate | Multiple binders | >= 1 design with pAE < 5 |
| TcdB scFv | CDR RMSD | 0.2-1.1 Å/loop | Mean CDR RMSD < 3 Å |
| RSV VHH | Hit rate | 2+ binders / 9000 | >= 1 design with pAE < 5 |
| SARS-CoV-2 | Binding mode | Framework-mediated (failure) | Document as baseline |

**Important caveat**: Exact reproduction requires identical random seeds and software versions. We validate that our pipeline produces designs in the same quality range, not identical designs.

### 12.3 New Target Validation

For new cancer/rare disease targets, since we lack experimental data:

1. **Computational sanity checks**:
   - Pilot run (20-100 designs) produces docked structures
   - Visual inspection in PyMOL: antibody contacts epitope, not floating
   - Hotspot contacts within 8 Å
   - CDR loop geometry looks reasonable
   - pAE distribution similar to paper targets

2. **Cross-validation with existing therapeutics**:
   - For targets with approved antibodies (PD-L1, HER2, etc.)
   - Compare designed epitope contacts vs known therapeutic antibody epitope
   - Should show overlap if targeting same site

3. **Negative control**:
   - Run SARS-CoV-2 RBD campaign (known to produce suboptimal orientation)
   - Verify we reproduce the framework-mediated binding issue
   - Confirms pipeline is not just producing random hits

---

## 13. Project Directory Structure

```
repro_rfantibody_for-cancer-targets/
├── README.md
├── docs_plans/
│   └── v0_harness+repor_detailed-specs.md    # This document
│
├── RFAntibody/                                # Git submodule or clone of upstream
│   ├── src/rfantibody/                        # Core RFdiffusion code
│   ├── scripts/                               # Pipeline scripts
│   ├── bin/                                   # Quiver tools
│   ├── weights/                               # Model weights (downloaded)
│   └── ...
│
├── harness/                                   # Our harness code
│   ├── __init__.py
│   ├── cli.py                                 # Entry point
│   ├── config/
│   │   ├── __init__.py
│   │   ├── schema.py                          # CampaignConfig dataclass/Pydantic
│   │   └── defaults.py
│   ├── target_prep/
│   │   ├── __init__.py
│   │   ├── fetch_pdb.py
│   │   ├── truncate.py
│   │   ├── convert_framework.py
│   │   └── validate.py
│   ├── pipeline/
│   │   ├── __init__.py
│   │   ├── orchestrator.py
│   │   ├── rfdiffusion.py
│   │   ├── proteinmpnn.py
│   │   ├── rf2.py
│   │   └── parallel.py
│   ├── analysis/
│   │   ├── __init__.py
│   │   ├── filter.py
│   │   ├── rank.py
│   │   ├── report.py
│   │   └── export.py
│   └── experimental/
│       ├── __init__.py
│       ├── synthesis_order.py
│       ├── yeast_display.py
│       ├── spr_protocol.py
│       └── maturation.py
│
├── campaigns/                                 # YAML campaign configs
│   ├── smoke_test.yaml
│   ├── paper_repro/
│   │   ├── ha_vhh.yaml
│   │   ├── tcdb_vhh.yaml
│   │   ├── tcdb_scfv.yaml
│   │   ├── rsv_vhh.yaml
│   │   ├── phox2b_scfv.yaml
│   │   └── sars_cov2_vhh.yaml
│   ├── cancer/
│   │   ├── pdl1_vhh.yaml
│   │   ├── her2_vhh.yaml
│   │   ├── egfr_vhh.yaml
│   │   ├── ctla4_vhh.yaml
│   │   ├── trop2_vhh.yaml
│   │   ├── cd20_vhh.yaml
│   │   ├── gpc3_vhh.yaml
│   │   ├── tigit_vhh.yaml
│   │   ├── cldn18_vhh.yaml
│   │   └── cd19_scfv.yaml
│   └── rare_disease/
│       ├── gne_stabilizer_vhh.yaml
│       ├── pcsk9_vhh.yaml
│       ├── complement_c5_vhh.yaml
│       └── il6r_vhh.yaml
│
├── results/                                   # Campaign output (gitignored)
│   └── {campaign_name}/
│       ├── pipeline/
│       │   ├── 01_backbones.qv
│       │   ├── 02_sequences.qv
│       │   ├── 03_predictions.qv
│       │   └── .checkpoint_*
│       ├── analysis/
│       │   ├── scores.tsv
│       │   ├── filtered.tsv
│       │   ├── ranked.tsv
│       │   └── report.html
│       ├── candidates/
│       │   ├── top_001.pdb
│       │   ├── ...
│       │   └── sequences.fasta
│       └── experimental/
│           ├── synthesis_order.csv
│           ├── ysd_protocol.md
│           └── spr_protocol.md
│
├── tests/
│   ├── test_config.py
│   ├── test_target_prep.py
│   ├── test_pipeline.py
│   └── test_analysis.py
│
├── papers/
│   └── Atomically accurate de novo design of antibodies with RFdiffusion_____s41586-025-09721-5.pdf
│
├── pyproject.toml
└── .gitignore
```

---

## 14. Implementation Roadmap

### Phase 1: Foundation (Week 1-2)
1. Clone RFAntibody repo, verify installation, run example workflows
2. Implement `config/schema.py` - CampaignConfig with validation
3. Implement `target_prep/` module - PDB fetch, truncation, framework conversion
4. Implement basic `pipeline/orchestrator.py` - sequential 3-stage execution
5. Smoke test with upstream examples

### Phase 2: Core Pipeline (Week 2-3)
6. Implement stage wrappers (`rfdiffusion.py`, `proteinmpnn.py`, `rf2.py`)
7. Implement Quiver I/O integration
8. Implement `analysis/` module - score extraction, filtering, ranking
9. Implement checkpoint/resume
10. End-to-end test with smoke_test.yaml

### Phase 3: Paper Reproduction (Week 3-4)
11. Create all paper_repro campaign YAML configs
12. Run HA VHH campaign (primary validation)
13. Run TcdB VHH campaign
14. Run RSV VHH campaign
15. Compare metrics against paper values
16. Run SARS-CoV-2 as negative control
17. Document reproduction results

### Phase 4: New Targets (Week 4-6)
18. Create cancer target campaign YAML configs (Tier 1 first)
19. Run PD-L1, HER2, EGFR campaigns
20. Run CTLA-4, TROP-2 campaigns
21. Create rare disease configs (PCSK9, C5, IL-6R, GNE)
22. Run rare disease campaigns
23. Analyze cross-target patterns

### Phase 5: Polish (Week 6-7)
24. Implement multi-GPU parallelization
25. Implement report generation (HTML)
26. Implement experimental planning modules
27. CLI interface (`rfab campaign run/status/report`)
28. Write tests
29. Documentation

### Phase 6: Experimental Bridge (Week 7+)
30. Generate synthesis orders for top candidates per target
31. Plan yeast display screening campaigns
32. SPR protocol generation
33. OrthoRep maturation planning for modest-affinity binders

---

## Appendix A: Paper Key Numbers Quick Reference

| Parameter | Value | Source |
|-----------|-------|--------|
| Designs per VHH campaign | ~9,000 | Paper methods |
| Designs for scFv library | ~10 million (combinatorial) | Paper methods |
| VHH hit rate | 0-2% experimental | Paper results |
| scFv AF3 pass rate | ~4% | Paper results |
| RF2 pAE threshold | < 10 | Paper + README |
| RF2 RMSD threshold | < 2 Å | Paper + README |
| RF2 pBind threshold | > 0.99 | Paper |
| AF3 iPTM threshold | > 0.6 | Paper |
| CDR RMSD accuracy | 0.2-1.1 Å/loop | Paper cryo-EM |
| Affinity range (initial) | 78 nM - 5.5 μM | Paper |
| Affinity range (matured) | ~10-20 nM | Paper (OrthoRep) |
| Cryo-EM resolution range | 3.0-5.7 Å | Paper |
| Neutralization EC50 | 460 nM (TcdB) | Paper |
| TM-score for scFv pairing | > 0.9-1.0 | Paper |
| Hotspot distance cutoff | Cβ < 8 Å to CDR | Paper methods |
| Target truncation buffer | ~10 Å | DeepWiki/README |
| CDR ranges: H3 (typical) | 5-13 residues | Paper/README |
| CDR ranges: H3 (extended) | up to 16 residues | DeepWiki |
| RF2 recycling iterations | 10 (default) | README |
| GPU VRAM requirement | 8-16 GB typical | DeepWiki |

## Appendix B: Framework Sequences

**VHH (NbBCII10 / NBBcH10F-GLA7)**:
- Single-domain, ~15 kDa
- 3 CDRs (H1, H2, H3)
- Well-characterized, stable, high expression
- Humanized scaffold (low immunogenicity)

**scFv (hu-4D5-8 Fv)**:
- VH-linker-VL format, ~27 kDa
- 6 CDRs (H1, H2, H3, L1, L2, L3)
- Trastuzumab-derived framework
- VH and VL paired via flexible linker

## Appendix C: Known Limitations

1. **Weak filtering**: RF2 filtering shows weak enrichment. AF3 may improve but is not yet integrated.
2. **Orientation prediction**: SARS-CoV-2 RBD showed framework-mediated binding (design failure mode).
3. **scFv complexity**: 6-CDR design has higher failure rate than VHH 3-CDR design.
4. **No structure for CLDN18.2**: Homology model introduces significant uncertainty.
5. **GNE stabilization**: No precedent for antibody-mediated enzyme stabilization.
6. **Experimental gap**: Pipeline produces candidates but cannot predict experimental success rate.
7. **Affinity ceiling**: Initial designs typically 78 nM - 5 μM; OrthoRep needed for tighter binding.
