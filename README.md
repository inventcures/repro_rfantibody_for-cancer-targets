# rfab-harness

[![Ask DeepWiki](https://deepwiki.com/badge.svg)](https://deepwiki.com/inventcures/repro_rfantibody_for-cancer-targets)

Campaign orchestration harness for [RFAntibody](https://github.com/RosettaCommons/RFAntibody) — the de novo antibody design pipeline from Bennett et al., *"Atomically accurate de novo design of single-domain antibodies"* (Nature, 2025).

Wraps the 3-stage pipeline (RFdiffusion → ProteinMPNN → RF2) with campaign configuration, target preparation, automated analysis, and 21 pre-built configs for cancer and rare disease targets.

**Homepage**: [inventcures.github.io/harness_for_rfantibody](https://inventcures.github.io/harness_for_rfantibody)

## Architecture

```
Campaign YAML ──▶ Target Prep ──▶ RFdiffusion ──▶ ProteinMPNN ──▶ RF2 ──▶ Analysis
                  (fetch/trunc)   (backbones)     (sequences)    (predict) (filter/rank)
```

- **Target Prep** — PDB fetch from RCSB, epitope-based truncation, hotspot validation, framework conversion to HLT format
- **Pipeline** — Sequential 3-stage execution with Quiver I/O, multi-GPU parallelization via `qvsplit`, checkpoint/resume
- **Analysis** — Score extraction, filtering (pAE < 10, RMSD < 2 Å, ddG < -20 REU), composite ranking, HTML/CSV reports, PDB export

## Prerequisites

- Python 3.10+
- [RFAntibody](https://github.com/RosettaCommons/RFAntibody) cloned locally with weights
- CUDA GPU for pipeline execution

## Installation

```bash
git clone https://github.com/inventcures/repro_rfantibody_for-cancer-targets.git
cd repro_rfantibody_for-cancer-targets
pip install -e .
```

## Quick Start

```bash
# Validate config
rfab validate campaigns/smoke_test.yaml

# Dry run — prepare inputs only, no GPU needed
rfab run campaigns/smoke_test.yaml --dry-run --rfantibody-root ./RFAntibody

# Full campaign (20 designs, ~10 min on 1 GPU)
rfab run campaigns/smoke_test.yaml --rfantibody-root ./RFAntibody

# Re-run analysis on existing predictions
rfab analyze campaigns/smoke_test.yaml
```

## Campaign Config

```yaml
campaign:
  name: "my_campaign"

target:
  pdb_id: "5N2C"
  chain_id: "A"
  epitope_residues: [54, 56, 58, 60, 62, 66, 68, 113, 115, 117]
  hotspot_residues: [56, 60, 115]

antibody:
  format: "vhh"                    # vhh | scfv
  framework: "builtin:NbBCII10"    # or path to custom PDB
  cdr_loops:
    H1: "7"                        # fixed length
    H2: "6"
    H3: "5-13"                     # variable range

pipeline:
  rfdiffusion:
    num_designs: 10000
  proteinmpnn:
    sequences_per_backbone: 5
    temperature: 0.2

filtering:
  pae_threshold: 10.0
  rmsd_threshold: 2.0
  ddg_threshold: -20.0
```

## Pre-built Campaigns

### Paper Reproductions

| Config | Target | Format | PDB |
|--------|--------|--------|-----|
| `paper_repro/ha_vhh.yaml` | Influenza HA stem | VHH | 4BGW |
| `paper_repro/tcdb_vhh.yaml` | C. difficile TcdB | VHH | 7UMN |
| `paper_repro/tcdb_scfv.yaml` | C. difficile TcdB | scFv | 7UMN |
| `paper_repro/rsv_vhh.yaml` | RSV Site III | VHH | 4JHW |
| `paper_repro/phox2b_scfv.yaml` | PHOX2B-HLA neoantigen | scFv | modeled |
| `paper_repro/sars_cov2_vhh.yaml` | SARS-CoV-2 RBD | VHH | 6M0J |

### Cancer Targets

| Config | Target | Indication | PDB |
|--------|--------|------------|-----|
| `cancer/pdl1_vhh.yaml` | PD-L1 | Immune checkpoint | 5N2C |
| `cancer/her2_vhh.yaml` | HER2 | Breast/gastric | 1N8Z |
| `cancer/egfr_vhh.yaml` | EGFR | NSCLC/colorectal | 1NQL |
| `cancer/ctla4_vhh.yaml` | CTLA-4 | Immune checkpoint | 1I8L |
| `cancer/trop2_vhh.yaml` | TROP-2 | Solid tumors | 7E5M |
| `cancer/cd20_vhh.yaml` | CD20 | Lymphoma | 6Y4A |
| `cancer/gpc3_vhh.yaml` | GPC3 | Hepatocellular carcinoma | 7YIO |
| `cancer/tigit_vhh.yaml` | TIGIT | Emerging checkpoint | 6V33 |
| `cancer/cldn18_vhh.yaml` | Claudin-18.2 | Gastric/pancreatic | 7RFB |
| `cancer/cd19_scfv.yaml` | CD19 | B-cell malignancies | 6AL5 |

### Rare Disease Targets

| Config | Target | Indication | PDB |
|--------|--------|------------|-----|
| `rare_disease/complement_c5_vhh.yaml` | Complement C5 | PNH/aHUS | 3CU7 |
| `rare_disease/pcsk9_vhh.yaml` | PCSK9 | Familial hypercholesterolemia | 3BPS |
| `rare_disease/il6r_vhh.yaml` | IL-6R | Systemic JIA | 1N26 |
| `rare_disease/gne_stabilizer_vhh.yaml` | GNE (enzyme stabilizer) | GNE myopathy | 4WMN |

## Project Structure

```
harness/
├── cli.py                    # rfab entry point (run / validate / analyze)
├── config/
│   ├── schema.py             # CampaignConfig dataclasses + 15 validation rules
│   └── defaults.py           # builtin frameworks, CDR ranges, thresholds
├── target_prep/
│   ├── fetch_pdb.py          # RCSB download + chain extraction
│   ├── truncate.py           # epitope-based target truncation
│   ├── convert_framework.py  # Chothia → HLT format conversion
│   └── validate.py           # hotspot + framework validation
├── pipeline/
│   ├── orchestrator.py       # 3-stage execution with checkpointing
│   ├── rfdiffusion.py        # Stage 1 wrapper
│   ├── proteinmpnn.py        # Stage 2 wrapper
│   ├── rf2.py                # Stage 3 wrapper
│   └── parallel.py           # multi-GPU via qvsplit
├── analysis/
│   ├── filter.py             # score extraction + threshold filtering
│   ├── rank.py               # composite scoring (0.4 pAE + 0.3 RMSD + 0.3 ddG)
│   ├── report.py             # HTML/CSV report generation
│   └── export.py             # PDB + FASTA export
└── experimental/
    ├── synthesis_order.py    # gene synthesis orders (yeast codon-optimized)
    ├── yeast_display.py      # YSD screening protocol
    ├── spr_protocol.py       # SPR binding characterization
    └── maturation.py         # OrthoRep affinity maturation plan
```

## Testing

```bash
python3 -m unittest discover -s tests -v
```

33 tests. No GPU or RFAntibody installation required. Analysis tests require pandas (gracefully skipped otherwise).

## CLI Reference

| Command | Description |
|---------|-------------|
| `rfab run <config> [--dry-run] [-v]` | Run full design campaign |
| `rfab validate <config>` | Validate campaign config |
| `rfab analyze <config> [-v]` | Re-run analysis on existing predictions |

Global option: `--rfantibody-root PATH` (default: `./RFAntibody`)

## Citation

```bibtex
@article{bennett2025rfantibody,
  title={Atomically accurate de novo design of single-domain antibodies},
  author={Bennett, Nathaniel R and Watson, Joseph L and Ragotte, Robert J and others},
  journal={Nature},
  year={2025},
  doi={10.1038/s41586-024-08480-5}
}
```

## License

This harness is provided under the MIT license. The upstream [RFAntibody pipeline](https://github.com/RosettaCommons/RFAntibody) has its own license — see that repository for details.
