"""Generate SPR (Surface Plasmon Resonance) binding characterization protocols."""

from __future__ import annotations

import logging
from pathlib import Path

from harness.config.schema import CampaignConfig

logger = logging.getLogger(__name__)

PROTOCOL_TEMPLATE = """\
# SPR Binding Characterization Protocol
## Campaign: {campaign_name}
## Target: {target_id}

### 1. Chip Preparation
- Sensor chip: CM5 (GE Healthcare / Cytiva)
- Immobilization: Amine coupling (EDC/NHS activation)
- Ligand: {target_id} target protein
- Immobilization level: 200-500 RU
- Reference channel: Blank (activated/deactivated, no protein)

### 2. Analyte Preparation
- Purify candidate antibodies by Ni-NTA or Protein A chromatography
- Buffer exchange into HBS-EP+ (10 mM HEPES pH 7.4, 150 mM NaCl, 3 mM EDTA, 0.05% P20)
- Prepare serial dilutions: {concentrations}
- All concentrations in nM

### 3. Kinetics Measurement
- Method: Multi-cycle kinetics (MCK)
- Association time: 180 seconds
- Dissociation time: 300 seconds
- Flow rate: 30 uL/min
- Temperature: 25C
- Regeneration: 10 mM glycine pH 2.0, 30 seconds

### 4. Data Analysis
- Software: Biacore Insight Evaluation (or equivalent)
- Fitting model: 1:1 Langmuir binding
- Report: ka (on-rate), kd (off-rate), KD (equilibrium)
- Quality criteria:
  - Chi2 < 1% of Rmax
  - U-value < 25 (uniqueness)
  - Residuals: random distribution

### 5. Expected Results
| Parameter | Acceptable Range | Therapeutic Target |
|-----------|------------------|--------------------|
| KD | < 1 uM | < 10 nM |
| ka | > 10^4 M-1s-1 | > 10^5 M-1s-1 |
| kd | < 10^-2 s-1 | < 10^-4 s-1 |

### 6. Concentration Series
{conc_table}

### 7. Run Order
1. Buffer blanks (2x) — baseline stability
2. Highest concentration (single inject) — verify binding
3. Full concentration series (low → high)
4. Buffer blank — check regeneration
5. Repeat lowest concentration — reproducibility check

### 8. Troubleshooting
| Issue | Likely Cause | Solution |
|-------|-------------|----------|
| No binding | Inactive analyte | Re-purify, check by SEC |
| Mass transport | Too much ligand | Reduce immobilization level |
| Incomplete regeneration | Strong binding | Try pH 1.5 or add NaCl |
| Bulk shift | Buffer mismatch | Filter/degas running buffer |
"""


def generate_spr_protocol(config: CampaignConfig, output_dir: Path) -> Path:
    """Generate an SPR binding characterization protocol."""
    output_dir.mkdir(parents=True, exist_ok=True)

    target_id = config.target.pdb_id or Path(config.target.pdb_file).stem
    concentrations = config.experimental.spr_concentrations
    conc_str = ", ".join(f"{c}" for c in concentrations) + " nM"

    conc_table = "| # | Concentration (nM) |\n|---|--------------------|\n"
    for i, c in enumerate(concentrations, 1):
        conc_table += f"| {i} | {c} |\n"

    content = PROTOCOL_TEMPLATE.format(
        campaign_name=config.campaign.name,
        target_id=target_id,
        concentrations=conc_str,
        conc_table=conc_table,
    )

    path = output_dir / "spr_protocol.md"
    path.write_text(content)
    logger.info("SPR protocol → %s", path)
    return path
