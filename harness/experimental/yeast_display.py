"""Generate yeast surface display (YSD) screening protocols."""

from __future__ import annotations

import logging
from pathlib import Path

from harness.config.schema import CampaignConfig

logger = logging.getLogger(__name__)

PROTOCOL_TEMPLATE = """\
# Yeast Surface Display Screening Protocol
## Campaign: {campaign_name}
## Target: {target_desc}

### 1. Construct Preparation
- Clone CDR sequences into Aga2p-fusion display vector ({vector})
- Transform into EBY100 yeast strain via electroporation
- Verify by colony PCR and Sanger sequencing
- Expected library size: {n_constructs} individual clones

### 2. Display Induction
- Grow transformants in SD-CAA media (30C, overnight)
- Transfer to SG-CAA media (galactose induction)
- Induce at 20C for 48 hours
- Verify display by anti-HA tag staining (FITC)
  - Expected display rate: 30-70% of cells

### 3. Target Labeling
- Biotinylate target protein ({target_id}) using NHS-PEG4-biotin
- Verify labeling by SDS-PAGE gel shift
- Prepare streptavidin-PE detection reagent

### 4. Binding Assay — Round 1 (Coarse Screen)
- Target concentration: {conc_high} nM
- Incubate 10^7 yeast cells with labeled target (1h, 4C)
- Wash 3x with PBS-BSA (0.1%)
- Stain with anti-HA-FITC (display) + SA-PE (binding)
- Flow cytometry: gate display+ / binding+ (top 5%)
- Sort into fresh SD-CAA media
- Regrow sorted population (48h)

### 5. Binding Assay — Round 2 (Stringent Screen)
- Target concentration: {conc_low} nM
- Same protocol as Round 1
- Sort top 1-3% of dual-positive population
- Regrow and plate individual colonies

### 6. Clone Characterization
- Pick 48-96 individual colonies
- Sequence by Sanger (CDR region)
- Re-test binding by flow cytometry (individual clones)
- Identify unique binders vs duplicates

### 7. Timeline
| Step | Duration |
|------|----------|
| Cloning + transformation | 1 week |
| Display induction | 2 days |
| Round 1 sort | 1 day |
| Regrowth | 2 days |
| Round 2 sort | 1 day |
| Colony picking + sequencing | 1 week |
| **Total** | **~3 weeks** |

### 8. Reagents
| Reagent | Quantity | Vendor |
|---------|----------|--------|
| {vector} vector | 1 ug | Addgene or in-house |
| EBY100 competent yeast | 10 vials | ATCC |
| Anti-HA-FITC | 50 uL | BioLegend |
| Streptavidin-PE | 50 uL | BioLegend |
| NHS-PEG4-biotin | 1 mg | Thermo |
| Target protein ({target_id}) | 500 ug | In-house or vendor |
| SD-CAA media | 2 L | In-house |
| SG-CAA media | 2 L | In-house |

### 9. Controls
- Positive: Known binder (if available for target)
- Negative: Empty vector (no insert) — display+, binding-
- Display control: Anti-HA only — confirms Aga2p surface display
"""


def generate_ysd_protocol(config: CampaignConfig, n_constructs: int, output_dir: Path) -> Path:
    """Generate a YSD screening protocol document."""
    output_dir.mkdir(parents=True, exist_ok=True)

    target_id = config.target.target_id
    concs = config.experimental.target_concentrations
    conc_high = concs[0] if concs else 1400.0
    conc_low = concs[1] if len(concs) > 1 else 78.0

    vector = "pYDS2.0" if config.antibody.format == "vhh" else "pYDS-scFv"

    content = PROTOCOL_TEMPLATE.format(
        campaign_name=config.campaign.name,
        target_desc=config.campaign.description,
        target_id=target_id,
        vector=vector,
        n_constructs=n_constructs,
        conc_high=conc_high,
        conc_low=conc_low,
    )

    path = output_dir / "ysd_protocol.md"
    path.write_text(content)
    logger.info("YSD protocol → %s", path)
    return path
