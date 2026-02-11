"""OrthoRep-based in vivo affinity maturation planning."""

from __future__ import annotations

import logging
from pathlib import Path

from harness.config.schema import CampaignConfig

logger = logging.getLogger(__name__)

PROTOCOL_TEMPLATE = """\
# OrthoRep Affinity Maturation Plan
## Campaign: {campaign_name}
## Starting Candidates: {n_candidates}

### Overview
OrthoRep enables continuous in vivo mutagenesis of a gene of interest using
an error-prone DNA polymerase (TP-DNAP1) in yeast. Combined with yeast
surface display and FACS selection, it achieves ~2 orders of magnitude
affinity improvement.

### 1. Candidate Selection Criteria
Select candidates for maturation based on:
- KD > 100 nM (modest binders with room for improvement)
- Correct binding mode (confirmed by competition assay or computational RMSD)
- Good expression / display levels on yeast surface
- Recommended: 3-5 candidates per maturation campaign

### 2. OrthoRep Vector Construction
- Clone candidate antibody gene into OrthoRep p1 plasmid
- Transform into OrthoRep-competent yeast strain (with TP-DNAP1 expression)
- Verify mutagenesis rate: ~10^-5 mutations/bp/generation (100x host genome)
- Expected diversity: ~1 mutation per antibody gene per generation

### 3. Selection Rounds
| Round | Target Conc (nM) | Sort Gate | Expected Improvement |
|-------|-------------------|-----------|---------------------|
| 1 | {conc_high} | Top 5% | 2-5x |
| 2 | {conc_mid} | Top 3% | 5-10x |
| 3 | {conc_low} | Top 1% | 10-50x |
| 4 | {conc_very_low} | Top 0.5% | 50-100x |

### 4. Protocol Per Round
1. Grow mutagenic culture in SG-CAA (galactose) at 30C for 48h
   - ~20 generations → ~20 mutations sampled per gene
2. Induce display: transfer to SG-CAA at 20C for 24h
3. Stain: anti-HA-FITC (display) + target-PE (binding)
4. FACS sort: gate as specified above
5. Recover sorted cells in SD-CAA
6. Passage 2-3 days before next round
7. **Duration per round**: 5-7 days

### 5. Analysis
After final round:
- Plate individual colonies (96-well)
- Sequence CDR regions (identify convergent mutations)
- Measure KD by flow cytometry titration or SPR
- Compare to parental clone

### 6. Expected Outcomes
- Starting KD: ~100-500 nM (typical for RFAntibody designs)
- Final KD: ~1-10 nM (after 3-4 rounds)
- From paper: TcdB_H2 improved from 262 nM to estimated ~10-20 nM
- Key: binding MODE should be preserved (verify by competition assay)

### 7. Timeline
| Step | Duration |
|------|----------|
| Vector construction | 2 weeks |
| Round 1-4 selection | 4 weeks |
| Clone characterization | 1 week |
| SPR validation | 1 week |
| **Total** | **~8 weeks** |

### 8. Critical Considerations
- Monitor for framework mutations (may destabilize fold)
- Periodically check display levels (avoid loss-of-display mutants)
- Cryo-EM recommended for top matured binder (confirm binding mode preserved)
- Binding mode validation is ESSENTIAL — improved affinity via altered binding
  mode is a failure (as observed with SARS-CoV-2 RBD designs in paper)
"""


def generate_maturation_plan(
    config: CampaignConfig, n_candidates: int, output_dir: Path
) -> Path:
    """Generate an OrthoRep affinity maturation plan."""
    output_dir.mkdir(parents=True, exist_ok=True)

    concs = config.experimental.spr_concentrations
    conc_high = concs[-1] if concs else 1000.0
    conc_mid = concs[-2] if len(concs) >= 2 else 100.0
    conc_low = concs[-3] if len(concs) >= 3 else 10.0
    conc_very_low = concs[-4] if len(concs) >= 4 else 1.0

    content = PROTOCOL_TEMPLATE.format(
        campaign_name=config.campaign.name,
        n_candidates=n_candidates,
        conc_high=conc_high,
        conc_mid=conc_mid,
        conc_low=conc_low,
        conc_very_low=conc_very_low,
    )

    path = output_dir / "maturation_plan.md"
    path.write_text(content)
    logger.info("Maturation plan → %s", path)
    return path
