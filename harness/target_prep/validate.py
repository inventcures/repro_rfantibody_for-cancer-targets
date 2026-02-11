"""Validate prepared inputs before handing off to the pipeline."""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path

from harness._bio_utils import THREE_TO_ONE, find_chain
from harness.config.defaults import HYDROPHOBIC_RESIDUES

logger = logging.getLogger(__name__)


@dataclass
class ValidationResult:
    valid: bool = True
    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)


def validate_hotspots(
    target_pdb: Path,
    hotspot_residues: list[int],
    chain_id: str,
) -> ValidationResult:
    """Check that the selected hotspot residues are sane."""
    from Bio.PDB import PDBParser

    result = ValidationResult()
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(target_pdb.stem, str(target_pdb))

    chain = find_chain(structure, chain_id)
    if chain is None:
        result.valid = False
        result.errors.append(f"Chain {chain_id} not found in {target_pdb}")
        return result

    residue_map = {res.id[1]: res for res in chain.get_residues()}

    missing = [r for r in hotspot_residues if r not in residue_map]
    if missing:
        result.valid = False
        result.errors.append(f"Hotspot residues not found in structure: {missing}")

    hydrophobic_count = 0
    for resid in hotspot_residues:
        res = residue_map.get(resid)
        if res is None:
            continue
        aa = THREE_TO_ONE.get(res.get_resname(), "")
        if aa in HYDROPHOBIC_RESIDUES:
            hydrophobic_count += 1

    if hydrophobic_count < 3:
        result.warnings.append(
            f"Only {hydrophobic_count} hydrophobic hotspot residues "
            f"(need >= 3 for stable binding)"
        )

    _check_patch_contiguity(residue_map, hotspot_residues, result)

    return result


def validate_framework_hlt(framework_pdb: Path, antibody_format: str) -> ValidationResult:
    """Validate that a framework PDB has the expected HLT chain labels."""
    from Bio.PDB import PDBParser

    result = ValidationResult()
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(framework_pdb.stem, str(framework_pdb))

    chain_ids = {chain.id for model in structure for chain in model}

    if "H" not in chain_ids:
        result.valid = False
        result.errors.append("Framework PDB missing heavy chain 'H'")

    if antibody_format == "scfv" and "L" not in chain_ids:
        result.valid = False
        result.errors.append("scFv framework PDB missing light chain 'L'")

    return result


def _check_patch_contiguity(
    residue_map: dict, hotspot_residues: list[int], result: ValidationResult
) -> None:
    """Warn if hotspot residues are far apart in sequence (likely scattered)."""
    if len(hotspot_residues) < 2:
        return

    sorted_ids = sorted(hotspot_residues)
    max_gap = max(
        sorted_ids[i + 1] - sorted_ids[i] for i in range(len(sorted_ids) - 1)
    )
    if max_gap > 30:
        result.warnings.append(
            f"Hotspot residues span a large sequence gap (max gap={max_gap}). "
            "Ensure they form a spatial patch on the target surface."
        )
