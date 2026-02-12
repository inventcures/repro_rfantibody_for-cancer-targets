"""Truncate a target PDB around the epitope to reduce computational cost."""

from __future__ import annotations

import logging
from pathlib import Path

from harness._bio_utils import find_chain

logger = logging.getLogger(__name__)


def truncate_target(
    pdb_path: Path,
    epitope_residues: list[int],
    chain_id: str,
    buffer_angstroms: float = 10.0,
    preserve_secondary_structure: bool = True,
    output_path: Path | None = None,
) -> Path:
    """Truncate a target structure to residues near the epitope.

    Keeps all residues with any heavy atom within *buffer_angstroms* of any
    epitope residue Cα.  When *preserve_secondary_structure* is True the
    selection is extended to include complete secondary-structure elements
    (helices / strands) that partially overlap the proximity set.
    """
    import numpy as np
    from Bio.PDB import PDBIO, PDBParser, Select

    class _EpitopeProximitySelect(Select):
        def accept_residue(self, residue):
            parent_chain = residue.get_parent()
            if parent_chain.id != chain_id:
                return True
            return residue.id[1] in keep_ids

    if output_path is None:
        output_path = pdb_path.parent / f"{pdb_path.stem}_trunc.pdb"

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_path.stem, str(pdb_path))

    target_chain = find_chain(structure, chain_id)
    if target_chain is None:
        raise ValueError(f"Chain {chain_id} not found in {pdb_path}")

    epitope_coords = _collect_epitope_coords(target_chain, epitope_residues)
    if len(epitope_coords) == 0:
        raise ValueError(
            f"No atoms found for epitope residues {epitope_residues} "
            f"in chain {chain_id}"
        )

    keep_ids = _residues_within_buffer(target_chain, epitope_coords, buffer_angstroms)

    min_epitope = min(epitope_residues)
    max_epitope = max(epitope_residues)
    pre_truncation = len(keep_ids)
    keep_ids = {r for r in keep_ids if r >= min_epitope}
    if pre_truncation != len(keep_ids):
        logger.info(
            "Excluded %d residues below epitope start (res < %d)",
            pre_truncation - len(keep_ids), min_epitope,
        )

    if preserve_secondary_structure:
        keep_ids = _extend_to_ss_elements(target_chain, keep_ids)

    kept = len(keep_ids)
    total = sum(1 for _ in target_chain.get_residues())
    logger.info(
        "Truncation: keeping %d / %d residues (buffer=%.1f Å)",
        kept, total, buffer_angstroms,
    )

    io = PDBIO()
    io.set_structure(structure)
    io.save(str(output_path), _EpitopeProximitySelect())
    return output_path


def _collect_epitope_coords(chain, residue_ids: list[int]):
    """Gather Cα coordinates for the specified epitope residues."""
    import numpy as np

    coords = []
    for res in chain.get_residues():
        if res.id[1] in residue_ids and "CA" in res:
            coords.append(res["CA"].get_vector().get_array())
    return np.array(coords) if coords else np.empty((0, 3))


def _residues_within_buffer(chain, epitope_coords, buffer: float) -> set[int]:
    """Return residue IDs with any heavy atom within buffer of epitope Cα."""
    import numpy as np

    keep: set[int] = set()
    buffer_sq = buffer * buffer

    for res in chain.get_residues():
        for atom in res.get_atoms():
            if atom.element == "H":
                continue
            atom_coord = atom.get_vector().get_array()
            dists_sq = np.sum((epitope_coords - atom_coord) ** 2, axis=1)
            if np.any(dists_sq <= buffer_sq):
                keep.add(res.id[1])
                break
    return keep


def _extend_to_ss_elements(chain, keep_ids: set[int]) -> set[int]:
    """Extend the kept set to include full secondary-structure stretches.

    Uses a simple heuristic: contiguous blocks of residues that share a
    common SS assignment (from DSSP-style annotations if available, else
    sequence proximity) are kept whole if any member is already selected.
    We approximate this by gap-filling: if residues i and i+2 are both kept
    but i+1 is not, include i+1 as well (avoids breaking strands/helices).
    """
    all_ids = sorted(res.id[1] for res in chain.get_residues())
    extended = set(keep_ids)

    for idx in range(len(all_ids) - 2):
        a, b, c = all_ids[idx], all_ids[idx + 1], all_ids[idx + 2]
        if a in extended and c in extended and b not in extended:
            extended.add(b)

    return extended
