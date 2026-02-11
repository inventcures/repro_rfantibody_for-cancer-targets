"""Download PDB structures from RCSB."""

from __future__ import annotations

import logging
from pathlib import Path
from urllib.request import urlretrieve

from harness.config.defaults import RCSB_BASE_URL

logger = logging.getLogger(__name__)


def fetch_pdb(pdb_id: str, output_dir: Path, chain_id: str | None = None) -> Path:
    """Download a PDB file from RCSB and optionally extract a single chain.

    Returns path to the downloaded (or chain-extracted) PDB file.
    """
    pdb_id = pdb_id.strip().upper()
    output_dir.mkdir(parents=True, exist_ok=True)

    url = f"{RCSB_BASE_URL}/{pdb_id}.pdb"
    raw_path = output_dir / f"{pdb_id}.pdb"

    if raw_path.exists():
        logger.info("PDB %s already downloaded: %s", pdb_id, raw_path)
    else:
        logger.info("Downloading %s from RCSB ...", pdb_id)
        urlretrieve(url, raw_path)
        logger.info("Saved %s", raw_path)

    if chain_id is None:
        return raw_path

    return extract_chain(raw_path, chain_id, output_dir)


def extract_chain(pdb_path: Path, chain_id: str, output_dir: Path) -> Path:
    """Extract a single chain from a PDB file using BioPython."""
    from Bio.PDB import PDBParser, PDBIO, Select

    class ChainSelect(Select):
        def accept_chain(self, chain):
            return chain.id == chain_id

    stem = pdb_path.stem
    out_path = output_dir / f"{stem}_chain{chain_id}.pdb"
    if out_path.exists():
        logger.info("Chain-extracted PDB already exists: %s", out_path)
        return out_path

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(stem, str(pdb_path))

    io = PDBIO()
    io.set_structure(structure)
    io.save(str(out_path), ChainSelect())
    logger.info("Extracted chain %s → %s", chain_id, out_path)
    return out_path
