"""Download PDB structures from RCSB."""

from __future__ import annotations

import logging
from pathlib import Path
from urllib.error import HTTPError
from urllib.request import urlretrieve

from harness.config.defaults import RCSB_BASE_URL

logger = logging.getLogger(__name__)


def fetch_pdb(
    pdb_id: str,
    output_dir: Path,
    chain_id: str | None = None,
    source_chain_id: str | None = None,
) -> Path:
    """Download a PDB file from RCSB and optionally extract a single chain.

    Falls back to mmCIF format when PDB format is unavailable (e.g. 7U8C).
    When *source_chain_id* differs from *chain_id*, the source chain is
    extracted from the structure and renamed to *chain_id* in the output PDB.

    Returns path to the downloaded (or chain-extracted) PDB file.
    """
    pdb_id = pdb_id.strip().upper()
    output_dir.mkdir(parents=True, exist_ok=True)

    raw_path = _download(pdb_id, output_dir)

    if chain_id is None:
        return raw_path

    actual_source = source_chain_id or chain_id
    return extract_chain(raw_path, chain_id, output_dir, source_chain_id=actual_source)


def _download(pdb_id: str, output_dir: Path) -> Path:
    """Download PDB file, falling back to mmCIF if PDB format unavailable."""
    pdb_path = output_dir / f"{pdb_id}.pdb"
    cif_path = output_dir / f"{pdb_id}.cif"

    if pdb_path.exists():
        logger.info("PDB %s already downloaded: %s", pdb_id, pdb_path)
        return pdb_path
    if cif_path.exists():
        logger.info("CIF %s already downloaded: %s", pdb_id, cif_path)
        return cif_path

    pdb_url = f"{RCSB_BASE_URL}/{pdb_id}.pdb"
    logger.info("Downloading %s from RCSB ...", pdb_id)
    try:
        urlretrieve(pdb_url, pdb_path)
        logger.info("Saved %s", pdb_path)
        return pdb_path
    except HTTPError as e:
        if e.code == 404:
            logger.info("PDB format unavailable for %s, trying mmCIF...", pdb_id)
        else:
            raise
    except Exception:
        if pdb_path.exists():
            pdb_path.unlink()
        raise

    cif_url = f"{RCSB_BASE_URL}/{pdb_id}.cif"
    try:
        urlretrieve(cif_url, cif_path)
        logger.info("Saved CIF %s", cif_path)
        return cif_path
    except Exception:
        if cif_path.exists():
            cif_path.unlink()
        raise


def _parse_structure(path: Path):
    """Parse a PDB or CIF file into a BioPython Structure."""
    if path.suffix == ".cif":
        from Bio.PDB import MMCIFParser
        parser = MMCIFParser(QUIET=True)
    else:
        from Bio.PDB import PDBParser
        parser = PDBParser(QUIET=True, PERMISSIVE=True)
    return parser.get_structure(path.stem, str(path))


def extract_chain(
    pdb_path: Path,
    chain_id: str,
    output_dir: Path,
    source_chain_id: str | None = None,
) -> Path:
    """Extract a single chain from a PDB/CIF file using BioPython.

    When *source_chain_id* differs from *chain_id*, the source chain is
    extracted and renamed to *chain_id* in the output PDB.
    """
    from Bio.PDB import PDBIO, Select

    actual_source = source_chain_id or chain_id

    class ChainSelect(Select):
        def accept_chain(self, chain):
            return chain.id == actual_source

        def accept_residue(self, residue):
            return residue.id[0] == " "

    stem = pdb_path.stem
    out_path = output_dir / f"{stem}_chain{chain_id}.pdb"

    structure = _parse_structure(pdb_path)

    # Rename source chain to target chain_id before saving
    if actual_source != chain_id:
        for model in structure:
            for chain in model:
                if chain.id == actual_source:
                    chain.id = chain_id
                    logger.info("Renamed chain %s → %s", actual_source, chain_id)
                    break
        # Update selector to use new ID
        class ChainSelect(Select):
            def accept_chain(self, chain):
                return chain.id == chain_id

            def accept_residue(self, residue):
                return residue.id[0] == " "

    io = PDBIO()
    io.set_structure(structure)
    io.save(str(out_path), ChainSelect())
    logger.info("Extracted chain %s → %s", actual_source, out_path)
    return out_path
