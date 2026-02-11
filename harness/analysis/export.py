"""Export top candidate designs as individual PDB files and FASTA."""

from __future__ import annotations

import logging
import subprocess
from pathlib import Path
from typing import TYPE_CHECKING

from harness._bio_utils import THREE_TO_ONE

if TYPE_CHECKING:
    import pandas as pd

logger = logging.getLogger(__name__)


def export_candidates(
    predictions_qv: Path,
    ranked: pd.DataFrame,
    top_n: int,
    output_dir: Path,
    rfantibody_root: Path,
) -> Path:
    """Export the top N candidates as individual PDB files and a FASTA summary."""
    output_dir.mkdir(parents=True, exist_ok=True)
    top = ranked.head(top_n)

    if "tag" not in top.columns:
        logger.warning("No 'tag' column in ranked data; cannot export PDBs")
        return output_dir

    tags = top["tag"].tolist()
    _extract_pdbs(predictions_qv, tags, output_dir, rfantibody_root)
    _write_summary_csv(top, output_dir)
    _write_fasta(output_dir, tags)

    logger.info("Exported %d candidates → %s", len(tags), output_dir)
    return output_dir


def _extract_pdbs(
    qv_path: Path, tags: list[str], output_dir: Path, rfantibody_root: Path
) -> None:
    """Extract individual PDBs from a Quiver file using qvextractspecific."""
    qvextract = rfantibody_root / "bin" / "qvextractspecific"

    if qvextract.exists():
        for tag in tags:
            out_pdb = output_dir / f"{tag}.pdb"
            result = subprocess.run(
                [str(qvextract), str(qv_path), tag],
                capture_output=True, text=True,
            )
            if result.returncode == 0:
                out_pdb.write_text(result.stdout)
            else:
                logger.warning("Failed to extract %s: %s", tag, result.stderr.strip())
    else:
        _extract_from_text(qv_path, tags, output_dir)


def _extract_from_text(qv_path: Path, tags: list[str], output_dir: Path) -> None:
    """Fallback: parse Quiver text to extract PDB blocks for specified tags."""
    tag_set = set(tags)
    current_tag = None
    current_lines: list[str] = []

    with open(qv_path) as f:
        for line in f:
            if line.startswith("QV_TAG"):
                if current_tag in tag_set and current_lines:
                    (output_dir / f"{current_tag}.pdb").write_text("".join(current_lines))
                current_tag = line.split(None, 1)[1].strip() if " " in line else None
                current_lines = []
            elif current_tag in tag_set:
                if not line.startswith("QV_SCORE"):
                    current_lines.append(line)

    if current_tag in tag_set and current_lines:
        (output_dir / f"{current_tag}.pdb").write_text("".join(current_lines))


def _write_summary_csv(top: pd.DataFrame, output_dir: Path) -> None:
    cols = [c for c in ["rank", "tag", "pae", "rmsd", "ddg", "composite_score"] if c in top.columns]
    top[cols].to_csv(output_dir / "summary.csv", index=False)


def _write_fasta(output_dir: Path, tags: list[str]) -> None:
    """Write a FASTA file with sequences extracted from exported PDBs."""
    fasta_lines: list[str] = []
    for tag in tags:
        pdb_path = output_dir / f"{tag}.pdb"
        if not pdb_path.exists():
            continue
        h_seq = _extract_sequence_from_pdb(pdb_path, chain="H")
        l_seq = _extract_sequence_from_pdb(pdb_path, chain="L")
        if h_seq:
            fasta_lines.append(f">{tag}_H\n{h_seq}\n")
        if l_seq:
            fasta_lines.append(f">{tag}_L\n{l_seq}\n")

    if fasta_lines:
        (output_dir / "sequences.fasta").write_text("".join(fasta_lines))


def _extract_sequence_from_pdb(pdb_path: Path, chain: str = "H") -> str:
    """Quick-and-dirty sequence extraction from ATOM records of a given chain."""
    seen: dict[int, str] = {}
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("ATOM") and len(line) > 26 and line[21] == chain:
                resnum = int(line[22:26])
                resname = line[17:20].strip()
                if resnum not in seen:
                    seen[resnum] = THREE_TO_ONE.get(resname, "X")
    return "".join(v for _, v in sorted(seen.items()))
