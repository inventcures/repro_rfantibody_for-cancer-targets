"""Apply filtering thresholds to scored designs."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import pandas as pd

logger = logging.getLogger(__name__)

PAE_ALIASES = ["pae", "pae_interaction", "pae_int", "PAE"]
RMSD_ALIASES = ["rmsd", "RMSD", "ca_rmsd", "bb_rmsd"]
DDG_ALIASES = ["ddg", "dG_separated", "ddG", "dG"]


def _find_column(df: pd.DataFrame, aliases: list[str]) -> str | None:
    for alias in aliases:
        if alias in df.columns:
            return alias
    return None


def extract_scores(predictions_qv: Path, rfantibody_root: Path) -> pd.DataFrame:
    """Extract QV_SCORE lines from a predictions Quiver file into a DataFrame.

    Falls back to parsing the Quiver text directly if ``qvscorefile`` is not
    available.
    """
    import shutil

    qvscorefile_bin = shutil.which("qvscorefile")
    if qvscorefile_bin:
        return _extract_via_tool(predictions_qv, Path(qvscorefile_bin))
    return _extract_from_text(predictions_qv)


def apply_filters(
    scores: pd.DataFrame,
    pae_threshold: float = 10.0,
    rmsd_threshold: float = 2.0,
    ddg_threshold: float | None = -20.0,
) -> pd.DataFrame:
    """Filter designs by metric thresholds. Returns rows that pass ALL filters."""
    import pandas as pd

    logger.info("Available score columns: %s", list(scores.columns))

    pae_col = _find_column(scores, PAE_ALIASES)
    rmsd_col = _find_column(scores, RMSD_ALIASES)
    ddg_col = _find_column(scores, DDG_ALIASES)

    mask = pd.Series(True, index=scores.index)

    if pae_col:
        mask &= scores[pae_col] < pae_threshold
    else:
        logger.warning("No PAE column found (tried %s) — skipping pAE filter", PAE_ALIASES)

    if rmsd_col:
        mask &= scores[rmsd_col] < rmsd_threshold
    else:
        logger.warning("No RMSD column found (tried %s) — skipping RMSD filter", RMSD_ALIASES)

    if ddg_threshold is not None and ddg_col:
        mask &= scores[ddg_col] < ddg_threshold

    filtered = scores[mask].copy()
    logger.info(
        "Filtering: %d / %d passed (pAE<%.1f [col=%s], RMSD<%.1f [col=%s]%s)",
        len(filtered),
        len(scores),
        pae_threshold,
        pae_col or "N/A",
        rmsd_threshold,
        rmsd_col or "N/A",
        f", ddG<{ddg_threshold} [col={ddg_col}]" if ddg_threshold is not None and ddg_col else "",
    )
    return filtered


def _extract_via_tool(qv_path: Path, tool_path: Path) -> pd.DataFrame:
    import subprocess

    import pandas as pd

    result = subprocess.run(
        [str(tool_path), str(qv_path)],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(f"qvscorefile failed:\n{result.stderr}")

    import io
    return pd.read_csv(io.StringIO(result.stdout), sep="\t")


def _extract_from_text(qv_path: Path) -> pd.DataFrame:
    """Parse QV_SCORE lines directly from Quiver text format."""
    import pandas as pd

    records: list[dict] = []
    current_tag = None

    with open(qv_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith("QV_TAG"):
                current_tag = line.split(None, 1)[1] if " " in line else None
            elif line.startswith("QV_SCORE"):
                parts = line.split()
                if len(parts) >= 3 and current_tag:
                    record: dict = {"tag": current_tag}
                    for part in parts[1:]:
                        if "=" in part:
                            k, v = part.split("=", 1)
                            try:
                                record[k] = float(v)
                            except ValueError:
                                record[k] = v
                    records.append(record)

    if not records:
        logger.warning("No QV_SCORE lines found in %s", qv_path)
        return pd.DataFrame(columns=["tag"])

    return pd.DataFrame(records)
