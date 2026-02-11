"""Apply filtering thresholds to scored designs."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import pandas as pd

logger = logging.getLogger(__name__)


def extract_scores(predictions_qv: Path, rfantibody_root: Path) -> pd.DataFrame:
    """Extract QV_SCORE lines from a predictions Quiver file into a DataFrame.

    Falls back to parsing the Quiver text directly if ``qvscorefile`` is not
    available.  Expected columns: tag, pae, rmsd, ddg (when present).
    """
    qvscorefile = rfantibody_root / "bin" / "qvscorefile"
    if qvscorefile.exists():
        return _extract_via_tool(predictions_qv, qvscorefile)
    return _extract_from_text(predictions_qv)


def apply_filters(
    scores: pd.DataFrame,
    pae_threshold: float = 10.0,
    rmsd_threshold: float = 2.0,
    ddg_threshold: float | None = -20.0,
) -> pd.DataFrame:
    """Filter designs by metric thresholds. Returns rows that pass ALL filters."""
    mask = (scores["pae"] < pae_threshold) & (scores["rmsd"] < rmsd_threshold)

    if ddg_threshold is not None and "ddg" in scores.columns:
        mask &= scores["ddg"] < ddg_threshold

    filtered = scores[mask].copy()
    logger.info(
        "Filtering: %d / %d passed (pAE<%.1f, RMSD<%.1f%s)",
        len(filtered),
        len(scores),
        pae_threshold,
        rmsd_threshold,
        f", ddG<{ddg_threshold}" if ddg_threshold is not None else "",
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
        return pd.DataFrame(columns=["tag", "pae", "rmsd"])

    return pd.DataFrame(records)
