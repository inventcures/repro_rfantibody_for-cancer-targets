"""Rank filtered candidates by composite score."""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import pandas as pd

logger = logging.getLogger(__name__)


def rank_candidates(
    filtered: pd.DataFrame,
    pae_weight: float = 0.4,
    rmsd_weight: float = 0.3,
    ddg_weight: float = 0.3,
) -> pd.DataFrame:
    """Rank candidates by a weighted composite of normalized metrics.

    All metrics are lower-is-better.  Each is min-max normalized to [0, 1],
    then combined with the specified weights.  The composite score is also
    lower-is-better; the returned DataFrame is sorted ascending.
    """
    df = filtered.copy()

    df["pae_norm"] = _normalize(df["pae"])
    df["rmsd_norm"] = _normalize(df["rmsd"])

    composite = pae_weight * df["pae_norm"] + rmsd_weight * df["rmsd_norm"]

    if "ddg" in df.columns and df["ddg"].notna().any():
        df["ddg_norm"] = _normalize(df["ddg"])
        composite += ddg_weight * df["ddg_norm"]
    else:
        # Redistribute ddg weight equally
        composite = (composite / (pae_weight + rmsd_weight))

    df["composite_score"] = composite
    df = df.sort_values("composite_score").reset_index(drop=True)
    df["rank"] = range(1, len(df) + 1)

    logger.info("Ranked %d candidates (best composite=%.4f)", len(df), df["composite_score"].iloc[0] if len(df) else float("nan"))
    return df


def _normalize(series: pd.Series) -> pd.Series:
    """Min-max normalize a series to [0, 1]. Returns 0 if constant."""
    import pandas as pd

    lo, hi = series.min(), series.max()
    if hi - lo < 1e-12:
        return pd.Series(0.0, index=series.index)
    return (series - lo) / (hi - lo)
