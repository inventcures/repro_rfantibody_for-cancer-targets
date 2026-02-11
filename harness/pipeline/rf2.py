"""Stage 3 wrapper: RF2 structure prediction and scoring."""

from __future__ import annotations

import logging
from pathlib import Path

from harness.pipeline._subprocess import run_pipeline_command

logger = logging.getLogger(__name__)


def run_rf2(
    input_qv: Path,
    output_qv: Path,
    rfantibody_root: Path,
    recycling_iterations: int = 10,
    env: dict[str, str] | None = None,
) -> Path:
    """Run RF2 to predict structures and score designed sequences."""
    script = rfantibody_root / "scripts" / "rf2_predict.py"

    cmd = [
        "python", str(script),
        f"input.quiver={input_qv}",
        f"output.quiver={output_qv}",
        f"predict.n_recycles={recycling_iterations}",
    ]

    logger.info("Stage 3 (RF2): recycling=%d", recycling_iterations)

    run_pipeline_command(cmd, "RF2", env=env)

    logger.info("Stage 3 complete → %s", output_qv)
    return output_qv
