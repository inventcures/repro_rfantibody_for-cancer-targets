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
    weights_path: Path | None = None,
    env: dict[str, str] | None = None,
) -> Path:
    """Run RF2 to predict structures and score designed sequences."""
    cmd = [
        "rf2",
        "--input-quiver", str(input_qv),
        "--output-quiver", str(output_qv),
        "--num-recycles", str(recycling_iterations),
    ]

    if weights_path and weights_path.exists():
        cmd.extend(["--weights", str(weights_path)])

    logger.info("Stage 3 (RF2): recycling=%d", recycling_iterations)

    run_pipeline_command(cmd, "RF2", env=env)

    logger.info("Stage 3 complete → %s", output_qv)
    return output_qv
