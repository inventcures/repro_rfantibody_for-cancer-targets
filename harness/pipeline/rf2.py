"""Stage 3 wrapper: RF2 structure prediction and scoring."""

from __future__ import annotations

import logging
import os
import subprocess
from pathlib import Path

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
    ]

    logger.info("Stage 3 (RF2): recycling=%d", recycling_iterations)
    logger.debug("Command: %s", " ".join(cmd))

    run_env = {**os.environ, **(env or {})}

    result = subprocess.run(cmd, capture_output=True, text=True, env=run_env)
    if result.returncode != 0:
        raise RuntimeError(
            f"RF2 failed (exit {result.returncode}):\n{result.stderr}"
        )

    logger.info("Stage 3 complete → %s", output_qv)
    return output_qv
