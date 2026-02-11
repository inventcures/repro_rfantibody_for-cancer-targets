"""Stage 2 wrapper: ProteinMPNN sequence design."""

from __future__ import annotations

import logging
import os
import subprocess
from pathlib import Path

logger = logging.getLogger(__name__)


def run_proteinmpnn(
    input_qv: Path,
    output_qv: Path,
    rfantibody_root: Path,
    sequences_per_backbone: int = 5,
    temperature: float = 0.2,
    env: dict[str, str] | None = None,
) -> Path:
    """Run ProteinMPNN to design CDR sequences for each backbone."""
    script = rfantibody_root / "scripts" / "proteinmpnn_interface_design.py"

    cmd = [
        "python", str(script),
        "-inquiver", str(input_qv),
        "-outquiver", str(output_qv),
    ]

    logger.info(
        "Stage 2 (ProteinMPNN): %d seqs/backbone, temp=%.2f",
        sequences_per_backbone, temperature,
    )
    logger.debug("Command: %s", " ".join(cmd))

    run_env = {**os.environ, **(env or {})}

    result = subprocess.run(cmd, capture_output=True, text=True, env=run_env)
    if result.returncode != 0:
        raise RuntimeError(
            f"ProteinMPNN failed (exit {result.returncode}):\n{result.stderr}"
        )

    logger.info("Stage 2 complete → %s", output_qv)
    return output_qv
