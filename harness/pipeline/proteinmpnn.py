"""Stage 2 wrapper: ProteinMPNN sequence design."""

from __future__ import annotations

import logging
from pathlib import Path

from harness.pipeline._subprocess import run_pipeline_command

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
        "--num_seq_per_target", str(sequences_per_backbone),
        "--sampling_temp", str(temperature),
    ]

    logger.info(
        "Stage 2 (ProteinMPNN): %d seqs/backbone, temp=%.2f",
        sequences_per_backbone, temperature,
    )

    run_pipeline_command(cmd, "ProteinMPNN", env=env)

    logger.info("Stage 2 complete → %s", output_qv)
    return output_qv
