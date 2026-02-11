"""Stage 1 wrapper: RFdiffusion backbone generation."""

from __future__ import annotations

import logging
import subprocess
from dataclasses import dataclass
from pathlib import Path

logger = logging.getLogger(__name__)


@dataclass
class RFdiffusionInputs:
    target_pdb: Path
    framework_hlt: Path
    hotspot_string: str
    cdr_loop_string: str
    num_designs: int
    weights_path: Path | None = None
    seed: int | None = None


def run_rfdiffusion(
    inputs: RFdiffusionInputs,
    output_qv: Path,
    rfantibody_root: Path,
    env: dict[str, str] | None = None,
) -> Path:
    """Run RFdiffusion to generate antibody backbone designs.

    Calls the upstream Hydra-style inference script and writes the output
    to a Quiver (.qv) file.
    """
    script = rfantibody_root / "src" / "rfantibody" / "scripts" / "rfdiffusion_inference.py"

    cmd = [
        "python", str(script),
        "--config-name", "antibody",
        f"antibody.target_pdb={inputs.target_pdb}",
        f"antibody.framework_pdb={inputs.framework_hlt}",
        f"ppi.hotspot_res=[{inputs.hotspot_string}]",
        f"antibody.design_loops=[{inputs.cdr_loop_string}]",
        f"inference.num_designs={inputs.num_designs}",
        f"inference.quiver={output_qv}",
    ]

    if inputs.weights_path:
        cmd.append(f"inference.ckpt_override_path={inputs.weights_path}")

    if inputs.seed is not None:
        cmd.append(f"inference.seed={inputs.seed}")

    logger.info("Stage 1 (RFdiffusion): generating %d backbones", inputs.num_designs)
    logger.debug("Command: %s", " ".join(cmd))

    import os
    run_env = {**os.environ, **(env or {})}

    result = subprocess.run(cmd, capture_output=True, text=True, env=run_env)
    if result.returncode != 0:
        raise RuntimeError(
            f"RFdiffusion failed (exit {result.returncode}):\n{result.stderr}"
        )

    logger.info("Stage 1 complete → %s", output_qv)
    return output_qv
