"""Prepare antibody framework structures in HLT format."""

from __future__ import annotations

import logging
import shutil
import subprocess
from pathlib import Path

from harness.config.defaults import BUILTIN_FRAMEWORKS

logger = logging.getLogger(__name__)


def prepare_framework(
    framework_source: str,
    antibody_format: str,
    rfantibody_root: Path,
    output_path: Path,
) -> Path:
    """Produce an HLT-format framework PDB ready for RFdiffusion.

    *framework_source* is either ``"builtin:<name>"`` (copies the pre-converted
    file shipped with RFAntibody) or a filesystem path to a Chothia-annotated
    PDB that will be converted via the upstream ``chothia_to_HLT.py`` script.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if framework_source.startswith("builtin:"):
        return _copy_builtin(framework_source, rfantibody_root, output_path)

    return _convert_custom(Path(framework_source), antibody_format, rfantibody_root, output_path)


def _copy_builtin(source_key: str, rfantibody_root: Path, output_path: Path) -> Path:
    name = source_key.split(":", 1)[1]
    meta = BUILTIN_FRAMEWORKS.get(name)
    if meta is None:
        raise ValueError(
            f"Unknown builtin framework '{name}'. "
            f"Available: {list(BUILTIN_FRAMEWORKS.keys())}"
        )

    src = rfantibody_root / meta["relative_path"]
    if not src.exists():
        raise FileNotFoundError(
            f"Builtin framework not found at {src}. "
            "Ensure RFAntibody repo is cloned at the expected location."
        )

    shutil.copy2(src, output_path)
    logger.info("Copied builtin framework %s → %s", name, output_path)
    return output_path


def _convert_custom(
    chothia_pdb: Path,
    antibody_format: str,
    rfantibody_root: Path,
    output_path: Path,
) -> Path:
    if not chothia_pdb.exists():
        raise FileNotFoundError(f"Custom framework PDB not found: {chothia_pdb}")

    converter = rfantibody_root / "scripts" / "util" / "chothia_to_HLT.py"
    if not converter.exists():
        raise FileNotFoundError(
            f"chothia_to_HLT.py not found at {converter}. "
            "Ensure RFAntibody repo is properly installed."
        )

    cmd = [
        "python", str(converter),
        "--input", str(chothia_pdb),
        "--output", str(output_path),
    ]
    if antibody_format == "vhh":
        cmd.append("--nanobody")

    logger.info("Converting framework: %s", " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"chothia_to_HLT.py failed (exit {result.returncode}):\n"
            f"{result.stderr}"
        )

    logger.info("Converted framework → %s", output_path)
    return output_path
