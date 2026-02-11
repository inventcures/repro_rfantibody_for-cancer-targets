"""Shared subprocess execution for pipeline stages."""

from __future__ import annotations

import logging
import os
import subprocess
from pathlib import Path

logger = logging.getLogger(__name__)


def run_pipeline_command(
    cmd: list[str],
    label: str,
    env: dict[str, str] | None = None,
    timeout: int | None = None,
) -> subprocess.CompletedProcess:
    """Run a subprocess command with merged environment and error handling.

    Args:
        cmd: Command and arguments.
        label: Human-readable label for error messages (e.g. "RFdiffusion").
        env: Extra environment variables merged on top of os.environ.
        timeout: Optional timeout in seconds.

    Returns:
        The completed process on success.

    Raises:
        RuntimeError: If the subprocess exits with a non-zero return code.
    """
    run_env = {**os.environ, **(env or {})}

    logger.debug("%s command: %s", label, " ".join(cmd))

    result = subprocess.run(
        cmd, capture_output=True, text=True, env=run_env, timeout=timeout,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"{label} failed (exit {result.returncode}):\n{result.stderr}"
        )
    return result
