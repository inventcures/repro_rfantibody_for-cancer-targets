"""Shared subprocess execution for pipeline stages."""

from __future__ import annotations

import logging
import os
import subprocess
import time

logger = logging.getLogger(__name__)


def run_pipeline_command(
    cmd: list[str],
    label: str,
    env: dict[str, str] | None = None,
    timeout: int | None = None,
) -> subprocess.CompletedProcess:
    """Run a subprocess command with merged environment and error handling.

    Streams stdout/stderr to the logger in real-time so per-design metrics
    and stage progress are visible in the log output.

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
    run_env = {**os.environ, "HYDRA_FULL_ERROR": "1", **(env or {})}

    logger.info("%s command: %s", label, " ".join(cmd))

    t0 = time.monotonic()
    proc = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        env=run_env,
        bufsize=1,
    )

    stdout_lines: list[str] = []
    try:
        for line in proc.stdout:
            line = line.rstrip("\n")
            stdout_lines.append(line)
            logger.info("[%s] %s", label, line)
        proc.wait(timeout=timeout)
    except subprocess.TimeoutExpired:
        proc.kill()
        proc.wait()
        raise

    elapsed = time.monotonic() - t0
    logger.info("%s finished in %.1fs (exit %d)", label, elapsed, proc.returncode)

    stdout_text = "\n".join(stdout_lines)

    if proc.returncode != 0:
        raise RuntimeError(
            f"{label} failed (exit {proc.returncode}):\n{stdout_text[-5000:]}"
        )
    return subprocess.CompletedProcess(cmd, proc.returncode, stdout_text, "")
