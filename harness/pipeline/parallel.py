"""Multi-GPU parallelization via Quiver splitting."""

from __future__ import annotations

import logging
import os
import subprocess
from pathlib import Path
from typing import Callable

logger = logging.getLogger(__name__)


def run_stage_parallel(
    input_qv: Path,
    stage_fn: Callable[[Path, Path, dict[str, str] | None], Path],
    num_gpus: int,
    work_dir: Path,
    rfantibody_root: Path,
) -> Path:
    """Split a Quiver file across GPUs, run a stage in parallel, merge results.

    Uses upstream ``qvsplit`` and concatenation to distribute work.
    """
    if num_gpus <= 1:
        output = work_dir / f"{input_qv.stem}_out.qv"
        return stage_fn(input_qv, output, None)

    work_dir.mkdir(parents=True, exist_ok=True)
    chunks = _split_quiver(input_qv, num_gpus, work_dir, rfantibody_root)

    outputs: list[Path] = []
    processes: list[subprocess.Popen] = []

    for gpu_id, chunk in enumerate(chunks):
        out = chunk.with_name(f"{chunk.stem}_out.qv")
        outputs.append(out)
        env = {**os.environ, "CUDA_VISIBLE_DEVICES": str(gpu_id)}

        # Each stage_fn returns a Path; here we launch them concurrently
        # by wrapping in a Popen-style approach via threading
        import threading

        def _run(fn, inp, outp, e):
            try:
                fn(inp, outp, e)
            except Exception:
                logger.exception("Parallel stage failed on GPU %d", gpu_id)

        t = threading.Thread(target=_run, args=(stage_fn, chunk, out, env))
        t.start()
        processes.append(t)

    for t in processes:
        t.join()

    merged = work_dir / f"{input_qv.stem}_merged.qv"
    _merge_quivers(outputs, merged)
    return merged


def _split_quiver(
    qv_path: Path, num_chunks: int, work_dir: Path, rfantibody_root: Path
) -> list[Path]:
    """Split a Quiver file into N roughly equal chunks."""
    qvsplit = rfantibody_root / "bin" / "qvsplit"
    if not qvsplit.exists():
        raise FileNotFoundError(f"qvsplit not found at {qvsplit}")

    cmd = [str(qvsplit), str(qv_path), str(num_chunks), str(work_dir / "chunk")]
    logger.info("Splitting %s into %d chunks", qv_path, num_chunks)
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"qvsplit failed:\n{result.stderr}")

    chunks = sorted(work_dir.glob("chunk*.qv"))
    if not chunks:
        raise RuntimeError(f"qvsplit produced no output in {work_dir}")
    return chunks


def _merge_quivers(qv_files: list[Path], output_path: Path) -> Path:
    """Concatenate multiple Quiver files into one."""
    with open(output_path, "w") as out:
        for qv in qv_files:
            if qv.exists():
                out.write(qv.read_text())
    logger.info("Merged %d chunks → %s", len(qv_files), output_path)
    return output_path
