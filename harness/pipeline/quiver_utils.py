"""Quiver (.qv) file utilities for batch splitting, counting, and merging.

Quiver is a text-based format used by RFAntibody to store collections of
protein structures. Each entry starts with a ``QV_TAG {name}`` line followed
by PDB content, optionally followed by ``QV_SCORE`` lines.

These utilities enable batch-level checkpointing: split a large quiver into
batches, process each batch independently, and merge results.
"""

from __future__ import annotations

import logging
from pathlib import Path

logger = logging.getLogger(__name__)

QV_TAG_PREFIX = "QV_TAG"


def count_entries(qv_path: Path) -> int:
    """Count the number of entries (designs) in a quiver file."""
    count = 0
    with open(qv_path) as f:
        for line in f:
            if line.startswith(QV_TAG_PREFIX):
                count += 1
    return count


def split(qv_path: Path, batch_size: int, batch_dir: Path) -> list[Path]:
    """Split a quiver file into batches of ``batch_size`` entries each.

    Returns the list of batch file paths in order.
    """
    batch_dir.mkdir(parents=True, exist_ok=True)
    batches: list[Path] = []
    current_lines: list[str] = []
    current_count = 0
    batch_idx = 0

    def _flush():
        nonlocal current_lines, current_count, batch_idx
        if not current_lines:
            return
        path = batch_dir / f"batch_{batch_idx:04d}.qv"
        path.write_text("".join(current_lines))
        batches.append(path)
        logger.info("  Split batch %d: %d entries → %s", batch_idx, current_count, path.name)
        current_lines = []
        current_count = 0
        batch_idx += 1

    with open(qv_path) as f:
        for line in f:
            if line.startswith(QV_TAG_PREFIX) and current_count >= batch_size:
                _flush()
            current_lines.append(line)
            if line.startswith(QV_TAG_PREFIX):
                current_count += 1

    _flush()
    return batches


def merge(qv_files: list[Path], output_path: Path) -> int:
    """Merge multiple quiver files into one by concatenation.

    Returns the total number of entries in the merged file.
    """
    total = 0
    with open(output_path, "w") as out:
        for qv in qv_files:
            if qv.exists():
                content = qv.read_text()
                out.write(content)
                total += content.count(f"\n{QV_TAG_PREFIX}") + (
                    1 if content.startswith(QV_TAG_PREFIX) else 0
                )
    logger.info("Merged %d files → %s (%d entries)", len(qv_files), output_path.name, total)
    return total
