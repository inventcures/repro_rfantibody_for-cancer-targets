"""Main pipeline orchestrator: runs all 3 stages with checkpointing.

Supports batch-level checkpointing for Stages 2 (ProteinMPNN) and 3 (RF2).
If a stage crashes mid-batch, re-running resumes from the last completed batch.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable

from harness.config.schema import CampaignConfig
from harness.pipeline import quiver_utils
from harness.pipeline.proteinmpnn import run_proteinmpnn
from harness.pipeline.rf2 import run_rf2
from harness.pipeline.rfdiffusion import RFdiffusionInputs, run_rfdiffusion

logger = logging.getLogger(__name__)

STAGE2_BATCH_SIZE = 1000
STAGE3_BATCH_SIZE = 5000


@dataclass
class PreparedInputs:
    target_pdb: Path
    framework_hlt: Path
    hotspot_string: str
    cdr_loop_string: str


@dataclass
class PipelineResults:
    backbones: Path
    sequences: Path
    predictions: Path


class PipelineOrchestrator:
    """Runs the 3-stage RFAntibody pipeline with checkpoint/resume support.

    Stages 2 and 3 use batch-level checkpointing: the input quiver is split
    into batches, each batch is processed independently, and a checkpoint
    is saved after each batch.  On resume, completed batches are skipped.
    """

    def __init__(self, config: CampaignConfig, rfantibody_root: Path):
        self.config = config
        self.rfantibody_root = rfantibody_root
        self.work_dir = config.pipeline_dir
        self.work_dir.mkdir(parents=True, exist_ok=True)

    def run(self, inputs: PreparedInputs) -> PipelineResults:
        """Execute all three pipeline stages sequentially."""
        backbones_qv = self.work_dir / "01_backbones.qv"
        sequences_qv = self.work_dir / "02_sequences.qv"
        predictions_qv = self.work_dir / "03_predictions.qv"

        # Stage 1: RFdiffusion (no batching — single monolithic run)
        if not self._checkpoint_exists("stage1"):
            logger.info("=== Stage 1: RFdiffusion ===")
            weights_path = self.rfantibody_root / "weights" / "RFdiffusion_Ab.pt"
            rfd_inputs = RFdiffusionInputs(
                target_pdb=inputs.target_pdb,
                framework_hlt=inputs.framework_hlt,
                hotspot_string=inputs.hotspot_string,
                cdr_loop_string=inputs.cdr_loop_string,
                num_designs=self.config.pipeline.rfdiffusion.num_designs,
                seed=self.config.pipeline.rfdiffusion.seed,
                weights_path=weights_path if weights_path.exists() else None,
            )
            run_rfdiffusion(rfd_inputs, backbones_qv, self.rfantibody_root)
            self._save_checkpoint("stage1")
        else:
            logger.info("Stage 1 checkpoint found — skipping RFdiffusion")

        # Stage 2: ProteinMPNN (batched)
        if not self._checkpoint_exists("stage2"):
            logger.info("=== Stage 2: ProteinMPNN ===")
            mpnn_weights = self.rfantibody_root / "weights" / "ProteinMPNN_v48_noise_0.2.pt"
            self._run_batched_stage(
                stage_name="stage2",
                input_qv=backbones_qv,
                output_qv=sequences_qv,
                batch_size=STAGE2_BATCH_SIZE,
                stage_fn=lambda inp, outp: run_proteinmpnn(
                    inp,
                    outp,
                    self.rfantibody_root,
                    sequences_per_backbone=self.config.pipeline.proteinmpnn.sequences_per_backbone,
                    temperature=self.config.pipeline.proteinmpnn.temperature,
                    weights_path=mpnn_weights if mpnn_weights.exists() else None,
                ),
            )
            self._save_checkpoint("stage2")
        else:
            logger.info("Stage 2 checkpoint found — skipping ProteinMPNN")

        # Stage 3: RF2 (batched)
        if not self._checkpoint_exists("stage3"):
            logger.info("=== Stage 3: RF2 ===")
            rf2_weights = self.rfantibody_root / "weights" / "RF2_ab.pt"
            self._run_batched_stage(
                stage_name="stage3",
                input_qv=sequences_qv,
                output_qv=predictions_qv,
                batch_size=STAGE3_BATCH_SIZE,
                stage_fn=lambda inp, outp: run_rf2(
                    inp,
                    outp,
                    self.rfantibody_root,
                    recycling_iterations=self.config.pipeline.rf2.recycling_iterations,
                    weights_path=rf2_weights if rf2_weights.exists() else None,
                ),
            )
            self._save_checkpoint("stage3")
        else:
            logger.info("Stage 3 checkpoint found — skipping RF2")

        return PipelineResults(
            backbones=backbones_qv,
            sequences=sequences_qv,
            predictions=predictions_qv,
        )

    def _run_batched_stage(
        self,
        stage_name: str,
        input_qv: Path,
        output_qv: Path,
        batch_size: int,
        stage_fn: Callable[[Path, Path], None],
    ) -> None:
        """Split input quiver into batches, process each, checkpoint between."""
        n_entries = quiver_utils.count_entries(input_qv)
        logger.info(
            "%s: %d entries in input, batch_size=%d",
            stage_name, n_entries, batch_size,
        )

        if n_entries <= batch_size:
            logger.info("%s: small enough for single run (no batching)", stage_name)
            stage_fn(input_qv, output_qv)
            return

        batch_dir = self.work_dir / f"{stage_name}_batches"

        # Split input (skip if already split from a previous partial run)
        split_marker = batch_dir / ".split_complete"
        if split_marker.exists():
            batch_inputs = sorted(batch_dir.glob("batch_*.qv"))
            logger.info(
                "%s: found existing split (%d batches)", stage_name, len(batch_inputs),
            )
        else:
            batch_inputs = quiver_utils.split(input_qv, batch_size, batch_dir)
            split_marker.write_text(f"entries:{n_entries} batches:{len(batch_inputs)}")
            logger.info(
                "%s: split into %d batches of ~%d entries",
                stage_name, len(batch_inputs), batch_size,
            )

        n_batches = len(batch_inputs)
        batch_outputs: list[Path] = []
        skipped = 0

        for i, batch_input in enumerate(batch_inputs):
            batch_output = batch_dir / f"output_{i:04d}.qv"
            batch_outputs.append(batch_output)

            batch_cp = batch_dir / f".checkpoint_batch_{i:04d}"
            if batch_cp.exists():
                skipped += 1
                continue

            logger.info(
                "%s batch %d/%d: processing %s",
                stage_name, i + 1, n_batches, batch_input.name,
            )
            stage_fn(batch_input, batch_output)

            ts = datetime.now(timezone.utc).isoformat()
            batch_cp.write_text(f"completed:{ts}")
            logger.info(
                "Batch checkpoint saved: %s_batch_%04d (%d/%d done)",
                stage_name, i, i + 1, n_batches,
            )

        if skipped:
            logger.info("%s: skipped %d already-completed batches", stage_name, skipped)

        quiver_utils.merge(batch_outputs, output_qv)

    def _checkpoint_exists(self, stage: str) -> bool:
        return (self.work_dir / f".checkpoint_{stage}").exists()

    def _save_checkpoint(self, stage: str) -> None:
        ts = datetime.now(timezone.utc).isoformat()
        (self.work_dir / f".checkpoint_{stage}").write_text(f"completed:{ts}")
        logger.info("Checkpoint saved: %s", stage)

    def clear_checkpoints(self) -> None:
        """Remove all checkpoint files to force a full re-run."""
        for cp in self.work_dir.glob(".checkpoint_*"):
            cp.unlink()
            logger.info("Removed checkpoint: %s", cp.name)
        for batch_dir in self.work_dir.glob("stage*_batches"):
            import shutil
            shutil.rmtree(batch_dir)
            logger.info("Removed batch dir: %s", batch_dir.name)
