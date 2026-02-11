"""Main pipeline orchestrator: runs all 3 stages with checkpointing."""

from __future__ import annotations

import logging
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

from harness.config.schema import CampaignConfig
from harness.pipeline.proteinmpnn import run_proteinmpnn
from harness.pipeline.rf2 import run_rf2
from harness.pipeline.rfdiffusion import RFdiffusionInputs, run_rfdiffusion

logger = logging.getLogger(__name__)


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
    """Runs the 3-stage RFAntibody pipeline with checkpoint/resume support."""

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

        # Stage 1: RFdiffusion
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

        # Stage 2: ProteinMPNN
        if not self._checkpoint_exists("stage2"):
            logger.info("=== Stage 2: ProteinMPNN ===")
            run_proteinmpnn(
                backbones_qv,
                sequences_qv,
                self.rfantibody_root,
                sequences_per_backbone=self.config.pipeline.proteinmpnn.sequences_per_backbone,
                temperature=self.config.pipeline.proteinmpnn.temperature,
            )
            self._save_checkpoint("stage2")
        else:
            logger.info("Stage 2 checkpoint found — skipping ProteinMPNN")

        # Stage 3: RF2
        if not self._checkpoint_exists("stage3"):
            logger.info("=== Stage 3: RF2 ===")
            run_rf2(
                sequences_qv,
                predictions_qv,
                self.rfantibody_root,
                recycling_iterations=self.config.pipeline.rf2.recycling_iterations,
            )
            self._save_checkpoint("stage3")
        else:
            logger.info("Stage 3 checkpoint found — skipping RF2")

        return PipelineResults(
            backbones=backbones_qv,
            sequences=sequences_qv,
            predictions=predictions_qv,
        )

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
