"""CLI entry point: ``rfab run config.yaml``."""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

from harness.config.schema import CampaignConfig, load_config
from harness.pipeline.orchestrator import PipelineOrchestrator, PreparedInputs
from harness.target_prep.convert_framework import prepare_framework
from harness.target_prep.fetch_pdb import fetch_pdb
from harness.target_prep.truncate import truncate_target
from harness.target_prep.validate import validate_framework_hlt, validate_hotspots

logger = logging.getLogger("harness")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="rfab", description="RFAntibody campaign harness")
    parser.add_argument("--rfantibody-root", type=Path, default=Path("RFAntibody"),
                        help="Path to cloned RFAntibody repository")

    sub = parser.add_subparsers(dest="command")

    # --- run ---
    run_p = sub.add_parser("run", help="Run a full design campaign")
    run_p.add_argument("config", type=Path, help="Campaign YAML config file")
    run_p.add_argument("--dry-run", action="store_true",
                       help="Validate config and prepare inputs only")
    run_p.add_argument("-v", "--verbose", action="store_true")

    # --- validate ---
    val_p = sub.add_parser("validate", help="Validate a campaign config")
    val_p.add_argument("config", type=Path)

    # --- analyze ---
    ana_p = sub.add_parser("analyze", help="Re-run analysis on existing predictions")
    ana_p.add_argument("config", type=Path, help="Campaign YAML config file")
    ana_p.add_argument("-v", "--verbose", action="store_true")

    args = parser.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if getattr(args, "verbose", False) else logging.INFO,
        format="%(asctime)s %(name)s %(levelname)s %(message)s",
    )

    handlers = {
        "validate": _cmd_validate,
        "run": _cmd_run,
        "analyze": _cmd_analyze,
    }

    if args.command is None:
        parser.print_help()
        return 1

    return handlers[args.command](args)


def _cmd_validate(args: argparse.Namespace) -> int:
    try:
        config = load_config(args.config)
    except ValueError as e:
        print(str(e), file=sys.stderr)
        return 1
    print(f"Config valid: {config.campaign.name}")
    return 0


def _cmd_run(args: argparse.Namespace) -> int:
    config = load_config(args.config)
    rfab_root = args.rfantibody_root.resolve()

    logger.info("Campaign: %s", config.campaign.name)

    inputs = _prepare_inputs(config, rfab_root)
    if inputs is None:
        return 1

    if args.dry_run:
        logger.info("Dry run — inputs prepared, stopping before pipeline execution")
        logger.info("  target:    %s", inputs.target_pdb)
        logger.info("  framework: %s", inputs.framework_hlt)
        logger.info("  hotspots:  %s", inputs.hotspot_string)
        logger.info("  CDR loops: %s", inputs.cdr_loop_string)
        return 0

    orchestrator = PipelineOrchestrator(config, rfab_root)
    results = orchestrator.run(inputs)

    logger.info("Pipeline complete → running analysis")
    return _run_analysis(config, results.predictions, rfab_root)


def _cmd_analyze(args: argparse.Namespace) -> int:
    config = load_config(args.config)
    rfab_root = args.rfantibody_root.resolve()

    predictions_qv = config.pipeline_dir / "03_predictions.qv"
    if not predictions_qv.exists():
        logger.error("Predictions file not found: %s", predictions_qv)
        logger.error("Run the pipeline first with: rfab run %s", args.config)
        return 1

    return _run_analysis(config, predictions_qv, rfab_root)


def _run_analysis(config: CampaignConfig, predictions_qv: Path, rfab_root: Path) -> int:
    from harness.analysis.export import export_candidates
    from harness.analysis.filter import apply_filters, extract_scores
    from harness.analysis.rank import rank_candidates
    from harness.analysis.report import generate_report

    logger.info("=== Analysis ===")

    all_scores = extract_scores(predictions_qv, rfab_root)
    logger.info("Extracted scores for %d designs", len(all_scores))

    filtered = apply_filters(
        all_scores,
        pae_threshold=config.filtering.pae_threshold,
        rmsd_threshold=config.filtering.rmsd_threshold,
        ddg_threshold=config.filtering.ddg_threshold,
    )

    if filtered.empty:
        logger.warning("No designs passed filters — report will be empty")

    ranked = rank_candidates(filtered) if not filtered.empty else filtered

    report_path = generate_report(ranked, all_scores, config, config.analysis_dir)
    logger.info("Report → %s", report_path)

    if config.output.export_pdbs and not ranked.empty:
        export_candidates(
            predictions_qv, ranked,
            top_n=config.output.top_n_candidates,
            output_dir=config.candidates_dir,
            rfantibody_root=rfab_root,
        )

    total = len(all_scores)
    passed = len(ranked)
    logger.info("Done: %d / %d designs passed (%.1f%%)", passed, total, passed / total * 100 if total else 0)
    return 0


def _prepare_inputs(config: CampaignConfig, rfab_root: Path) -> PreparedInputs | None:
    """Prepare target + framework inputs, returning None on validation failure."""
    prep_dir = config.output_dir / "prep"
    prep_dir.mkdir(parents=True, exist_ok=True)

    if config.target.pdb_id:
        target_pdb = fetch_pdb(
            config.target.pdb_id, prep_dir, config.target.chain_id,
            source_chain_id=config.target.source_chain_id,
        )
    else:
        target_pdb = Path(config.target.pdb_file)

    if config.target.truncation.enabled and config.target.epitope_residues:
        target_pdb = truncate_target(
            target_pdb,
            config.target.epitope_residues,
            config.target.chain_id,
            buffer_angstroms=config.target.truncation.buffer_angstroms,
            preserve_secondary_structure=config.target.truncation.preserve_secondary_structure,
            output_path=prep_dir / "target_truncated.pdb",
        )

    hotspot_result = validate_hotspots(
        target_pdb, config.target.hotspot_residues, config.target.chain_id,
    )
    if not hotspot_result.valid:
        for err in hotspot_result.errors:
            logger.error("Hotspot validation: %s", err)
        return None
    for w in hotspot_result.warnings:
        logger.warning("Hotspot validation: %s", w)

    framework_hlt = prepare_framework(
        config.antibody.framework,
        config.antibody.format,
        rfab_root,
        prep_dir / "framework.pdb",
    )

    fw_result = validate_framework_hlt(framework_hlt, config.antibody.format)
    if not fw_result.valid:
        for err in fw_result.errors:
            logger.error("Framework validation: %s", err)
        return None

    hotspot_str = ",".join(
        f"{config.target.chain_id}{r}" for r in config.target.hotspot_residues
    )
    cdr_str = ",".join(
        f"{loop}:{spec}" for loop, spec in config.antibody.cdr_loops.items()
    )
    return PreparedInputs(
        target_pdb=target_pdb,
        framework_hlt=framework_hlt,
        hotspot_string=hotspot_str,
        cdr_loop_string=cdr_str,
    )


if __name__ == "__main__":
    sys.exit(main())
