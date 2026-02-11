#!/usr/bin/env python3
"""Batch runner for cancer driver target campaigns.

Discovers campaign configs in campaigns/cancer_drivers/, validates them,
runs them sequentially or in parallel, and aggregates cross-campaign results.

Usage:
    python scripts/run_cancer_drivers.py [OPTIONS]

Options:
    --config-dir PATH       Campaign config directory (default: campaigns/cancer_drivers/)
    --rfantibody-root PATH  Path to RFAntibody repo (default: ./RFAntibody)
    --parallel N            Run N campaigns concurrently (default: 1)
    --dry-run               Validate and prepare inputs only
    --continue-on-error     Continue remaining campaigns if one fails
    --campaigns LIST        Comma-separated campaign names to run (default: all)
    --summary-dir PATH      Cross-campaign summary output dir
    -v, --verbose           Verbose logging
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict, dataclass
from pathlib import Path

logger = logging.getLogger("cancer_drivers")

REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_CONFIG_DIR = REPO_ROOT / "campaigns" / "cancer_drivers"
DEFAULT_SUMMARY_DIR = REPO_ROOT / "results" / "cancer_drivers_summary"


@dataclass
class CampaignResult:
    name: str
    config_path: str
    returncode: int
    elapsed_seconds: float
    stdout: str
    stderr: str

    @property
    def success(self) -> bool:
        return self.returncode == 0


def discover_configs(config_dir: Path, filter_names: list[str] | None = None) -> list[Path]:
    configs = sorted(config_dir.glob("*.yaml"))
    if not configs:
        logger.error("No YAML configs found in %s", config_dir)
        sys.exit(1)

    if filter_names:
        allowed = set(filter_names)
        configs = [c for c in configs if c.stem in allowed]
        if not configs:
            logger.error("No configs matched filter: %s", filter_names)
            sys.exit(1)

    return configs


def validate_all(configs: list[Path], continue_on_error: bool) -> list[Path]:
    valid = []
    for config_path in configs:
        proc = subprocess.run(
            ["rfab", "validate", str(config_path)],
            capture_output=True,
            text=True,
        )
        if proc.returncode == 0:
            logger.info("  PASS: %s", config_path.stem)
            valid.append(config_path)
        else:
            logger.error("  FAIL: %s — %s", config_path.stem, proc.stderr.strip())
            if not continue_on_error:
                logger.error("Aborting. Fix config or use --continue-on-error.")
                sys.exit(1)
    return valid


def estimate_compute(configs: list[Path]) -> None:
    try:
        import yaml
    except ImportError:
        logger.warning("PyYAML not available; skipping compute estimate")
        return

    total_designs = 0
    for config_path in configs:
        with open(config_path) as f:
            raw = yaml.safe_load(f)
        num = raw.get("pipeline", {}).get("rfdiffusion", {}).get("num_designs", 10000)
        total_designs += num

    hrs_per_10k = 47  # midpoint estimate per campaign
    total_hrs = (total_designs / 10000) * hrs_per_10k
    logger.info("Compute estimate: %d total designs across %d campaigns", total_designs, len(configs))
    logger.info("  Sequential (1 GPU):  ~%.0f GPU-hours (%.1f days)", total_hrs, total_hrs / 24)
    logger.info("  Parallel (4 GPUs):   ~%.0f wall-hours (%.1f days)", total_hrs / 4, total_hrs / 96)
    logger.info("  Parallel (8 GPUs):   ~%.0f wall-hours (%.1f days)", total_hrs / 8, total_hrs / 192)


def run_single_campaign(
    config_path: Path,
    rfantibody_root: Path,
    dry_run: bool,
    verbose: bool,
    env: dict | None = None,
) -> CampaignResult:
    cmd = ["rfab", "run", str(config_path), "--rfantibody-root", str(rfantibody_root)]
    if dry_run:
        cmd.append("--dry-run")
    if verbose:
        cmd.append("-v")

    start = time.monotonic()
    proc = subprocess.run(cmd, capture_output=True, text=True, env=env)
    elapsed = time.monotonic() - start

    return CampaignResult(
        name=config_path.stem,
        config_path=str(config_path),
        returncode=proc.returncode,
        elapsed_seconds=elapsed,
        stdout=proc.stdout,
        stderr=proc.stderr,
    )


def run_sequential(
    configs: list[Path],
    rfantibody_root: Path,
    dry_run: bool,
    verbose: bool,
    continue_on_error: bool,
) -> list[CampaignResult]:
    results = []
    for i, config_path in enumerate(configs):
        logger.info("[%d/%d] Starting: %s", i + 1, len(configs), config_path.stem)

        result = run_single_campaign(config_path, rfantibody_root, dry_run, verbose)
        results.append(result)

        status = "OK" if result.success else "FAILED"
        elapsed_str = format_elapsed(result.elapsed_seconds)
        logger.info("[%d/%d] %s: %s (%s)", i + 1, len(configs), config_path.stem, status, elapsed_str)

        if not result.success:
            logger.error("  stderr: %s", result.stderr.strip()[:500])
            if not continue_on_error:
                logger.error("Aborting. Use --continue-on-error to skip failures.")
                break

    return results


def run_parallel(
    configs: list[Path],
    max_workers: int,
    rfantibody_root: Path,
    dry_run: bool,
    verbose: bool,
) -> list[CampaignResult]:
    results = []

    with ProcessPoolExecutor(max_workers=max_workers) as pool:
        futures = {}
        for i, config_path in enumerate(configs):
            gpu_id = i % max_workers
            env = {**os.environ, "CUDA_VISIBLE_DEVICES": str(gpu_id)}
            future = pool.submit(
                run_single_campaign, config_path, rfantibody_root, dry_run, verbose, env,
            )
            futures[future] = config_path.stem
            logger.info("Submitted: %s → GPU %d", config_path.stem, gpu_id)

        for future in as_completed(futures):
            name = futures[future]
            try:
                result = future.result()
                results.append(result)
                status = "OK" if result.success else "FAILED"
                logger.info("Completed: %s — %s (%s)", name, status, format_elapsed(result.elapsed_seconds))
            except Exception as exc:
                logger.error("Campaign %s raised exception: %s", name, exc)
                results.append(CampaignResult(
                    name=name, config_path="", returncode=1,
                    elapsed_seconds=0, stdout="", stderr=str(exc),
                ))

    return results


def aggregate_cross_campaign(results: list[CampaignResult], summary_dir: Path) -> None:
    summary_dir.mkdir(parents=True, exist_ok=True)

    results_log = []
    for r in results:
        results_log.append({
            "campaign": r.name,
            "success": r.success,
            "elapsed_hours": round(r.elapsed_seconds / 3600, 2),
        })

    log_path = summary_dir / "batch_run.json"
    with open(log_path, "w") as f:
        json.dump(results_log, f, indent=2)
    logger.info("Batch log → %s", log_path)

    try:
        import pandas as pd
    except ImportError:
        logger.warning("pandas not available; skipping CSV aggregation")
        return

    all_candidates = []
    for r in results:
        if not r.success:
            continue
        try:
            import yaml
            with open(r.config_path) as f:
                raw = yaml.safe_load(f)
            output_dir = Path(raw.get("output", {}).get("directory", ""))
            candidates_csv = output_dir / "candidates" / "summary.csv"
            if candidates_csv.exists():
                df = pd.read_csv(candidates_csv)
                df["campaign"] = r.name
                all_candidates.append(df)
        except Exception as exc:
            logger.warning("Could not read candidates for %s: %s", r.name, exc)

    if not all_candidates:
        logger.info("No candidate CSVs found (expected for --dry-run)")
        return

    combined = pd.concat(all_candidates, ignore_index=True)
    if "composite_score" in combined.columns:
        combined.sort_values("composite_score", ascending=True, inplace=True)

    combined.to_csv(summary_dir / "cross_campaign_comparison.csv", index=False)
    logger.info("Cross-campaign comparison → %s", summary_dir / "cross_campaign_comparison.csv")

    stats = combined.groupby("campaign").agg(
        total=("campaign", "count"),
        best_pae=("pae", "min") if "pae" in combined.columns else ("campaign", "count"),
        best_ddg=("ddg", "min") if "ddg" in combined.columns else ("campaign", "count"),
    )
    stats.to_csv(summary_dir / "campaign_stats.csv")
    logger.info("Campaign stats → %s", summary_dir / "campaign_stats.csv")


def format_elapsed(seconds: float) -> str:
    if seconds < 60:
        return f"{seconds:.0f}s"
    if seconds < 3600:
        return f"{seconds / 60:.1f}m"
    return f"{seconds / 3600:.1f}h"


def print_summary(results: list[CampaignResult]) -> None:
    print("\n" + "=" * 72)
    print("BATCH RUN SUMMARY")
    print("=" * 72)
    print(f"{'Campaign':<30} {'Status':<10} {'Time':<12}")
    print("-" * 72)

    succeeded = 0
    for r in results:
        status = "OK" if r.success else "FAILED"
        elapsed = format_elapsed(r.elapsed_seconds)
        print(f"{r.name:<30} {status:<10} {elapsed:<12}")
        if r.success:
            succeeded += 1

    print("-" * 72)
    print(f"Total: {succeeded}/{len(results)} succeeded")
    print("=" * 72)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="run_cancer_drivers",
        description="Batch runner for cancer driver target campaigns",
    )
    parser.add_argument(
        "--config-dir", type=Path, default=DEFAULT_CONFIG_DIR,
        help="Campaign config directory",
    )
    parser.add_argument(
        "--rfantibody-root", type=Path, default=REPO_ROOT / "RFAntibody",
        help="Path to RFAntibody repo",
    )
    parser.add_argument(
        "--parallel", type=int, default=1,
        help="Number of concurrent campaigns (default: 1 = sequential)",
    )
    parser.add_argument("--dry-run", action="store_true", help="Validate and prepare only")
    parser.add_argument("--continue-on-error", action="store_true", help="Skip failed campaigns")
    parser.add_argument(
        "--campaigns", type=str, default=None,
        help="Comma-separated campaign names to run (default: all)",
    )
    parser.add_argument(
        "--summary-dir", type=Path, default=DEFAULT_SUMMARY_DIR,
        help="Cross-campaign summary output directory",
    )
    parser.add_argument("-v", "--verbose", action="store_true")

    args = parser.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(name)s %(levelname)s %(message)s",
        handlers=[
            logging.StreamHandler(),
        ],
    )

    filter_names = args.campaigns.split(",") if args.campaigns else None

    # Phase 1: Discover
    logger.info("=== Phase 1: Discovering campaigns in %s ===", args.config_dir)
    configs = discover_configs(args.config_dir, filter_names)
    logger.info("Found %d campaign configs", len(configs))
    for c in configs:
        logger.info("  - %s", c.stem)

    # Phase 2: Validate
    logger.info("=== Phase 2: Validating all configs ===")
    valid_configs = validate_all(configs, args.continue_on_error)
    logger.info("%d/%d configs valid", len(valid_configs), len(configs))

    if not valid_configs:
        logger.error("No valid configs to run")
        return 1

    # Phase 3: Compute estimate
    logger.info("=== Phase 3: Compute estimate ===")
    estimate_compute(valid_configs)

    if args.dry_run:
        logger.info("=== Dry run — running rfab dry-run for each config ===")

    # Phase 4: Execute
    logger.info("=== Phase 4: Executing %d campaigns (%s) ===",
                len(valid_configs), f"parallel={args.parallel}" if args.parallel > 1 else "sequential")

    if args.parallel > 1:
        results = run_parallel(
            valid_configs, args.parallel, args.rfantibody_root, args.dry_run, args.verbose,
        )
    else:
        results = run_sequential(
            valid_configs, args.rfantibody_root, args.dry_run, args.verbose, args.continue_on_error,
        )

    # Phase 5: Aggregate
    if not args.dry_run:
        logger.info("=== Phase 5: Aggregating cross-campaign results ===")
        aggregate_cross_campaign(results, args.summary_dir)

    # Summary
    print_summary(results)

    # Add file handler for the log
    args.summary_dir.mkdir(parents=True, exist_ok=True)
    log_data = [{"campaign": r.name, "success": r.success, "elapsed_s": r.elapsed_seconds} for r in results]
    with open(args.summary_dir / "batch_run.json", "w") as f:
        json.dump(log_data, f, indent=2)

    failed = sum(1 for r in results if not r.success)
    return 1 if failed > 0 else 0


if __name__ == "__main__":
    sys.exit(main())
