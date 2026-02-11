#!/usr/bin/env python3
"""Modal.com batch runner for cancer driver target campaigns.

Runs campaign configs in parallel on Modal cloud GPUs. Each campaign
gets a dedicated A100-80GB container running the full RFAntibody
3-stage pipeline (RFdiffusion -> ProteinMPNN -> RF2).

Features:
- Idempotent: re-running resumes from last checkpoint (skips completed campaigns)
- Detailed logging with timestamps (local + per-container)
- Progress bar via tqdm
- Batch state persisted to Modal Volume for crash recovery
- Pipeline-level checkpointing (RFAntibody orchestrator saves per-stage checkpoints)

Setup:
    pip install modal tqdm pyyaml

    # Auth (one-time):
    modal token set --token-id <ID> --token-secret <SECRET>

    # GitHub token for private RFAntibody repo (one-time):
    modal secret create github-token GITHUB_TOKEN=ghp_YOUR_TOKEN_HERE

Usage:
    # Local dry-run (no Modal infra needed):
    python scripts/run_cancer_drivers_modal.py --dry-run

    # Full run on Modal GPUs:
    modal run scripts/run_cancer_drivers_modal.py

    # Run subset:
    modal run scripts/run_cancer_drivers_modal.py --campaigns b7h3_vhh,egfrviii_vhh

    # Force re-run (clear checkpoints):
    modal run scripts/run_cancer_drivers_modal.py --reset
"""

from __future__ import annotations

import json
import logging
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

RFANTIBODY_PATH = "/opt/RFAntibody"
RESULTS_MOUNT = "/results"
HARNESS_PATH = "/harness"
BATCH_STATE_FILE = f"{RESULTS_MOUNT}/.batch_state.json"

log = logging.getLogger("rfab-modal")


def _utcnow() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _setup_logging() -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(name)s %(levelname)s %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def _discover_configs(
    config_dir: str, campaigns_filter: str
) -> list[Path]:
    config_path = Path(config_dir)
    if not config_path.exists():
        log.error("Config directory not found: %s", config_path)
        return []

    configs = sorted(config_path.glob("*.yaml"))
    if not configs:
        log.error("No YAML configs in %s", config_path)
        return []

    if campaigns_filter:
        allowed = set(campaigns_filter.split(","))
        configs = [c for c in configs if c.stem in allowed]
        if not configs:
            log.error("No configs matched: %s", campaigns_filter)
            return []

    return configs


def _local_dry_run(config_dir: str, campaigns_filter: str) -> None:
    """Pure local dry-run: discover and validate configs without touching Modal."""
    _setup_logging()
    configs = _discover_configs(config_dir, campaigns_filter)
    if not configs:
        return

    import yaml

    log.info("DRY RUN: %d campaign configs found", len(configs))
    log.info("=" * 60)

    total_designs = 0
    for c in configs:
        raw = yaml.safe_load(c.read_text())
        name = raw.get("campaign", {}).get("name", c.stem)
        pdb = raw.get("target", {}).get("pdb_id", "?")
        n_designs = raw.get("pipeline", {}).get("rfdiffusion", {}).get("num_designs", 10000)
        total_designs += n_designs
        log.info("  %-30s PDB:%-6s designs:%d", name, pdb, n_designs)

    hrs_per_10k = 47
    total_hrs = (total_designs / 10000) * hrs_per_10k
    log.info("=" * 60)
    log.info("Total: %d designs across %d campaigns", total_designs, len(configs))
    log.info("Estimated GPU-hours: %.0f (parallel wall: %.1f h on %d GPUs)",
             total_hrs, total_hrs / len(configs), len(configs))
    log.info("")
    log.info("To run on Modal: modal run scripts/run_cancer_drivers_modal.py")
    log.info("Prerequisites:")
    log.info("  1. modal token set --token-id <ID> --token-secret <SECRET>")
    log.info("  2. modal secret create github-token GITHUB_TOKEN=ghp_...")


# ---------------------------------------------------------------------------
# Modal App, Image, Volume — defined at module level for Modal discovery
# ---------------------------------------------------------------------------
import modal

app = modal.App("rfantibody-harness-cancer-targets")

repo_root = Path(__file__).resolve().parent.parent

rfab_image = (
    modal.Image.from_registry(
        "nvidia/cuda:12.4.0-devel-ubuntu22.04", add_python="3.10"
    )
    .apt_install("git", "wget", "build-essential")
    .pip_install(
        "torch==2.3.1",
        "torchdata==0.7.1",
        "pydantic",
        "psutil>=5.8",
        "scipy>=1.1",
        "requests",
        "pyyaml>=6.0",
        "pandas>=2.0",
        "numpy<2.0",
        "biopython>=1.80",
        "hydra-core>=1.3",
        "tqdm>=4.60",
    )
    .run_commands(
        "pip install dgl --no-deps --no-index"
        " -f https://data.dgl.ai/wheels/torch-2.3/cu121/repo.html",
    )
    .run_commands(
        "git clone https://${GITHUB_TOKEN}@github.com/RosettaCommons/RFantibody.git "
        + RFANTIBODY_PATH,
        f"sed -i 's/cuda-python==11.8/cuda-python>=11.8,<13/' {RFANTIBODY_PATH}/pyproject.toml",
        f"cd {RFANTIBODY_PATH} && pip install --no-deps -e .",
        "pip install biotite click e3nn icecream opt-einsum pyrsistent torch-utils "
        "torchaudio==2.3.1 torchvision==0.18.1 'numpy<2.0'",
        f"mkdir -p {RFANTIBODY_PATH}/weights && cd {RFANTIBODY_PATH}/weights"
        " && wget -q https://files.ipd.uw.edu/pub/RFantibody/RFdiffusion_Ab.pt"
        " && wget -q https://files.ipd.uw.edu/pub/RFantibody/ProteinMPNN_v48_noise_0.2.pt"
        " && wget -q https://files.ipd.uw.edu/pub/RFantibody/RF2_ab.pt",
        secrets=[modal.Secret.from_name("github-token")],
    )
    .add_local_dir(str(repo_root / "harness"), remote_path=f"{HARNESS_PATH}/harness", copy=True)
    .add_local_file(
        str(repo_root / "pyproject.toml"),
        remote_path=f"{HARNESS_PATH}/pyproject.toml",
        copy=True,
    )
    .run_commands(f"pip install {HARNESS_PATH}")
)

results_volume = modal.Volume.from_name(
    "rfab-cancer-drivers-results", create_if_missing=True
)


BATCH_STATE_VOLUME_PATH = ".batch_state.json"


def _load_batch_state_local() -> dict:
    """Load batch state from Volume via client API (for local entrypoint)."""
    try:
        data = b"".join(results_volume.read_file(BATCH_STATE_VOLUME_PATH))
        return json.loads(data.decode())
    except Exception:
        return {"campaigns": {}, "created_at": _utcnow()}


def _reset_batch_state_local() -> None:
    """Delete batch state from Volume via client API (for local entrypoint)."""
    try:
        results_volume.remove_file(BATCH_STATE_VOLUME_PATH)
    except Exception:
        pass


def _load_batch_state_remote() -> dict:
    """Load batch state from mounted Volume filesystem (inside container)."""
    state_path = Path(BATCH_STATE_FILE)
    try:
        results_volume.reload()
        if state_path.exists():
            return json.loads(state_path.read_text())
    except Exception:
        pass
    return {"campaigns": {}, "created_at": _utcnow()}


def _save_batch_state_remote(state: dict) -> None:
    """Persist batch state to mounted Volume filesystem (inside container)."""
    state["updated_at"] = _utcnow()
    Path(BATCH_STATE_FILE).write_text(json.dumps(state, indent=2))
    results_volume.commit()


@app.function(
    image=rfab_image,
    gpu="A100-80GB",
    volumes={RESULTS_MOUNT: results_volume},
    timeout=86400,
)
def run_campaign(config_yaml: str, campaign_name: str) -> dict:
    """Run a single campaign inside a Modal GPU container."""
    import subprocess
    import tempfile
    import traceback

    import yaml

    logging.basicConfig(
        level=logging.INFO,
        format=f"%(asctime)s [{campaign_name}] %(levelname)s %(message)s",
        datefmt="%H:%M:%S",
    )
    clog = logging.getLogger(campaign_name)

    start = time.monotonic()

    try:
        config = yaml.safe_load(config_yaml)
        original_output = config.get("output", {}).get("directory", "")
        patched_output = original_output.replace("./results", RESULTS_MOUNT, 1)
        config["output"]["directory"] = patched_output

        clog.info("Output: %s", patched_output)
        clog.info("PDB: %s | Chain: %s | Epitope: %d residues | Hotspots: %s",
                   config.get("target", {}).get("pdb_id"),
                   config.get("target", {}).get("chain_id"),
                   len(config.get("target", {}).get("epitope_residues", [])),
                   config.get("target", {}).get("hotspot_residues"))
        clog.info("Designs: %d | Format: %s",
                   config.get("pipeline", {}).get("rfdiffusion", {}).get("num_designs", 0),
                   config.get("antibody", {}).get("format"))

        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".yaml", delete=False, prefix=f"{campaign_name}_"
        ) as f:
            yaml.dump(config, f)
            config_path = f.name

        cmd = ["rfab", "--rfantibody-root", RFANTIBODY_PATH, "run", config_path]
        clog.info("Stage 1/3: RFdiffusion (backbone generation)")
        clog.info("Stage 2/3: ProteinMPNN (sequence design)")
        clog.info("Stage 3/3: RF2 (structure prediction)")
        clog.info("Pipeline started...")

        proc = subprocess.run(cmd, text=True, capture_output=True)
        elapsed = time.monotonic() - start

        if proc.stdout:
            for line in proc.stdout.strip().split("\n")[-30:]:
                clog.info("  %s", line)
        if proc.returncode != 0 and proc.stderr:
            for line in proc.stderr.strip().split("\n")[-30:]:
                clog.error("  %s", line)

        result = {
            "campaign": campaign_name,
            "success": proc.returncode == 0,
            "returncode": proc.returncode,
            "elapsed_seconds": elapsed,
            "completed_at": _utcnow(),
            "stdout_tail": proc.stdout[-3000:] if proc.stdout else "",
            "stderr_tail": proc.stderr[-3000:] if proc.stderr else "",
        }

    except Exception as exc:
        elapsed = time.monotonic() - start
        clog.error("Exception: %s", exc)
        result = {
            "campaign": campaign_name,
            "success": False,
            "returncode": -1,
            "elapsed_seconds": elapsed,
            "completed_at": _utcnow(),
            "error": f"{exc}\n{traceback.format_exc()}",
        }

    # Update batch state on Volume
    try:
        state = _load_batch_state_remote()
        state["campaigns"][campaign_name] = {
            "status": "completed" if result["success"] else "failed",
            "elapsed_seconds": result["elapsed_seconds"],
            "completed_at": result["completed_at"],
            "returncode": result.get("returncode", -1),
        }
        _save_batch_state_remote(state)
        clog.info("Checkpoint saved to Volume")
    except Exception as exc:
        clog.warning("Could not save checkpoint: %s", exc)

    results_volume.commit()

    status = "OK" if result["success"] else "FAILED"
    hours = elapsed / 3600
    clog.info("DONE: %s — %s (%.1fh)", campaign_name, status, hours)

    return result


@app.local_entrypoint()
def main(
    config_dir: str = "campaigns/cancer_drivers",
    campaigns: str = "",
    dry_run: bool = False,
    reset: bool = False,
):
    """Discover configs locally, check for prior progress, dispatch to Modal GPUs."""
    from tqdm import tqdm

    _setup_logging()

    configs = _discover_configs(config_dir, campaigns)
    if not configs:
        return

    log.info("Discovered %d campaign configs:", len(configs))
    for c in configs:
        log.info("  - %s", c.stem)

    if dry_run:
        _local_dry_run(config_dir, campaigns)
        return

    # Check for prior progress on Volume (idempotent resume)
    log.info("Checking Volume for prior batch state...")
    prior_state = _load_batch_state_local()
    prior_completed = set()

    if reset:
        log.info("--reset: clearing all prior batch state")
        _reset_batch_state_local()
        prior_state = {"campaigns": {}, "created_at": _utcnow()}
    else:
        for name, info in prior_state.get("campaigns", {}).items():
            if info.get("status") == "completed":
                prior_completed.add(name)
        if prior_completed:
            log.info(
                "Resuming: %d/%d campaigns already completed:",
                len(prior_completed),
                len(configs),
            )
            for name in sorted(prior_completed):
                elapsed_h = prior_state["campaigns"][name].get("elapsed_seconds", 0) / 3600
                log.info("  [SKIP] %-30s (%.1fh)", name, elapsed_h)

    pending_configs = [c for c in configs if c.stem not in prior_completed]

    if not pending_configs:
        log.info("All %d campaigns already completed.", len(configs))
        log.info("Use --reset to force a full re-run.")
        _print_summary(prior_state, 0.0)
        return

    log.info(
        "%d campaigns to run (%d skipped)",
        len(pending_configs),
        len(prior_completed),
    )

    config_args = [(c.read_text(), c.stem) for c in pending_configs]

    log.info("Launching %d campaigns on Modal (A100-80GB)...", len(config_args))
    print("=" * 72)

    wall_start = time.monotonic()

    results = []
    with tqdm(
        total=len(config_args),
        desc="Campaigns",
        unit="campaign",
        bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}]",
    ) as pbar:
        for result in run_campaign.starmap(config_args):
            results.append(result)
            status = "OK" if result["success"] else "FAIL"
            elapsed_h = result["elapsed_seconds"] / 3600
            pbar.set_postfix_str(f"{result['campaign']}: {status} ({elapsed_h:.1f}h)")
            pbar.update(1)

    wall_elapsed = time.monotonic() - wall_start

    final_state = _load_batch_state_local()
    _print_summary(final_state, wall_elapsed)

    # Write local log
    summary_dir = Path("results/cancer_drivers_summary")
    summary_dir.mkdir(parents=True, exist_ok=True)
    log_path = summary_dir / "modal_batch_run.json"
    all_results = {
        "wall_seconds": wall_elapsed,
        "completed_at": _utcnow(),
        "campaigns_run": [r["campaign"] for r in results],
        "campaigns_skipped": sorted(prior_completed),
        "results": results,
    }
    with open(log_path, "w") as f:
        json.dump(all_results, f, indent=2)
    log.info("Local log: %s", log_path)
    log.info("Download: modal volume get rfab-cancer-drivers-results results/")


def _print_summary(state: dict, wall_elapsed: float) -> None:
    campaigns = state.get("campaigns", {})

    print("\n" + "=" * 72)
    print("MODAL BATCH RUN SUMMARY")
    print("=" * 72)
    print(f"{'Campaign':<30} {'Status':<12} {'GPU-hours':<12} {'Completed'}")
    print("-" * 72)

    succeeded = 0
    total_gpu_hours = 0.0
    for name in sorted(campaigns.keys()):
        info = campaigns[name]
        status = info.get("status", "unknown").upper()
        elapsed_s = info.get("elapsed_seconds", 0)
        gpu_h = elapsed_s / 3600
        total_gpu_hours += gpu_h
        completed_at = info.get("completed_at", "—")
        print(f"{name:<30} {status:<12} {gpu_h:<12.1f} {completed_at}")
        if status == "COMPLETED":
            succeeded += 1

    print("-" * 72)
    wall_h = wall_elapsed / 3600 if wall_elapsed > 0 else 0
    print(f"Total: {succeeded}/{len(campaigns)} campaigns succeeded")
    if wall_elapsed > 0:
        print(f"  GPU-hours: {total_gpu_hours:.1f}h | Wall time: {wall_h:.1f}h")
    print("=" * 72)


if __name__ == "__main__":
    if "--dry-run" in sys.argv:
        config_dir = "campaigns/cancer_drivers"
        campaigns_filter = ""
        for i, arg in enumerate(sys.argv):
            if arg == "--config-dir" and i + 1 < len(sys.argv):
                config_dir = sys.argv[i + 1]
            if arg == "--campaigns" and i + 1 < len(sys.argv):
                campaigns_filter = sys.argv[i + 1]
        _local_dry_run(config_dir, campaigns_filter)
    else:
        print("Use: modal run scripts/run_cancer_drivers_modal.py")
        print("  or: python scripts/run_cancer_drivers_modal.py --dry-run")
        sys.exit(1)
