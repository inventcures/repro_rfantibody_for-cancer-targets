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

    configs = sorted(
        c for c in config_path.glob("*.yaml") if not c.name.startswith("_")
    )
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
        "scipy>=1.7.3",
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
    .add_local_file(
        str(repo_root / "scripts" / "patch_rfantibody.py"),
        remote_path="/tmp/patch_rfantibody.py",
        copy=True,
    )
    .run_commands("python /tmp/patch_rfantibody.py")
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
    timeout=3600,
)
def test_rfantibody_example(use_9lme: bool = False) -> dict:
    """Run RFAntibody example or a custom 9LME test."""
    import subprocess

    env = {**__import__("os").environ, "HYDRA_FULL_ERROR": "1"}

    if use_9lme:
        from urllib.request import urlretrieve
        import os

        os.makedirs(f"{RESULTS_MOUNT}/test_9lme", exist_ok=True)
        raw_pdb = f"{RESULTS_MOUNT}/test_9lme/9LME.pdb"
        urlretrieve("https://files.rcsb.org/download/9LME.pdb", raw_pdb)

        from Bio.PDB import PDBParser, PDBIO, Select

        class CleanChainA(Select):
            def accept_chain(self, chain):
                return chain.id == "A"
            def accept_residue(self, residue):
                return residue.id[0] == " "

        parser = PDBParser(QUIET=True)
        struct = parser.get_structure("9LME", raw_pdb)
        io = PDBIO()
        io.set_structure(struct)
        clean_path = f"{RESULTS_MOUNT}/test_9lme/9LME_chainA_clean.pdb"
        io.save(clean_path, CleanChainA())
        print(f"Extracted clean chain A: {clean_path}")

        cmd = [
            "rfdiffusion",
            "--target", clean_path,
            "--framework", f"{RFANTIBODY_PATH}/scripts/examples/example_inputs/h-NbBCII10.pdb",
            "--output-quiver", f"{RESULTS_MOUNT}/test_9lme/designs.qv",
            "--design-loops", "H1:7,H2:6,H3:5-13",
            "--hotspots", "A33,A54,A76,A101",
            "--diffuser-t", "50",
            "--num-designs", "2",
            "--weights", f"{RFANTIBODY_PATH}/weights/RFdiffusion_Ab.pt",
        ]
    else:
        examples = f"{RFANTIBODY_PATH}/scripts/examples/example_inputs"
        cmd = [
            "python", f"{RFANTIBODY_PATH}/scripts/rfdiffusion_inference.py",
            "--config-name", "antibody",
            f"antibody.target_pdb={examples}/flu_HA.pdb",
            f"antibody.framework_pdb={examples}/h-NbBCII10.pdb",
            f"inference.quiver={RESULTS_MOUNT}/test_example.qv",
            "antibody.design_loops=[H1:7,H2:6,H3:5-13]",
            "ppi.hotspot_res=[B146,B170,B177]",
            "diffuser.T=50",
            "inference.num_designs=2",
            f"inference.ckpt_override_path={RFANTIBODY_PATH}/weights/RFdiffusion_Ab.pt",
        ]

    print(f"CMD: {' '.join(cmd)}")
    proc = subprocess.run(cmd, text=True, capture_output=True, env=env)
    print("STDOUT:", proc.stdout[-5000:] if proc.stdout else "(empty)")
    print("STDERR:", proc.stderr[-5000:] if proc.stderr else "(empty)")
    return {"success": proc.returncode == 0, "returncode": proc.returncode}


@app.function(
    image=rfab_image,
    gpu="A100-80GB",
    volumes={RESULTS_MOUNT: results_volume},
    timeout=86400,
)
def run_campaign(config_yaml: str, campaign_name: str) -> dict:
    """Run a single campaign inside a Modal GPU container.

    Streams full subprocess output to a log file on the Volume and parses
    per-stage timing for post-hoc analysis.
    """
    import os
    import re
    import subprocess
    import tempfile
    import traceback

    import yaml

    logging.basicConfig(
        level=logging.DEBUG,
        format=f"%(asctime)s [{campaign_name}] %(levelname)s %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    clog = logging.getLogger(campaign_name)

    start = time.monotonic()
    stage_times: dict[str, dict] = {}

    try:
        config = yaml.safe_load(config_yaml)
        original_output = config.get("output", {}).get("directory", "")
        patched_output = original_output.replace("./results", RESULTS_MOUNT, 1)
        config["output"]["directory"] = patched_output

        log_dir = Path(patched_output) / "logs"
        log_dir.mkdir(parents=True, exist_ok=True)
        log_path = log_dir / f"{campaign_name}.log"

        clog.info("Output: %s", patched_output)
        clog.info("Log file: %s", log_path)
        clog.info("PDB: %s | Chain: %s | Epitope: %d residues | Hotspots: %s",
                   config.get("target", {}).get("pdb_id"),
                   config.get("target", {}).get("chain_id"),
                   len(config.get("target", {}).get("epitope_residues", [])),
                   config.get("target", {}).get("hotspot_residues"))
        clog.info("Designs: %d | Seqs/backbone: %d | Format: %s",
                   config.get("pipeline", {}).get("rfdiffusion", {}).get("num_designs", 0),
                   config.get("pipeline", {}).get("proteinmpnn", {}).get("sequences_per_backbone", 5),
                   config.get("antibody", {}).get("format"))

        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".yaml", delete=False, prefix=f"{campaign_name}_"
        ) as f:
            yaml.dump(config, f)
            config_path = f.name

        cmd = ["rfab", "--rfantibody-root", RFANTIBODY_PATH, "run", "-v", config_path]
        env = {**os.environ, "PYTHONUNBUFFERED": "1"}

        clog.info("CMD: %s", " ".join(cmd))
        clog.info("Pipeline started...")

        stage_pattern = re.compile(r"=== Stage (\d): (\w+)")
        stage_complete_pattern = re.compile(r"Stage (\d) complete")
        checkpoint_pattern = re.compile(r"Checkpoint saved: stage(\d)")
        current_stage_start: float | None = None
        current_stage_name: str | None = None

        proc = subprocess.Popen(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            text=True, env=env, bufsize=1,
        )

        with open(log_path, "w") as logf:
            logf.write(f"# Campaign: {campaign_name}\n")
            logf.write(f"# Started: {_utcnow()}\n")
            logf.write(f"# Config: {config_path}\n")
            logf.write(f"# Command: {' '.join(cmd)}\n\n")

            for line in proc.stdout:
                line = line.rstrip("\n")
                logf.write(line + "\n")
                clog.info("  %s", line)

                stage_match = stage_pattern.search(line)
                if stage_match:
                    if current_stage_name and current_stage_start is not None:
                        stage_times[current_stage_name]["elapsed_s"] = (
                            time.monotonic() - current_stage_start
                        )
                    stage_num = stage_match.group(1)
                    stage_label = stage_match.group(2)
                    current_stage_name = f"stage{stage_num}_{stage_label}"
                    current_stage_start = time.monotonic()
                    stage_times[current_stage_name] = {
                        "started_at": _utcnow(),
                    }
                    clog.info(">>> STAGE %s (%s) STARTED", stage_num, stage_label)

                if checkpoint_pattern.search(line):
                    if current_stage_name and current_stage_start is not None:
                        elapsed = time.monotonic() - current_stage_start
                        stage_times[current_stage_name]["elapsed_s"] = elapsed
                        stage_times[current_stage_name]["completed_at"] = _utcnow()
                        clog.info(
                            ">>> STAGE %s COMPLETED in %.1fs",
                            current_stage_name, elapsed,
                        )
                        current_stage_name = None
                        current_stage_start = None

            logf.write(f"\n# Finished: {_utcnow()}\n")

        proc.wait()
        elapsed = time.monotonic() - start

        if current_stage_name and current_stage_start is not None:
            stage_times[current_stage_name]["elapsed_s"] = (
                time.monotonic() - current_stage_start
            )

        design_count = _count_designs_on_volume(patched_output)

        result = {
            "campaign": campaign_name,
            "success": proc.returncode == 0,
            "returncode": proc.returncode,
            "elapsed_seconds": elapsed,
            "completed_at": _utcnow(),
            "stage_times": stage_times,
            "log_path": str(log_path),
            "design_count": design_count,
        }

    except Exception as exc:
        elapsed = time.monotonic() - start
        clog.error("Exception: %s\n%s", exc, traceback.format_exc())
        result = {
            "campaign": campaign_name,
            "success": False,
            "returncode": -1,
            "elapsed_seconds": elapsed,
            "completed_at": _utcnow(),
            "stage_times": stage_times,
            "error": f"{exc}\n{traceback.format_exc()}",
        }

    # Save per-campaign result JSON
    try:
        result_json_path = Path(patched_output) / "logs" / f"{campaign_name}_result.json"
        result_json_path.parent.mkdir(parents=True, exist_ok=True)
        result_json_path.write_text(json.dumps(result, indent=2, default=str))
    except Exception as exc:
        clog.warning("Could not write result JSON: %s", exc)

    # Update batch state on Volume
    try:
        state = _load_batch_state_remote()
        state["campaigns"][campaign_name] = {
            "status": "completed" if result["success"] else "failed",
            "elapsed_seconds": result["elapsed_seconds"],
            "completed_at": result["completed_at"],
            "returncode": result.get("returncode", -1),
            "stage_times": stage_times,
            "design_count": result.get("design_count", {}),
        }
        _save_batch_state_remote(state)
        clog.info("Checkpoint saved to Volume")
    except Exception as exc:
        clog.warning("Could not save checkpoint: %s", exc)

    results_volume.commit()

    status = "OK" if result["success"] else "FAILED"
    hours = elapsed / 3600
    clog.info("DONE: %s — %s (%.1fh)", campaign_name, status, hours)
    for sname, sinfo in stage_times.items():
        s_elapsed = sinfo.get("elapsed_s", 0)
        clog.info("  %s: %.1fs (%.1f min)", sname, s_elapsed, s_elapsed / 60)

    return result


def _count_designs_on_volume(output_dir: str) -> dict:
    """Count quiver files and score data in the output directory."""
    from pathlib import Path
    counts: dict[str, int | str] = {}
    try:
        pipeline_dir = Path(output_dir) / "pipeline"
        for qv in pipeline_dir.glob("*.qv"):
            counts[qv.name] = qv.stat().st_size
        analysis_dir = Path(output_dir) / "analysis"
        for f in analysis_dir.glob("*"):
            counts[f"analysis/{f.name}"] = f.stat().st_size
    except Exception:
        pass
    return counts


@app.local_entrypoint()
def main(
    config_dir: str = "campaigns/cancer_drivers",
    campaigns: str = "",
    dry_run: bool = False,
    reset: bool = False,
    test: bool = False,
):
    """Discover configs locally, check for prior progress, dispatch to Modal GPUs."""
    from tqdm import tqdm

    _setup_logging()

    if test:
        use_9lme = campaigns == "9lme"
        label = "9LME chain A" if use_9lme else "flu_HA built-in"
        log.info("Running test: %s + NbBCII10, 2 designs...", label)
        result = test_rfantibody_example.remote(use_9lme=use_9lme)
        status = "PASSED" if result["success"] else "FAILED"
        log.info("Test: %s (exit %d)", status, result["returncode"])
        return

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

    # Write detailed local summary
    summary_dir = Path("results/cancer_drivers_summary")
    summary_dir.mkdir(parents=True, exist_ok=True)

    ts = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S")
    log_path = summary_dir / f"batch_run_{ts}.json"
    all_results = {
        "wall_seconds": wall_elapsed,
        "started_at": datetime.now(timezone.utc).isoformat(),
        "completed_at": _utcnow(),
        "total_campaigns": len(configs),
        "campaigns_run": len(results),
        "campaigns_skipped": len(prior_completed),
        "campaigns_succeeded": sum(1 for r in results if r.get("success")),
        "campaigns_failed": sum(1 for r in results if not r.get("success")),
        "total_gpu_seconds": sum(r.get("elapsed_seconds", 0) for r in results),
        "per_campaign": results,
        "skipped_campaigns": sorted(prior_completed),
    }
    with open(log_path, "w") as f:
        json.dump(all_results, f, indent=2, default=str)

    latest_path = summary_dir / "latest_batch_run.json"
    with open(latest_path, "w") as f:
        json.dump(all_results, f, indent=2, default=str)

    log.info("Local summary: %s", log_path)
    log.info("Download logs: modal volume get rfab-cancer-drivers-results results/")


def _print_summary(state: dict, wall_elapsed: float) -> None:
    campaigns = state.get("campaigns", {})

    print("\n" + "=" * 90)
    print("MODAL BATCH RUN SUMMARY")
    print("=" * 90)
    print(f"{'Campaign':<25} {'Status':<10} {'Total':<8} {'RFdiff':<8} {'MPNN':<8} {'RF2':<8} {'Completed'}")
    print("-" * 90)

    succeeded = 0
    total_gpu_hours = 0.0
    for name in sorted(campaigns.keys()):
        info = campaigns[name]
        status = info.get("status", "unknown").upper()
        elapsed_s = info.get("elapsed_seconds", 0)
        gpu_h = elapsed_s / 3600
        total_gpu_hours += gpu_h

        stages = info.get("stage_times", {})
        rfd_m = stages.get("stage1_RFdiffusion", {}).get("elapsed_s", 0) / 60
        mpnn_m = stages.get("stage2_ProteinMPNN", {}).get("elapsed_s", 0) / 60
        rf2_m = stages.get("stage3_RF2", {}).get("elapsed_s", 0) / 60

        completed_at = info.get("completed_at", "—")
        if completed_at != "—" and len(completed_at) > 19:
            completed_at = completed_at[:19]

        print(
            f"{name:<25} {status:<10} {gpu_h:<8.1f} "
            f"{rfd_m:<8.0f} {mpnn_m:<8.0f} {rf2_m:<8.0f} {completed_at}"
        )
        if status == "COMPLETED":
            succeeded += 1

    print("-" * 90)
    wall_h = wall_elapsed / 3600 if wall_elapsed > 0 else 0
    print(f"Total: {succeeded}/{len(campaigns)} campaigns succeeded")
    if wall_elapsed > 0:
        print(f"  GPU-hours: {total_gpu_hours:.1f}h | Wall time: {wall_h:.1f}h")
    print("=" * 90)


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
