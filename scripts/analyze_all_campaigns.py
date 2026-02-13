#!/usr/bin/env python3
"""Download prediction quivers from Modal Volume and run full analysis for all 10 campaigns.

Outputs:
  results/analysis/per_campaign/{name}/scores.csv     - all design scores
  results/analysis/per_campaign/{name}/ranked.csv     - filtered + ranked candidates
  results/analysis/per_campaign/{name}/summary.json   - stats summary
  results/analysis/cross_campaign_summary.csv         - one row per campaign
  results/analysis/all_scores.csv                     - every design across all campaigns

Usage:
  python scripts/analyze_all_campaigns.py [--skip-download]
"""

import argparse
import json
import subprocess
import sys
from pathlib import Path

VOLUME = "rfab-cancer-drivers-results"
CAMPAIGNS = [
    "b7h3_vhh", "cd47_vhh", "ceacam5_vhh", "egfr_cetuximab_vhh",
    "egfrviii_vhh", "epha2_vhh", "gpc2_vhh", "her2_domIV_vhh",
    "msln_cterm_vhh", "msln_nterm_vhh",
]
LOCAL_BASE = Path("results/analysis")
QV_CACHE = Path("results/qv_cache")

PAE_THRESHOLD = 10.0
RMSD_THRESHOLD = 2.0
DDG_THRESHOLD = -20.0

PAE_WEIGHT = 0.4
RMSD_WEIGHT = 0.3
DDG_WEIGHT = 0.3


def download_quivers(campaigns: list[str]) -> dict[str, Path]:
    QV_CACHE.mkdir(parents=True, exist_ok=True)
    paths = {}
    for name in campaigns:
        remote = f"cancer_drivers/{name}_v1/pipeline/03_predictions.qv"
        local = QV_CACHE / f"{name}_predictions.qv"
        if local.exists() and local.stat().st_size > 0:
            print(f"  {name}: cached ({local.stat().st_size / 1e6:.1f} MB)")
            paths[name] = local
            continue
        print(f"  {name}: downloading from Volume...")
        r = subprocess.run(
            ["modal", "volume", "get", VOLUME, remote, str(local), "--force"],
            capture_output=True, text=True,
        )
        if r.returncode == 0 and local.exists():
            print(f"  {name}: OK ({local.stat().st_size / 1e6:.1f} MB)")
            paths[name] = local
        else:
            print(f"  {name}: FAILED — {r.stderr.strip()}")
    return paths


def download_logs(campaigns: list[str]) -> dict[str, Path]:
    log_dir = LOCAL_BASE / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    paths = {}
    for name in campaigns:
        remote = f"cancer_drivers/{name}_v1/logs/{name}.log"
        local = log_dir / f"{name}.log"
        subprocess.run(
            ["modal", "volume", "get", VOLUME, remote, str(local), "--force"],
            capture_output=True, text=True,
        )
        if local.exists():
            paths[name] = local
    return paths


def parse_scores(qv_path: Path) -> "pd.DataFrame":
    import pandas as pd

    records = []
    current_tag = None
    with open(qv_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith("QV_TAG"):
                current_tag = line.split(None, 1)[1] if " " in line else None
            elif line.startswith("QV_SCORE") and current_tag:
                parts = line.split()
                record = {"tag": current_tag}
                for part in parts[1:]:
                    if "=" in part:
                        k, v = part.split("=", 1)
                        try:
                            record[k] = float(v)
                        except ValueError:
                            record[k] = v
                records.append(record)
    return pd.DataFrame(records) if records else pd.DataFrame(columns=["tag"])


def normalize(series: "pd.Series") -> "pd.Series":
    import pandas as pd
    lo, hi = series.min(), series.max()
    if hi - lo < 1e-12:
        return pd.Series(0.0, index=series.index)
    return (series - lo) / (hi - lo)


def analyze_campaign(name: str, qv_path: Path) -> dict:
    import pandas as pd
    import numpy as np

    out_dir = LOCAL_BASE / "per_campaign" / name
    out_dir.mkdir(parents=True, exist_ok=True)

    scores = parse_scores(qv_path)
    if scores.empty:
        return {"campaign": name, "total_designs": 0, "error": "no scores parsed"}

    pae_col = next((c for c in ["pae", "pae_interaction", "pae_int"] if c in scores.columns), None)
    rmsd_col = next((c for c in ["rmsd", "ca_rmsd", "bb_rmsd"] if c in scores.columns), None)
    ddg_col = next((c for c in ["ddg", "dG_separated", "dG"] if c in scores.columns), None)

    # Standardize column names
    col_map = {}
    if pae_col and pae_col != "pae":
        col_map[pae_col] = "pae"
    if rmsd_col and rmsd_col != "rmsd":
        col_map[rmsd_col] = "rmsd"
    if ddg_col and ddg_col != "ddg":
        col_map[ddg_col] = "ddg"
    if col_map:
        scores = scores.rename(columns=col_map)
        pae_col = "pae" if pae_col else None
        rmsd_col = "rmsd" if rmsd_col else None
        ddg_col = "ddg" if ddg_col else None

    scores.to_csv(out_dir / "scores.csv", index=False)

    total = len(scores)
    summary = {"campaign": name, "total_designs": total, "score_columns": list(scores.columns)}

    for metric in ["pae", "rmsd", "ddg"]:
        if metric in scores.columns and scores[metric].notna().any():
            s = scores[metric].dropna()
            summary[f"{metric}_min"] = float(s.min())
            summary[f"{metric}_q25"] = float(s.quantile(0.25))
            summary[f"{metric}_median"] = float(s.median())
            summary[f"{metric}_mean"] = float(s.mean())
            summary[f"{metric}_q75"] = float(s.quantile(0.75))
            summary[f"{metric}_max"] = float(s.max())
            summary[f"{metric}_std"] = float(s.std())

    # Filtering
    mask = pd.Series(True, index=scores.index)
    if "pae" in scores.columns:
        mask &= scores["pae"] < PAE_THRESHOLD
    if "rmsd" in scores.columns:
        mask &= scores["rmsd"] < RMSD_THRESHOLD
    if ddg_col and "ddg" in scores.columns and scores["ddg"].notna().any():
        mask &= scores["ddg"] < DDG_THRESHOLD

    filtered = scores[mask].copy()
    summary["passed_filters"] = len(filtered)
    summary["pass_rate_pct"] = round(len(filtered) / total * 100, 1) if total else 0

    # Ranking
    if not filtered.empty:
        if "pae" in filtered.columns:
            filtered["pae_norm"] = normalize(filtered["pae"])
        if "rmsd" in filtered.columns:
            filtered["rmsd_norm"] = normalize(filtered["rmsd"])

        composite = pd.Series(0.0, index=filtered.index)
        denom = 0.0
        if "pae_norm" in filtered.columns:
            composite += PAE_WEIGHT * filtered["pae_norm"]
            denom += PAE_WEIGHT
        if "rmsd_norm" in filtered.columns:
            composite += RMSD_WEIGHT * filtered["rmsd_norm"]
            denom += RMSD_WEIGHT
        if "ddg" in filtered.columns and filtered["ddg"].notna().any():
            filtered["ddg_norm"] = normalize(filtered["ddg"])
            composite += DDG_WEIGHT * filtered["ddg_norm"]
            denom += DDG_WEIGHT

        if denom > 0:
            filtered["composite_score"] = composite / denom * (PAE_WEIGHT + RMSD_WEIGHT + DDG_WEIGHT)
        else:
            filtered["composite_score"] = 0.0

        filtered = filtered.sort_values("composite_score").reset_index(drop=True)
        filtered["rank"] = range(1, len(filtered) + 1)
        filtered.to_csv(out_dir / "ranked.csv", index=False)

        top = filtered.head(10)
        summary["top10"] = top[["rank", "tag"] + [c for c in ["pae", "rmsd", "ddg", "composite_score"] if c in top.columns]].to_dict("records")

        if len(filtered) > 0:
            best = filtered.iloc[0]
            summary["best_tag"] = best.get("tag", "")
            summary["best_composite"] = float(best.get("composite_score", 0))

    with open(out_dir / "summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=str)

    return summary


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--skip-download", action="store_true")
    parser.add_argument("--campaigns", nargs="*", default=CAMPAIGNS)
    args = parser.parse_args()

    LOCAL_BASE.mkdir(parents=True, exist_ok=True)

    print("=== Downloading prediction quivers ===")
    if args.skip_download:
        qv_paths = {n: QV_CACHE / f"{n}_predictions.qv" for n in args.campaigns if (QV_CACHE / f"{n}_predictions.qv").exists()}
        print(f"  Skipping download, using {len(qv_paths)} cached files")
    else:
        qv_paths = download_quivers(args.campaigns)

    print(f"\n=== Downloading logs ===")
    log_paths = download_logs(args.campaigns)

    print(f"\n=== Analyzing {len(qv_paths)} campaigns ===")
    import pandas as pd

    all_summaries = []
    all_scores_frames = []

    for name in args.campaigns:
        if name not in qv_paths:
            print(f"  {name}: SKIPPED (no quiver)")
            continue
        print(f"  {name}: analyzing...")
        summary = analyze_campaign(name, qv_paths[name])
        all_summaries.append(summary)
        scores_csv = LOCAL_BASE / "per_campaign" / name / "scores.csv"
        if scores_csv.exists():
            df = pd.read_csv(scores_csv)
            df["campaign"] = name
            all_scores_frames.append(df)
        print(f"    {summary.get('total_designs', 0)} designs, "
              f"{summary.get('passed_filters', 0)} passed ({summary.get('pass_rate_pct', 0)}%), "
              f"best={summary.get('best_composite', 'N/A')}")

    # Cross-campaign summary
    summary_df = pd.DataFrame(all_summaries)
    summary_df.to_csv(LOCAL_BASE / "cross_campaign_summary.csv", index=False)
    with open(LOCAL_BASE / "cross_campaign_summary.json", "w") as f:
        json.dump(all_summaries, f, indent=2, default=str)

    if all_scores_frames:
        all_scores = pd.concat(all_scores_frames, ignore_index=True)
        all_scores.to_csv(LOCAL_BASE / "all_scores.csv", index=False)
        print(f"\n=== Global summary ===")
        print(f"  Total designs scored: {len(all_scores)}")
        total_passed = sum(s.get("passed_filters", 0) for s in all_summaries)
        print(f"  Total passed filters: {total_passed}")
        print(f"  Global pass rate: {total_passed / len(all_scores) * 100:.1f}%")

    print(f"\n=== Outputs ===")
    print(f"  Per-campaign:  {LOCAL_BASE}/per_campaign/*/")
    print(f"  Cross-campaign: {LOCAL_BASE}/cross_campaign_summary.csv")
    print(f"  All scores:     {LOCAL_BASE}/all_scores.csv")


if __name__ == "__main__":
    main()
