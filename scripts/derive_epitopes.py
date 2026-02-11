#!/usr/bin/env python3
"""Epitope derivation and verification for cancer driver campaign configs.

Downloads PDB structures, identifies antibody-antigen interface contacts at
a 4.5A heavy-atom distance cutoff, and compares against the hardcoded epitope
residues in the campaign YAML configs.

Usage:
    python scripts/derive_epitopes.py [OPTIONS]

Options:
    --config-dir PATH   Campaign configs to verify (default: campaigns/cancer_drivers/)
    --cutoff FLOAT      Contact distance cutoff in Angstroms (default: 4.5)
    --output PATH       Write results to JSON file
    -v, --verbose       Verbose output
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
import tempfile
import urllib.request
from pathlib import Path

logger = logging.getLogger("derive_epitopes")

REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_CONFIG_DIR = REPO_ROOT / "campaigns" / "cancer_drivers"

RCSB_URL = "https://files.rcsb.org/download"

TARGET_AB_CHAIN_MAP = {
    "b7h3_vhh": {"target_chains": ["A"], "ab_chains": ["B"]},
    "egfrviii_vhh": {"target_chains": ["A"], "ab_chains": ["B"]},
    "egfr_cetuximab_vhh": {"target_chains": ["A"], "ab_chains": ["H", "L"]},
    "gpc2_vhh": {"target_chains": ["A"], "ab_chains": ["H", "L"]},
    "msln_nterm_vhh": {"target_chains": ["A"], "ab_chains": ["H", "L"]},
    "cd47_vhh": {"target_chains": ["C"], "ab_chains": ["A", "B"]},
    "epha2_vhh": {"target_chains": ["E"], "ab_chains": ["H", "L"]},
    "ceacam5_vhh": {"target_chains": ["C"], "ab_chains": ["H", "L"]},
    "her2_domIV_vhh": {"target_chains": ["A"], "ab_chains": ["H", "L"]},
    "msln_cterm_vhh": {"target_chains": ["A"], "ab_chains": ["H", "L"]},
}


def fetch_pdb(pdb_id: str, output_dir: Path) -> Path:
    url = f"{RCSB_URL}/{pdb_id}.pdb"
    output_path = output_dir / f"{pdb_id}.pdb"
    if output_path.exists():
        return output_path
    logger.info("Downloading %s → %s", url, output_path)
    urllib.request.urlretrieve(url, output_path)
    return output_path


def compute_contacts(
    pdb_path: Path,
    target_chains: list[str],
    ab_chains: list[str],
    cutoff: float = 4.5,
) -> dict:
    try:
        from Bio.PDB import PDBParser, NeighborSearch
    except ImportError:
        logger.error("BioPython required. Install with: pip install biopython")
        sys.exit(1)

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("complex", str(pdb_path))
    model = structure[0]

    target_atoms = []
    ab_atoms = []
    for chain in model:
        chain_id = chain.get_id()
        for residue in chain:
            if residue.get_id()[0] != " ":
                continue
            for atom in residue:
                if chain_id in target_chains:
                    target_atoms.append(atom)
                elif chain_id in ab_chains:
                    ab_atoms.append(atom)

    if not target_atoms:
        logger.warning("No target atoms found for chains %s", target_chains)
        return {"epitope_residues": [], "contact_details": []}
    if not ab_atoms:
        logger.warning("No antibody atoms found for chains %s", ab_chains)
        return {"epitope_residues": [], "contact_details": []}

    ns = NeighborSearch(target_atoms)

    contact_residues = set()
    contact_details = []

    for ab_atom in ab_atoms:
        nearby = ns.search(ab_atom.get_vector().get_array(), cutoff, level="R")
        for target_res in nearby:
            res_id = target_res.get_id()[1]
            chain_id = target_res.get_parent().get_id()
            if chain_id in target_chains:
                contact_residues.add(res_id)
                contact_details.append({
                    "residue": res_id,
                    "resname": target_res.get_resname(),
                    "chain": chain_id,
                })

    unique_contacts = sorted(contact_residues)

    unique_details = {}
    for d in contact_details:
        key = d["residue"]
        if key not in unique_details:
            unique_details[key] = d

    return {
        "epitope_residues": unique_contacts,
        "contact_details": [unique_details[r] for r in unique_contacts],
    }


def suggest_hotspots(contact_details: list[dict], max_hotspots: int = 5) -> list[int]:
    HYDROPHOBIC = {"PHE", "TRP", "TYR", "LEU", "ILE", "VAL", "ALA", "MET", "PRO"}

    hydrophobic_contacts = [d for d in contact_details if d["resname"] in HYDROPHOBIC]

    if len(hydrophobic_contacts) >= 3:
        selected = hydrophobic_contacts[:max_hotspots]
    else:
        selected = contact_details[:max_hotspots]

    return sorted(d["residue"] for d in selected)


def load_yaml_epitopes(config_path: Path) -> dict:
    try:
        import yaml
    except ImportError:
        logger.error("PyYAML required")
        sys.exit(1)

    with open(config_path) as f:
        raw = yaml.safe_load(f)

    target = raw.get("target", {})
    return {
        "pdb_id": target.get("pdb_id"),
        "chain_id": target.get("chain_id", "A"),
        "epitope_residues": target.get("epitope_residues", []),
        "hotspot_residues": target.get("hotspot_residues", []),
    }


def compare_epitopes(yaml_residues: list[int], computed_residues: list[int]) -> dict:
    yaml_set = set(yaml_residues)
    computed_set = set(computed_residues)
    overlap = yaml_set & computed_set
    yaml_only = yaml_set - computed_set
    computed_only = computed_set - yaml_set
    jaccard = len(overlap) / len(yaml_set | computed_set) if (yaml_set | computed_set) else 0.0

    return {
        "yaml_count": len(yaml_set),
        "computed_count": len(computed_set),
        "overlap_count": len(overlap),
        "jaccard_similarity": round(jaccard, 3),
        "yaml_only": sorted(yaml_only),
        "computed_only": sorted(computed_only),
        "status": "MATCH" if jaccard > 0.5 else "MISMATCH",
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="derive_epitopes",
        description="Verify campaign config epitopes against PDB structures",
    )
    parser.add_argument("--config-dir", type=Path, default=DEFAULT_CONFIG_DIR)
    parser.add_argument("--cutoff", type=float, default=4.5)
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(name)s %(levelname)s %(message)s",
    )

    configs = sorted(args.config_dir.glob("*.yaml"))
    if not configs:
        logger.error("No configs found in %s", args.config_dir)
        return 1

    logger.info("Verifying %d campaign configs (cutoff=%.1fA)", len(configs), args.cutoff)

    results = {}

    with tempfile.TemporaryDirectory(prefix="epitope_verify_") as tmpdir:
        pdb_dir = Path(tmpdir)

        for config_path in configs:
            campaign_name = config_path.stem
            logger.info("--- %s ---", campaign_name)

            yaml_data = load_yaml_epitopes(config_path)
            pdb_id = yaml_data["pdb_id"]
            if not pdb_id:
                logger.warning("  No pdb_id; skipping")
                continue

            chain_map = TARGET_AB_CHAIN_MAP.get(campaign_name)
            if not chain_map:
                logger.warning("  No chain map for %s; skipping", campaign_name)
                continue

            try:
                pdb_path = fetch_pdb(pdb_id, pdb_dir)
            except Exception as exc:
                logger.error("  Failed to download %s: %s", pdb_id, exc)
                continue

            contacts = compute_contacts(
                pdb_path,
                chain_map["target_chains"],
                chain_map["ab_chains"],
                args.cutoff,
            )

            comparison = compare_epitopes(
                yaml_data["epitope_residues"],
                contacts["epitope_residues"],
            )

            hotspot_suggestions = suggest_hotspots(contacts["contact_details"])

            results[campaign_name] = {
                "pdb_id": pdb_id,
                "target_chain": yaml_data["chain_id"],
                "yaml_epitope": yaml_data["epitope_residues"],
                "computed_epitope": contacts["epitope_residues"],
                "comparison": comparison,
                "yaml_hotspots": yaml_data["hotspot_residues"],
                "suggested_hotspots": hotspot_suggestions,
            }

            status = comparison["status"]
            jaccard = comparison["jaccard_similarity"]
            logger.info("  PDB: %s | YAML: %d res | Computed: %d res | Overlap: %d | Jaccard: %.2f | %s",
                        pdb_id,
                        comparison["yaml_count"],
                        comparison["computed_count"],
                        comparison["overlap_count"],
                        jaccard,
                        status)

            if comparison["yaml_only"]:
                logger.info("  YAML-only residues (may need update): %s", comparison["yaml_only"])
            if comparison["computed_only"]:
                logger.info("  Computed-only (not in YAML): %s", comparison["computed_only"])

    # Summary
    print("\n" + "=" * 72)
    print("EPITOPE VERIFICATION SUMMARY")
    print("=" * 72)
    print(f"{'Campaign':<25} {'PDB':<8} {'YAML':<6} {'Computed':<10} {'Jaccard':<9} {'Status'}")
    print("-" * 72)

    mismatches = 0
    for name, data in results.items():
        comp = data["comparison"]
        print(f"{name:<25} {data['pdb_id']:<8} {comp['yaml_count']:<6} "
              f"{comp['computed_count']:<10} {comp['jaccard_similarity']:<9.3f} {comp['status']}")
        if comp["status"] == "MISMATCH":
            mismatches += 1

    print("-" * 72)
    if mismatches:
        print(f"WARNING: {mismatches} campaigns have epitope mismatches. Update YAMLs before running.")
    else:
        print("All campaigns have acceptable epitope overlap.")
    print("=" * 72)

    if args.output:
        with open(args.output, "w") as f:
            json.dump(results, f, indent=2)
        logger.info("Results written to %s", args.output)

    return 1 if mismatches > 0 else 0


if __name__ == "__main__":
    sys.exit(main())
