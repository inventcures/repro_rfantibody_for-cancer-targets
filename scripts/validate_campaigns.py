#!/usr/bin/env python3
"""Validate all campaign YAML configs against actual PDB structures.

Downloads each PDB from RCSB, checks:
1. chain_id exists in the PDB
2. All epitope_residues exist in that chain
3. All hotspot_residues exist in that chain
4. hotspot_residues are a subset of epitope_residues
5. >= 3 hydrophobic hotspot residues
6. PDB format is available (not CIF-only)
7. Residues have Cα atoms (for truncation)

Usage:
    python scripts/validate_campaigns.py [campaigns/cancer_drivers/]
"""
from __future__ import annotations

import sys
from pathlib import Path
from urllib.request import urlretrieve

import yaml

RCSB_PDB_URL = "https://files.rcsb.org/download/{}.pdb"
RCSB_CIF_URL = "https://files.rcsb.org/download/{}.cif"
HYDROPHOBIC = set("AVILMFWP")
THREE_TO_ONE = {
    "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F",
    "GLY": "G", "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L",
    "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q", "ARG": "R",
    "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y",
}


def fetch_pdb(pdb_id: str, cache_dir: Path) -> Path:
    cache_dir.mkdir(parents=True, exist_ok=True)
    path = cache_dir / f"{pdb_id}.pdb"
    if not path.exists():
        url = RCSB_PDB_URL.format(pdb_id)
        try:
            urlretrieve(url, path)
        except Exception:
            cif_path = cache_dir / f"{pdb_id}.cif"
            if not cif_path.exists():
                try:
                    urlretrieve(RCSB_CIF_URL.format(pdb_id), cif_path)
                except Exception as e2:
                    print(f"  ERROR: Could not download {pdb_id} from RCSB: {e2}")
                    return None
            print(f"  NOTE: Using CIF format for {pdb_id} (PDB format unavailable)")
            return cif_path
    return path


def get_structure(pdb_path: Path):
    from Bio.PDB import PDBParser
    import warnings
    warnings.filterwarnings("ignore")
    parser = PDBParser(QUIET=True, PERMISSIVE=True)
    try:
        return parser.get_structure(pdb_path.stem, str(pdb_path))
    except Exception as e:
        # Fall back to mmCIF if PDB parsing fails
        print(f"  WARNING: PDB parse failed ({e}), trying mmCIF...")
        from Bio.PDB import MMCIFParser
        cif_path = pdb_path.parent / f"{pdb_path.stem}.cif"
        if not cif_path.exists():
            from urllib.request import urlretrieve
            urlretrieve(f"https://files.rcsb.org/download/{pdb_path.stem}.cif", cif_path)
        parser2 = MMCIFParser(QUIET=True)
        return parser2.get_structure(pdb_path.stem, str(cif_path))


def validate_config(yaml_path: Path, cache_dir: Path) -> list[str]:
    errors = []
    warnings = []

    raw = yaml.safe_load(yaml_path.read_text())
    name = raw.get("campaign", {}).get("name", yaml_path.stem)
    target = raw.get("target", {})
    antibody = raw.get("antibody", {})
    pipeline = raw.get("pipeline", {})

    pdb_id = target.get("pdb_id", "")
    chain_id = target.get("chain_id", "A")
    source_chain_id = target.get("source_chain_id")
    lookup_chain = source_chain_id or chain_id
    epitope = target.get("epitope_residues", [])
    hotspots = target.get("hotspot_residues", [])
    truncation = target.get("truncation", {})

    print(f"\n{'='*70}")
    print(f"Campaign: {name}")
    chain_label = f"{chain_id} (source: {source_chain_id})" if source_chain_id else chain_id
    print(f"  PDB: {pdb_id} | Chain: {chain_label}")
    print(f"  Epitope: {len(epitope)} residues {epitope[:5]}{'...' if len(epitope) > 5 else ''}")
    print(f"  Hotspots: {hotspots}")
    print(f"  Format: {antibody.get('format')} | CDR: {antibody.get('cdr_loops')}")
    print(f"  Designs: {pipeline.get('rfdiffusion', {}).get('num_designs')}")
    print(f"  Truncation: {'ON' if truncation.get('enabled', True) else 'OFF'}")

    # Check hotspots subset of epitope
    hotset = set(hotspots)
    epitset = set(epitope)
    extra = hotset - epitset
    if extra:
        errors.append(f"Hotspot residues {extra} not in epitope_residues")

    if len(hotspots) < 3:
        errors.append(f"Need >= 3 hotspot residues, got {len(hotspots)}")

    # Download and check PDB
    pdb_path = fetch_pdb(pdb_id, cache_dir)
    if pdb_path is None:
        errors.append(f"PDB {pdb_id} not downloadable in PDB format (may be CIF-only)")
        return errors

    structure = get_structure(pdb_path)

    # List all chains
    all_chains = {}
    for model in structure:
        for chain in model:
            residues = [r for r in chain.get_residues() if r.id[0] == " "]
            if residues:
                res_ids = sorted(r.id[1] for r in residues)
                all_chains[chain.id] = {
                    "count": len(residues),
                    "range": f"{res_ids[0]}-{res_ids[-1]}",
                    "ids": set(res_ids),
                }

    print(f"  PDB chains:")
    for cid, info in sorted(all_chains.items()):
        marker = " <-- TARGET" if cid == lookup_chain else ""
        print(f"    Chain {cid}: {info['count']} residues ({info['range']}){marker}")

    if lookup_chain not in all_chains:
        errors.append(
            f"Chain {lookup_chain} NOT FOUND in {pdb_id}. "
            f"Available: {sorted(all_chains.keys())}"
        )
        return errors

    chain_info = all_chains[lookup_chain]
    chain_res_ids = chain_info["ids"]

    # Check epitope residues exist
    missing_epitope = [r for r in epitope if r not in chain_res_ids]
    present_epitope = [r for r in epitope if r in chain_res_ids]
    if missing_epitope:
        errors.append(
            f"{len(missing_epitope)}/{len(epitope)} epitope residues MISSING from chain {lookup_chain}: "
            f"{missing_epitope}"
        )
        if present_epitope:
            print(f"  Present epitope residues: {present_epitope}")
    else:
        print(f"  All {len(epitope)} epitope residues found in chain {lookup_chain}")

    # Check hotspot residues exist
    missing_hotspots = [r for r in hotspots if r not in chain_res_ids]
    if missing_hotspots:
        errors.append(
            f"Hotspot residues MISSING from chain {lookup_chain}: {missing_hotspots}"
        )

    # Check Cα atoms for epitope residues (needed for truncation)
    chain_obj = None
    for model in structure:
        for c in model:
            if c.id == lookup_chain:
                chain_obj = c
                break

    if chain_obj and present_epitope:
        no_ca = []
        for res in chain_obj.get_residues():
            if res.id[1] in epitset and "CA" not in res:
                no_ca.append(res.id[1])
        if no_ca:
            warnings.append(f"Epitope residues without Cα (may fail truncation): {no_ca}")

    # Check hotspot hydrophobicity
    if chain_obj:
        hydro_count = 0
        hotspot_types = []
        for res in chain_obj.get_residues():
            if res.id[1] in hotset:
                aa = THREE_TO_ONE.get(res.get_resname(), "?")
                hotspot_types.append(f"{res.id[1]}{aa}({res.get_resname()})")
                if aa in HYDROPHOBIC:
                    hydro_count += 1
        if hotspot_types:
            print(f"  Hotspot types: {', '.join(hotspot_types)}")
            if hydro_count < 3:
                warnings.append(
                    f"Only {hydro_count}/4 hydrophobic hotspots "
                    f"(need >= 3 for stable binding)"
                )

    # Print results
    for w in warnings:
        print(f"  WARNING: {w}")
    for e in errors:
        print(f"  ERROR: {e}")

    if not errors:
        print(f"  PASS")
    else:
        print(f"  FAIL ({len(errors)} errors)")

    return errors


def main():
    config_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("campaigns/cancer_drivers")
    cache_dir = Path("/tmp/pdb_cache")

    configs = sorted(c for c in config_dir.glob("*.yaml") if not c.name.startswith("_"))
    print(f"Validating {len(configs)} campaign configs from {config_dir}")

    results = {}
    for config_path in configs:
        errors = validate_config(config_path, cache_dir)
        results[config_path.stem] = errors

    print(f"\n{'='*70}")
    print("SUMMARY")
    print(f"{'='*70}")
    passed = sum(1 for e in results.values() if not e)
    failed = sum(1 for e in results.values() if e)
    for name, errors in sorted(results.items()):
        status = "PASS" if not errors else f"FAIL ({len(errors)} errors)"
        print(f"  {name:<30} {status}")
    print(f"\n{passed}/{len(results)} passed, {failed}/{len(results)} failed")


if __name__ == "__main__":
    main()
