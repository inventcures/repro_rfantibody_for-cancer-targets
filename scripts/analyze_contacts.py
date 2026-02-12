#!/usr/bin/env python3
"""Analyze antibody-antigen contacts in PDB structures.

For each PDB, identifies which residues on the target chain are
within contact distance of the antibody chain(s), and reports
hydrophobic residues suitable as hotspots.

Usage:
    python3 scripts/analyze_contacts.py
"""
from __future__ import annotations

import warnings
from pathlib import Path
from urllib.request import urlretrieve

import numpy as np

warnings.filterwarnings("ignore")

THREE_TO_ONE = {
    "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F",
    "GLY": "G", "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L",
    "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q", "ARG": "R",
    "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y",
}
HYDROPHOBIC = set("AVILMFWP")
CACHE = Path("/tmp/pdb_cache")


def fetch(pdb_id: str, fmt: str = "pdb") -> Path | None:
    CACHE.mkdir(parents=True, exist_ok=True)
    ext = "pdb" if fmt == "pdb" else "cif"
    path = CACHE / f"{pdb_id}.{ext}"
    if not path.exists():
        url = f"https://files.rcsb.org/download/{pdb_id}.{ext}"
        try:
            urlretrieve(url, path)
        except Exception:
            return None
    return path


def load_structure(pdb_id: str):
    from Bio.PDB import PDBParser, MMCIFParser

    pdb_path = fetch(pdb_id, "pdb")
    if pdb_path:
        try:
            parser = PDBParser(QUIET=True, PERMISSIVE=True)
            return parser.get_structure(pdb_id, str(pdb_path))
        except Exception:
            pass

    cif_path = fetch(pdb_id, "cif")
    if cif_path:
        parser = MMCIFParser(QUIET=True)
        return parser.get_structure(pdb_id, str(cif_path))
    return None


def get_contacts(structure, target_chain_id: str, ab_chain_ids: list[str],
                 cutoff: float = 4.5) -> list[dict]:
    """Find target residues within cutoff of any antibody atom."""
    target_chain = None
    ab_chains = []
    for model in structure:
        for chain in model:
            if chain.id == target_chain_id:
                target_chain = chain
            if chain.id in ab_chain_ids:
                ab_chains.append(chain)
        break  # first model only

    if not target_chain or not ab_chains:
        return []

    # Collect all antibody atom coords
    ab_coords = []
    for chain in ab_chains:
        for res in chain.get_residues():
            if res.id[0] != " ":
                continue
            for atom in res.get_atoms():
                if atom.element == "H":
                    continue
                ab_coords.append(atom.get_vector().get_array())
    ab_coords = np.array(ab_coords)

    contacts = []
    for res in target_chain.get_residues():
        if res.id[0] != " ":
            continue
        resid = res.id[1]
        resname = res.get_resname()
        aa = THREE_TO_ONE.get(resname, "?")

        min_dist = float("inf")
        for atom in res.get_atoms():
            if atom.element == "H":
                continue
            coord = atom.get_vector().get_array()
            dists = np.sqrt(np.sum((ab_coords - coord) ** 2, axis=1))
            d = dists.min()
            if d < min_dist:
                min_dist = d

        if min_dist <= cutoff:
            contacts.append({
                "resid": resid,
                "resname": resname,
                "aa": aa,
                "min_dist": min_dist,
                "hydrophobic": aa in HYDROPHOBIC,
            })

    return sorted(contacts, key=lambda x: x["resid"])


def analyze_pdb(pdb_id: str, target_chain: str, ab_chains: list[str],
                description: str = ""):
    print(f"\n{'='*70}")
    print(f"PDB: {pdb_id} | Target: chain {target_chain} | Ab: chains {ab_chains}")
    if description:
        print(f"  {description}")

    structure = load_structure(pdb_id)
    if not structure:
        print(f"  ERROR: Could not load {pdb_id}")
        return

    # Show all chains
    for model in structure:
        for chain in model:
            residues = [r for r in chain.get_residues() if r.id[0] == " "]
            if residues:
                ids = sorted(r.id[1] for r in residues)
                marker = " <-- TARGET" if chain.id == target_chain else ""
                marker += " <-- AB" if chain.id in ab_chains else ""
                print(f"  Chain {chain.id}: {len(residues)} res ({ids[0]}-{ids[-1]}){marker}")
        break

    contacts = get_contacts(structure, target_chain, ab_chains)
    if not contacts:
        print("  No contacts found!")
        return

    epitope_ids = [c["resid"] for c in contacts]
    hydro_contacts = [c for c in contacts if c["hydrophobic"]]

    print(f"\n  Contact residues ({len(contacts)} at 4.5A):")
    for c in contacts:
        h = " *HYDRO*" if c["hydrophobic"] else ""
        print(f"    {c['resid']:4d} {c['aa']}({c['resname']}) dist={c['min_dist']:.1f}A{h}")

    print(f"\n  epitope_residues: {epitope_ids}")
    print(f"  Hydrophobic contacts ({len(hydro_contacts)}):")
    for c in hydro_contacts:
        print(f"    {c['resid']} {c['aa']}({c['resname']}) dist={c['min_dist']:.1f}A")

    # Suggest 4 best hotspots (closest hydrophobic first, then closest polar)
    hydro_sorted = sorted(hydro_contacts, key=lambda x: x["min_dist"])
    if len(hydro_sorted) >= 4:
        suggested = hydro_sorted[:4]
    else:
        all_sorted = sorted(contacts, key=lambda x: (not x["hydrophobic"], x["min_dist"]))
        suggested = all_sorted[:4]

    print(f"\n  Suggested hotspots: {[c['resid'] for c in suggested]}")
    for c in suggested:
        h = " *HYDRO*" if c["hydrophobic"] else ""
        print(f"    {c['resid']} {c['aa']}({c['resname']}) dist={c['min_dist']:.1f}A{h}")


def main():
    # 1. EGFRvIII (8UKV): Chain A target, nanobody chain D (34E5 at junction)
    analyze_pdb("8UKV", "A", ["D"],
                "EGFRvIII + nanobody 34E5 (junction epitope)")

    # 2. GPC2 (6WJL): Chain E target, D3 Fab chains F+J
    analyze_pdb("6WJL", "E", ["F", "J"],
                "GPC2 + D3 Fab")

    # 3. HER2 (1N8Z): Chain C target, trastuzumab Fab chains A+B
    analyze_pdb("1N8Z", "C", ["A", "B"],
                "HER2 domain IV + trastuzumab Fab")

    # 4. Mesothelin N-term (4F3F): Chain C target, amatuximab Fab chains A+B
    analyze_pdb("4F3F", "C", ["A", "B"],
                "Mesothelin N-term peptide + amatuximab Fab")

    # 5. Mesothelin C-term (7U8C): need to identify chains
    analyze_pdb("7U8C", "A", ["H", "L"],
                "Mesothelin C-term + 15B6 Fab (CIF-only)")

    # Also check the PASSING configs for hotspot improvement
    print("\n\n--- HOTSPOT IMPROVEMENT CHECK ---")

    # CD47 (5IWL): Chain C target
    analyze_pdb("5IWL", "C", ["A"],
                "CD47 + SIRPalpha (check hydrophobic hotspots)")

    # CEACAM5 (8BW0): Chain C target
    analyze_pdb("8BW0", "C", ["H", "L"],
                "CEACAM5 + antibody (check hydrophobic hotspots)")

    # EphA2 (3SKJ): Chain E target
    analyze_pdb("3SKJ", "E", ["H", "L"],
                "EphA2 + antibody (check hydrophobic hotspots)")


if __name__ == "__main__":
    main()
