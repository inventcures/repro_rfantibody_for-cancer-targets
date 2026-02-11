"""Shared BioPython utilities used across target_prep and analysis modules."""

from __future__ import annotations

THREE_TO_ONE = {
    "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F",
    "GLY": "G", "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L",
    "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q", "ARG": "R",
    "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y",
}


def find_chain(structure, chain_id: str):
    """Find a chain by ID in a BioPython Structure, or None."""
    for model in structure:
        for chain in model:
            if chain.id == chain_id:
                return chain
    return None
