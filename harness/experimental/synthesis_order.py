"""Generate gene synthesis orders from top candidate designs."""

from __future__ import annotations

import csv
import logging
from dataclasses import dataclass
from pathlib import Path

logger = logging.getLogger(__name__)

YEAST_CODON_TABLE = {
    "A": "GCT", "C": "TGT", "D": "GAT", "E": "GAA", "F": "TTT",
    "G": "GGT", "H": "CAT", "I": "ATT", "K": "AAG", "L": "TTG",
    "M": "ATG", "N": "AAT", "P": "CCA", "Q": "CAA", "R": "AGA",
    "S": "TCT", "T": "ACT", "V": "GTT", "W": "TGG", "Y": "TAT",
}


@dataclass
class SynthesisEntry:
    name: str
    protein_sequence: str
    dna_sequence: str
    format: str
    vector: str


def generate_synthesis_order(
    fasta_path: Path,
    output_dir: Path,
    antibody_format: str = "vhh",
    vector: str = "pYDS2.0",
) -> Path:
    """Generate a synthesis order spreadsheet from a FASTA of candidate sequences.

    Produces codon-optimized DNA sequences for yeast expression and a CSV
    ready for vendor upload (Twist, IDT, etc.).
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    entries = _parse_fasta(fasta_path)
    orders: list[SynthesisEntry] = []

    for name, seq in entries:
        dna = codon_optimize(seq)
        orders.append(SynthesisEntry(
            name=name,
            protein_sequence=seq,
            dna_sequence=dna,
            format=antibody_format,
            vector=vector,
        ))

    csv_path = output_dir / "synthesis_order.csv"
    _write_csv(orders, csv_path)

    dna_fasta_path = output_dir / "codon_optimized.fasta"
    _write_dna_fasta(orders, dna_fasta_path)

    logger.info("Synthesis order: %d constructs → %s", len(orders), csv_path)
    return csv_path


def codon_optimize(protein_seq: str) -> str:
    """Simple yeast (S. cerevisiae) codon optimization using preferred codons."""
    return "".join(YEAST_CODON_TABLE.get(aa, "NNN") for aa in protein_seq.upper())


def _parse_fasta(path: Path) -> list[tuple[str, str]]:
    entries: list[tuple[str, str]] = []
    name = ""
    seq_parts: list[str] = []

    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if name and seq_parts:
                    entries.append((name, "".join(seq_parts)))
                name = line[1:].split()[0]
                seq_parts = []
            elif line:
                seq_parts.append(line)
    if name and seq_parts:
        entries.append((name, "".join(seq_parts)))

    return entries


def _write_csv(orders: list[SynthesisEntry], path: Path) -> None:
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["name", "format", "vector", "protein_sequence", "dna_sequence", "dna_length_bp"])
        for o in orders:
            writer.writerow([o.name, o.format, o.vector, o.protein_sequence, o.dna_sequence, len(o.dna_sequence)])


def _write_dna_fasta(orders: list[SynthesisEntry], path: Path) -> None:
    with open(path, "w") as f:
        for o in orders:
            f.write(f">{o.name}_codon_opt\n{o.dna_sequence}\n")
