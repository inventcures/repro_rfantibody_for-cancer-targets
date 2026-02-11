"""Default parameters and built-in resources for RFAntibody campaigns."""

HYDROPHOBIC_RESIDUES = set("AVLIMFWPY")

DEFAULT_CDR_RANGES: dict[str, str] = {
    "H1": "7",
    "H2": "6",
    "H3": "5-13",
    "L1": "8-13",
    "L2": "7",
    "L3": "9-11",
}

BUILTIN_FRAMEWORKS: dict[str, dict] = {
    "NbBCII10": {
        "filename": "h-NbBCII10.pdb",
        "relative_path": "scripts/examples/example_inputs/h-NbBCII10.pdb",
        "format": "vhh",
        "description": "Humanized VHH scaffold (NBBcH10F-GLA7)",
    },
    "hu4D5-8": {
        "filename": "hu-4D5-8_Fv.pdb",
        "relative_path": "scripts/examples/example_inputs/hu-4D5-8_Fv.pdb",
        "format": "scfv",
        "description": "Trastuzumab-derived scFv framework",
    },
}

RCSB_BASE_URL = "https://files.rcsb.org/download"
