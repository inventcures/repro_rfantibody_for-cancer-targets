from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import yaml

from harness.config.defaults import (
    BUILTIN_FRAMEWORKS,
    DEFAULT_CDR_RANGES,
    HYDROPHOBIC_RESIDUES,
)

logger = logging.getLogger(__name__)


@dataclass
class CampaignMeta:
    name: str
    description: str = ""
    version: str = "1.0"


@dataclass
class TruncationConfig:
    enabled: bool = True
    buffer_angstroms: float = 10.0
    preserve_secondary_structure: bool = True


@dataclass
class TargetConfig:
    pdb_id: Optional[str] = None
    pdb_file: Optional[str] = None
    chain_id: str = "A"
    epitope_residues: list[int] = field(default_factory=list)
    hotspot_residues: list[int] = field(default_factory=list)
    truncation: TruncationConfig = field(default_factory=TruncationConfig)


@dataclass
class AntibodyConfig:
    format: str = "vhh"
    framework: str = "builtin:NbBCII10"
    cdr_loops: dict[str, str] = field(default_factory=dict)


@dataclass
class RFdiffusionConfig:
    num_designs: int = 10000
    weights: str = "default"
    seed: Optional[int] = None


@dataclass
class ProteinMPNNConfig:
    sequences_per_backbone: int = 5
    temperature: float = 0.2


@dataclass
class RF2Config:
    recycling_iterations: int = 10
    hotspot_percentage: float = 0.1


@dataclass
class PipelineConfig:
    rfdiffusion: RFdiffusionConfig = field(default_factory=RFdiffusionConfig)
    proteinmpnn: ProteinMPNNConfig = field(default_factory=ProteinMPNNConfig)
    rf2: RF2Config = field(default_factory=RF2Config)


@dataclass
class FilteringConfig:
    pae_threshold: float = 10.0
    rmsd_threshold: float = 2.0
    ddg_threshold: Optional[float] = -20.0


@dataclass
class ComputeConfig:
    gpus: int = 1
    memory_gb: int = 16
    container: str = "local"


@dataclass
class OutputConfig:
    directory: str = "./results/default"
    top_n_candidates: int = 50
    report_format: str = "html"
    export_pdbs: bool = True
    keep_intermediates: bool = False


@dataclass
class ExperimentalConfig:
    enabled: bool = False
    screening_method: str = "yeast_display"
    target_concentrations: list[float] = field(default_factory=lambda: [1400.0, 78.0])
    spr_concentrations: list[float] = field(
        default_factory=lambda: [0.1, 1.0, 10.0, 100.0, 1000.0]
    )


@dataclass
class CampaignConfig:
    campaign: CampaignMeta
    target: TargetConfig
    antibody: AntibodyConfig
    pipeline: PipelineConfig = field(default_factory=PipelineConfig)
    filtering: FilteringConfig = field(default_factory=FilteringConfig)
    compute: ComputeConfig = field(default_factory=ComputeConfig)
    output: OutputConfig = field(default_factory=OutputConfig)
    experimental: ExperimentalConfig = field(default_factory=ExperimentalConfig)

    def validate(self) -> list[str]:
        """Validate config and return list of errors (empty = valid)."""
        errors: list[str] = []

        # Target: exactly one of pdb_id or pdb_file
        has_id = self.target.pdb_id is not None
        has_file = self.target.pdb_file is not None
        if has_id == has_file:
            errors.append("Exactly one of target.pdb_id or target.pdb_file must be set")

        if has_file and not Path(self.target.pdb_file).exists():
            errors.append(f"target.pdb_file not found: {self.target.pdb_file}")

        # Hotspots must be subset of epitope
        if self.target.hotspot_residues and self.target.epitope_residues:
            hotset = set(self.target.hotspot_residues)
            epitset = set(self.target.epitope_residues)
            extra = hotset - epitset
            if extra:
                errors.append(
                    f"Hotspot residues {extra} not in epitope_residues"
                )

        # Minimum hotspot count
        if len(self.target.hotspot_residues) < 3:
            errors.append(
                "target.hotspot_residues must contain >= 3 residues "
                "(need >= 3 hydrophobic for stable binding)"
            )

        # Antibody format validation
        fmt = self.antibody.format.lower()
        if fmt not in ("vhh", "scfv"):
            errors.append(f"antibody.format must be 'vhh' or 'scfv', got '{fmt}'")

        if fmt == "vhh":
            for loop in self.antibody.cdr_loops:
                if loop.startswith("L"):
                    errors.append(
                        f"VHH format cannot have light chain loop {loop}"
                    )

        if fmt == "scfv":
            has_heavy = any(k.startswith("H") for k in self.antibody.cdr_loops)
            has_light = any(k.startswith("L") for k in self.antibody.cdr_loops)
            if not (has_heavy and has_light):
                errors.append("scFv format requires both H and L CDR loops")

        # Framework validation
        fw = self.antibody.framework
        if fw.startswith("builtin:"):
            name = fw.split(":", 1)[1]
            if name not in BUILTIN_FRAMEWORKS:
                errors.append(
                    f"Unknown builtin framework '{name}'. "
                    f"Available: {list(BUILTIN_FRAMEWORKS.keys())}"
                )

        # Pipeline bounds
        if self.pipeline.rfdiffusion.num_designs < 20:
            errors.append("pipeline.rfdiffusion.num_designs must be >= 20")

        if not (0.0 < self.pipeline.proteinmpnn.temperature <= 1.0):
            errors.append("pipeline.proteinmpnn.temperature must be in (0, 1]")

        # CDR loop range validation
        for loop_name, loop_spec in self.antibody.cdr_loops.items():
            if loop_name not in DEFAULT_CDR_RANGES:
                errors.append(f"Unknown CDR loop: {loop_name}")
                continue
            err = _validate_cdr_spec(loop_name, loop_spec)
            if err:
                errors.append(err)

        return errors

    @property
    def output_dir(self) -> Path:
        return Path(self.output.directory)

    @property
    def pipeline_dir(self) -> Path:
        return self.output_dir / "pipeline"

    @property
    def analysis_dir(self) -> Path:
        return self.output_dir / "analysis"

    @property
    def candidates_dir(self) -> Path:
        return self.output_dir / "candidates"


def _validate_cdr_spec(loop_name: str, spec: str) -> Optional[str]:
    """Validate a CDR loop length specification like '7' or '5-13'."""
    spec = spec.strip()
    if "-" in spec:
        parts = spec.split("-")
        if len(parts) != 2:
            return f"{loop_name}: invalid range '{spec}' (use 'min-max')"
        try:
            lo, hi = int(parts[0]), int(parts[1])
        except ValueError:
            return f"{loop_name}: non-integer range '{spec}'"
        if lo < 1 or hi < lo:
            return f"{loop_name}: invalid range {lo}-{hi}"
        if hi > 25:
            return f"{loop_name}: max length {hi} exceeds biological limit (25)"
    else:
        try:
            val = int(spec)
        except ValueError:
            return f"{loop_name}: non-integer length '{spec}'"
        if val < 1 or val > 25:
            return f"{loop_name}: length {val} out of range [1, 25]"
    return None


def _build_dataclass(cls, data: dict):
    """Recursively build a dataclass from a dict, ignoring unknown keys."""
    if not isinstance(data, dict):
        return data
    import dataclasses

    field_names = {f.name for f in dataclasses.fields(cls)}
    filtered = {}
    for k, v in data.items():
        if k not in field_names:
            logger.warning("Ignoring unknown config key: %s", k)
            continue
        fld = next(f for f in dataclasses.fields(cls) if f.name == k)
        if dataclasses.is_dataclass(fld.type):
            filtered[k] = _build_dataclass(fld.type, v)
        elif hasattr(fld.type, "__origin__"):
            filtered[k] = v
        else:
            # Try to resolve string type annotations
            type_map = {
                "TruncationConfig": TruncationConfig,
                "RFdiffusionConfig": RFdiffusionConfig,
                "ProteinMPNNConfig": ProteinMPNNConfig,
                "RF2Config": RF2Config,
            }
            resolved = type_map.get(fld.type if isinstance(fld.type, str) else "", None)
            if resolved and isinstance(v, dict):
                filtered[k] = _build_dataclass(resolved, v)
            else:
                filtered[k] = v
    return cls(**filtered)


def load_config(path: str | Path) -> CampaignConfig:
    """Load and validate a campaign YAML config file."""
    path = Path(path)
    with open(path) as f:
        raw = yaml.safe_load(f)

    campaign_data = raw.get("campaign", {})
    target_data = raw.get("target", {})
    antibody_data = raw.get("antibody", {})
    pipeline_data = raw.get("pipeline", {})
    filtering_data = raw.get("filtering", {})
    compute_data = raw.get("compute", {})
    output_data = raw.get("output", {})
    experimental_data = raw.get("experimental", {})

    config = CampaignConfig(
        campaign=_build_dataclass(CampaignMeta, campaign_data),
        target=_build_dataclass(TargetConfig, target_data),
        antibody=_build_dataclass(AntibodyConfig, antibody_data),
        pipeline=_build_dataclass(PipelineConfig, pipeline_data),
        filtering=_build_dataclass(FilteringConfig, filtering_data),
        compute=_build_dataclass(ComputeConfig, compute_data),
        output=_build_dataclass(OutputConfig, output_data),
        experimental=_build_dataclass(ExperimentalConfig, experimental_data),
    )

    errors = config.validate()
    if errors:
        msg = "Campaign config validation failed:\n" + "\n".join(f"  - {e}" for e in errors)
        raise ValueError(msg)

    return config
