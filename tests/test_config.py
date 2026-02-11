"""Tests for campaign config loading and validation."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import yaml

from harness.config.defaults import BUILTIN_FRAMEWORKS, DEFAULT_CDR_RANGES
from harness.config.schema import (
    AntibodyConfig,
    CampaignConfig,
    CampaignMeta,
    FilteringConfig,
    PipelineConfig,
    TargetConfig,
    load_config,
)


def _minimal_config(**overrides) -> CampaignConfig:
    defaults = dict(
        campaign=CampaignMeta(name="test"),
        target=TargetConfig(
            pdb_id="5N2C",
            chain_id="A",
            epitope_residues=[10, 20, 30, 40, 50],
            hotspot_residues=[10, 20, 30],
        ),
        antibody=AntibodyConfig(
            format="vhh",
            framework="builtin:NbBCII10",
            cdr_loops={"H1": "7", "H2": "6", "H3": "5-13"},
        ),
    )
    defaults.update(overrides)
    return CampaignConfig(**defaults)


class TestConfigValidation(unittest.TestCase):
    def test_valid_minimal(self):
        cfg = _minimal_config()
        self.assertEqual(cfg.validate(), [])

    def test_target_requires_exactly_one_source(self):
        cfg = _minimal_config(target=TargetConfig(
            pdb_id="5N2C", pdb_file="/some/file.pdb",
            epitope_residues=[10, 20, 30], hotspot_residues=[10, 20, 30],
        ))
        errors = cfg.validate()
        self.assertTrue(any("Exactly one" in e for e in errors))

    def test_target_no_source(self):
        cfg = _minimal_config(target=TargetConfig(
            epitope_residues=[10, 20, 30], hotspot_residues=[10, 20, 30],
        ))
        errors = cfg.validate()
        self.assertTrue(any("Exactly one" in e for e in errors))

    def test_hotspot_subset_of_epitope(self):
        cfg = _minimal_config(target=TargetConfig(
            pdb_id="5N2C",
            epitope_residues=[10, 20, 30],
            hotspot_residues=[10, 20, 999],
        ))
        errors = cfg.validate()
        self.assertTrue(any("999" in str(e) for e in errors))

    def test_min_hotspot_count(self):
        cfg = _minimal_config(target=TargetConfig(
            pdb_id="5N2C",
            epitope_residues=[10, 20],
            hotspot_residues=[10, 20],
        ))
        errors = cfg.validate()
        self.assertTrue(any(">= 3" in e for e in errors))

    def test_vhh_cannot_have_light_chain(self):
        cfg = _minimal_config(antibody=AntibodyConfig(
            format="vhh",
            framework="builtin:NbBCII10",
            cdr_loops={"H1": "7", "H2": "6", "H3": "5-13", "L1": "8"},
        ))
        errors = cfg.validate()
        self.assertTrue(any("light chain" in e.lower() for e in errors))

    def test_scfv_needs_both_chains(self):
        cfg = _minimal_config(antibody=AntibodyConfig(
            format="scfv",
            framework="builtin:hu4D5-8",
            cdr_loops={"H1": "7", "H2": "6", "H3": "10"},
        ))
        errors = cfg.validate()
        self.assertTrue(any("both H and L" in e for e in errors))

    def test_scfv_valid(self):
        cfg = _minimal_config(antibody=AntibodyConfig(
            format="scfv",
            framework="builtin:hu4D5-8",
            cdr_loops={"H1": "7", "H2": "6", "H3": "10", "L1": "8", "L2": "7", "L3": "9"},
        ))
        errors = cfg.validate()
        self.assertFalse(any("both H and L" in e for e in errors))

    def test_unknown_framework(self):
        cfg = _minimal_config(antibody=AntibodyConfig(
            format="vhh",
            framework="builtin:Nonexistent",
            cdr_loops={"H1": "7", "H2": "6", "H3": "10"},
        ))
        errors = cfg.validate()
        self.assertTrue(any("Unknown builtin" in e for e in errors))

    def test_num_designs_too_low(self):
        cfg = _minimal_config(pipeline=PipelineConfig())
        cfg.pipeline.rfdiffusion.num_designs = 5
        errors = cfg.validate()
        self.assertTrue(any("num_designs" in e for e in errors))

    def test_temperature_out_of_range(self):
        cfg = _minimal_config(pipeline=PipelineConfig())
        cfg.pipeline.proteinmpnn.temperature = 0.0
        errors = cfg.validate()
        self.assertTrue(any("temperature" in e for e in errors))

    def test_cdr_range_too_large(self):
        cfg = _minimal_config(antibody=AntibodyConfig(
            format="vhh",
            framework="builtin:NbBCII10",
            cdr_loops={"H1": "7", "H2": "6", "H3": "5-30"},
        ))
        errors = cfg.validate()
        self.assertTrue(any("biological limit" in e for e in errors))

    def test_cdr_invalid_range(self):
        cfg = _minimal_config(antibody=AntibodyConfig(
            format="vhh",
            framework="builtin:NbBCII10",
            cdr_loops={"H1": "7", "H2": "6", "H3": "13-5"},
        ))
        errors = cfg.validate()
        self.assertTrue(any("invalid range" in e for e in errors))


class TestConfigLoading(unittest.TestCase):
    def test_load_smoke_test(self):
        cfg = load_config("campaigns/smoke_test.yaml")
        self.assertEqual(cfg.campaign.name, "smoke_test")
        self.assertEqual(cfg.target.pdb_id, "5N2C")
        self.assertEqual(cfg.antibody.format, "vhh")

    def test_load_roundtrip(self):
        data = {
            "campaign": {"name": "roundtrip_test", "description": "test"},
            "target": {
                "pdb_id": "5N2C", "chain_id": "A",
                "epitope_residues": [10, 20, 30, 40, 50],
                "hotspot_residues": [10, 20, 30],
            },
            "antibody": {
                "format": "vhh", "framework": "builtin:NbBCII10",
                "cdr_loops": {"H1": "7", "H2": "6", "H3": "10"},
            },
        }
        with tempfile.NamedTemporaryFile(suffix=".yaml", mode="w", delete=False) as f:
            yaml.dump(data, f)
            f.flush()
            cfg = load_config(f.name)
        self.assertEqual(cfg.campaign.name, "roundtrip_test")
        self.assertEqual(cfg.pipeline.rfdiffusion.num_designs, 10000)


class TestDefaults(unittest.TestCase):
    def test_builtin_frameworks_exist(self):
        self.assertIn("NbBCII10", BUILTIN_FRAMEWORKS)
        self.assertIn("hu4D5-8", BUILTIN_FRAMEWORKS)
        self.assertEqual(BUILTIN_FRAMEWORKS["NbBCII10"]["format"], "vhh")
        self.assertEqual(BUILTIN_FRAMEWORKS["hu4D5-8"]["format"], "scfv")

    def test_default_cdr_ranges(self):
        self.assertIn("H1", DEFAULT_CDR_RANGES)
        self.assertIn("H3", DEFAULT_CDR_RANGES)
        self.assertIn("L1", DEFAULT_CDR_RANGES)


class TestOutputPaths(unittest.TestCase):
    def test_directory_properties(self):
        cfg = _minimal_config()
        self.assertEqual(cfg.output_dir, Path("./results/default"))
        self.assertEqual(cfg.pipeline_dir, Path("./results/default/pipeline"))
        self.assertEqual(cfg.analysis_dir, Path("./results/default/analysis"))
        self.assertEqual(cfg.candidates_dir, Path("./results/default/candidates"))


if __name__ == "__main__":
    unittest.main()
