"""Tests for issues uncovered by the code audit."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import yaml

try:
    import pandas as pd
    HAS_PANDAS = True
except ImportError:
    HAS_PANDAS = False


class TestEmptyYamlLoading(unittest.TestCase):
    def test_empty_file_raises_valueerror(self):
        from harness.config.schema import load_config

        with tempfile.NamedTemporaryFile(suffix=".yaml", mode="w", delete=False) as f:
            f.write("")
            f.flush()
            with self.assertRaises(ValueError) as ctx:
                load_config(f.name)
            self.assertIn("expected YAML mapping", str(ctx.exception))

    def test_non_dict_yaml_raises_valueerror(self):
        from harness.config.schema import load_config

        with tempfile.NamedTemporaryFile(suffix=".yaml", mode="w", delete=False) as f:
            f.write("- just\n- a\n- list\n")
            f.flush()
            with self.assertRaises(ValueError) as ctx:
                load_config(f.name)
            self.assertIn("expected YAML mapping", str(ctx.exception))


class TestTargetIdProperty(unittest.TestCase):
    def test_target_id_from_pdb_id(self):
        from harness.config.schema import TargetConfig

        t = TargetConfig(pdb_id="5N2C")
        self.assertEqual(t.target_id, "5N2C")

    def test_target_id_from_pdb_file(self):
        from harness.config.schema import TargetConfig

        t = TargetConfig(pdb_file="/some/path/my_target.pdb")
        self.assertEqual(t.target_id, "my_target")

    def test_target_id_unknown(self):
        from harness.config.schema import TargetConfig

        t = TargetConfig()
        self.assertEqual(t.target_id, "unknown")


class TestSubprocessRunner(unittest.TestCase):
    def test_success(self):
        from harness.pipeline._subprocess import run_pipeline_command

        result = run_pipeline_command(["echo", "hello"], "TestCmd")
        self.assertEqual(result.returncode, 0)
        self.assertIn("hello", result.stdout)

    def test_failure_raises_runtime_error(self):
        from harness.pipeline._subprocess import run_pipeline_command

        with self.assertRaises(RuntimeError) as ctx:
            run_pipeline_command(["false"], "TestFail")
        self.assertIn("TestFail failed", str(ctx.exception))

    def test_env_merge(self):
        from harness.pipeline._subprocess import run_pipeline_command

        result = run_pipeline_command(
            ["env"], "TestEnv", env={"MY_TEST_VAR": "42"}
        )
        self.assertIn("MY_TEST_VAR=42", result.stdout)


class TestCLIValidate(unittest.TestCase):
    def test_validate_smoke_test(self):
        from harness.cli import main

        ret = main(["validate", "campaigns/smoke_test.yaml"])
        self.assertEqual(ret, 0)

    def test_validate_invalid_file(self):
        with tempfile.NamedTemporaryFile(suffix=".yaml", mode="w", delete=False) as f:
            f.write("")
            f.flush()
            from harness.cli import main
            ret = main(["validate", f.name])
            self.assertEqual(ret, 1)

    def test_no_command_returns_1(self):
        from harness.cli import main

        ret = main([])
        self.assertEqual(ret, 1)


class TestFetchPdbCleanup(unittest.TestCase):
    @patch("harness.target_prep.fetch_pdb.urlretrieve")
    def test_partial_download_cleaned_up(self, mock_retrieve):
        from harness.target_prep.fetch_pdb import fetch_pdb

        def fake_download(url, path):
            Path(path).write_text("partial")
            raise ConnectionError("network failure")

        mock_retrieve.side_effect = fake_download
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(ConnectionError):
                fetch_pdb("1ABC", Path(tmp))
            self.assertFalse((Path(tmp) / "1ABC.pdb").exists())


class TestBioUtils(unittest.TestCase):
    def test_three_to_one_complete(self):
        from harness._bio_utils import THREE_TO_ONE

        self.assertEqual(len(THREE_TO_ONE), 20)
        self.assertEqual(THREE_TO_ONE["ALA"], "A")
        self.assertEqual(THREE_TO_ONE["TYR"], "Y")

    def test_find_chain_returns_none_for_missing(self):
        from harness._bio_utils import find_chain

        class MockChain:
            def __init__(self, cid):
                self.id = cid

        class MockModel:
            def __init__(self):
                self.chains = [MockChain("A"), MockChain("B")]

            def __iter__(self):
                return iter(self.chains)

        result = find_chain([MockModel()], "Z")
        self.assertIsNone(result)


@unittest.skipUnless(HAS_PANDAS, "pandas not installed")
class TestRankingWithDdg(unittest.TestCase):
    def setUp(self):
        from harness.analysis.rank import rank_candidates
        self.rank_candidates = rank_candidates

    def test_rank_with_ddg(self):
        scores = pd.DataFrame({
            "tag": ["d001", "d002", "d003"],
            "pae": [5.0, 8.0, 3.0],
            "rmsd": [1.0, 1.5, 0.5],
            "ddg": [-25.0, -15.0, -30.0],
        })
        ranked = self.rank_candidates(scores)
        self.assertIn("ddg_norm", ranked.columns)
        self.assertEqual(ranked.iloc[0]["tag"], "d003")

    def test_rank_empty_dataframe(self):
        scores = pd.DataFrame(columns=["tag", "pae", "rmsd"])
        ranked = self.rank_candidates(scores)
        self.assertEqual(len(ranked), 0)


@unittest.skipUnless(HAS_PANDAS, "pandas not installed")
class TestReportGeneration(unittest.TestCase):
    def test_csv_report(self):
        from harness.analysis.report import generate_report
        from harness.config.schema import (
            AntibodyConfig, CampaignConfig, CampaignMeta, TargetConfig,
        )

        config = CampaignConfig(
            campaign=CampaignMeta(name="test_report"),
            target=TargetConfig(
                pdb_id="5N2C", epitope_residues=[10, 20, 30, 40, 50],
                hotspot_residues=[10, 20, 30],
            ),
            antibody=AntibodyConfig(
                format="vhh", framework="builtin:NbBCII10",
                cdr_loops={"H1": "7", "H2": "6", "H3": "10"},
            ),
        )

        ranked = pd.DataFrame({
            "tag": ["d001"], "pae": [5.0], "rmsd": [1.0], "rank": [1],
            "composite_score": [0.0],
        })
        all_scores = pd.DataFrame({
            "tag": ["d001", "d002"], "pae": [5.0, 8.0], "rmsd": [1.0, 2.0],
        })

        with tempfile.TemporaryDirectory() as tmp:
            config.output.report_format = "csv"
            result = generate_report(ranked, all_scores, config, Path(tmp))
            self.assertTrue(result.exists())
            self.assertTrue(result.name.endswith(".csv"))

    def test_html_escapes_special_chars(self):
        from harness.analysis.report import generate_report
        from harness.config.schema import (
            AntibodyConfig, CampaignConfig, CampaignMeta, TargetConfig,
        )

        config = CampaignConfig(
            campaign=CampaignMeta(name="<script>alert('xss')</script>"),
            target=TargetConfig(
                pdb_id="5N2C", epitope_residues=[10, 20, 30, 40, 50],
                hotspot_residues=[10, 20, 30],
            ),
            antibody=AntibodyConfig(
                format="vhh", framework="builtin:NbBCII10",
                cdr_loops={"H1": "7", "H2": "6", "H3": "10"},
            ),
        )

        ranked = pd.DataFrame({
            "tag": ["d001"], "pae": [5.0], "rmsd": [1.0], "rank": [1],
            "composite_score": [0.0],
        })
        all_scores = ranked.copy()

        with tempfile.TemporaryDirectory() as tmp:
            config.output.report_format = "html"
            result = generate_report(ranked, all_scores, config, Path(tmp))
            html_content = result.read_text()
            self.assertNotIn("<script>", html_content)
            self.assertIn("&lt;script&gt;", html_content)


if __name__ == "__main__":
    unittest.main()
