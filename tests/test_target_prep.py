"""Tests for target preparation module (mocked I/O)."""

from __future__ import annotations

import unittest
from pathlib import Path
from unittest.mock import patch

from harness.config.defaults import BUILTIN_FRAMEWORKS
from harness.target_prep.convert_framework import prepare_framework
from harness.target_prep.fetch_pdb import fetch_pdb
from harness.target_prep.validate import ValidationResult, _check_patch_contiguity


class TestFetchPdb(unittest.TestCase):
    def test_already_downloaded(self):
        import tempfile
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            pdb_path = tmp / "5N2C.pdb"
            pdb_path.write_text("ATOM mock")
            result = fetch_pdb("5N2C", tmp)
            self.assertEqual(result, pdb_path)

    @patch("harness.target_prep.fetch_pdb.urlretrieve")
    def test_download(self, mock_retrieve):
        import tempfile

        def fake_download(url, path):
            Path(path).write_text("ATOM mock download")

        mock_retrieve.side_effect = fake_download
        with tempfile.TemporaryDirectory() as tmp:
            result = fetch_pdb("1ABC", Path(tmp))
            self.assertTrue(result.exists())
            self.assertEqual(result.name, "1ABC.pdb")
            mock_retrieve.assert_called_once()

    def test_uppercase_normalization(self):
        import tempfile
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            pdb_path = tmp / "5N2C.pdb"
            pdb_path.write_text("ATOM mock")
            result = fetch_pdb("  5n2c  ", tmp)
            self.assertEqual(result, pdb_path)


class TestPrepareFramework(unittest.TestCase):
    def test_unknown_builtin_raises(self):
        import tempfile
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            with self.assertRaises(ValueError):
                prepare_framework("builtin:Nonexistent", "vhh", tmp, tmp / "out.pdb")

    def test_builtin_file_missing_raises(self):
        import tempfile
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(FileNotFoundError):
                prepare_framework("builtin:NbBCII10", "vhh", Path(tmp), Path(tmp) / "out.pdb")

    def test_builtin_copy(self):
        import tempfile
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            rfab_root = tmp / "RFAntibody"
            meta = BUILTIN_FRAMEWORKS["NbBCII10"]
            src = rfab_root / meta["relative_path"]
            src.parent.mkdir(parents=True)
            src.write_text("ATOM mock framework")

            out = tmp / "output" / "framework.pdb"
            result = prepare_framework("builtin:NbBCII10", "vhh", rfab_root, out)
            self.assertTrue(result.exists())
            self.assertEqual(result.read_text(), "ATOM mock framework")

    def test_custom_file_missing_raises(self):
        import tempfile
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(FileNotFoundError):
                prepare_framework("/nonexistent/path.pdb", "vhh", Path(tmp), Path(tmp) / "out.pdb")


class TestPatchContiguity(unittest.TestCase):
    def test_small_gap_no_warning(self):
        result = ValidationResult()
        _check_patch_contiguity({}, [10, 12, 15], result)
        self.assertEqual(result.warnings, [])

    def test_large_gap_warns(self):
        result = ValidationResult()
        _check_patch_contiguity({}, [10, 50, 100], result)
        self.assertTrue(any("large sequence gap" in w for w in result.warnings))

    def test_single_residue_no_warning(self):
        result = ValidationResult()
        _check_patch_contiguity({}, [10], result)
        self.assertEqual(result.warnings, [])


if __name__ == "__main__":
    unittest.main()
