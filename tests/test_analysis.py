"""Tests for analysis module (filtering, ranking). Requires pandas."""

from __future__ import annotations

import unittest

try:
    import pandas as pd
    HAS_PANDAS = True
except ImportError:
    HAS_PANDAS = False


@unittest.skipUnless(HAS_PANDAS, "pandas not installed")
class TestFiltering(unittest.TestCase):
    def setUp(self):
        from harness.analysis.filter import apply_filters
        self.apply_filters = apply_filters
        self.scores = pd.DataFrame({
            "tag": ["d001", "d002", "d003", "d004", "d005"],
            "pae": [5.0, 8.0, 12.0, 3.0, 9.5],
            "rmsd": [1.0, 1.5, 3.0, 0.5, 1.8],
            "ddg": [-25.0, -15.0, -30.0, -22.0, -18.0],
        })

    def test_default_thresholds(self):
        filtered = self.apply_filters(self.scores)
        self.assertEqual(len(filtered), 2)
        self.assertEqual(set(filtered["tag"]), {"d001", "d004"})

    def test_no_ddg_filter(self):
        filtered = self.apply_filters(self.scores, ddg_threshold=None)
        self.assertEqual(len(filtered), 4)

    def test_empty_result(self):
        filtered = self.apply_filters(self.scores, pae_threshold=1.0)
        self.assertEqual(len(filtered), 0)


@unittest.skipUnless(HAS_PANDAS, "pandas not installed")
class TestRanking(unittest.TestCase):
    def setUp(self):
        from harness.analysis.rank import rank_candidates
        self.rank_candidates = rank_candidates

    def test_rank_order(self):
        scores = pd.DataFrame({
            "tag": ["d001", "d002", "d004"],
            "pae": [5.0, 8.0, 3.0],
            "rmsd": [1.0, 1.5, 0.5],
        })
        ranked = self.rank_candidates(scores)
        self.assertIn("rank", ranked.columns)
        self.assertIn("composite_score", ranked.columns)
        self.assertEqual(ranked.iloc[0]["rank"], 1)
        self.assertEqual(ranked.iloc[0]["tag"], "d004")

    def test_single_candidate(self):
        scores = pd.DataFrame({
            "tag": ["only_one"],
            "pae": [5.0],
            "rmsd": [1.0],
        })
        ranked = self.rank_candidates(scores)
        self.assertEqual(len(ranked), 1)
        self.assertEqual(ranked.iloc[0]["composite_score"], 0.0)


if __name__ == "__main__":
    unittest.main()
