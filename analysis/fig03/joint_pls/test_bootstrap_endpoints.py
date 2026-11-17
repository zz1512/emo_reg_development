"""Rebuilding M0 must recover the original absolute bootstrap model scores."""
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

from .bootstrap_endpoints import baseline_r2, recover, read_intervals
from .io import MODELS, SCHEMES
from .pipeline import nested_loso
from .resampling import bootstrap_index
from .run import mock
from .test_joint_pls import fixture


class BootstrapEndpointTests(unittest.TestCase):
    def test_first_reference_rebuild_matches_full_nested_targets(self):
        data = fixture()
        draw = data.take(bootstrap_index(data.sites, 88, 2))
        with threadpool_limits(limits=1):
            baseline = baseline_r2(draw)
            expected = nested_loso(draw, 3, details=False)["performance"].set_index(["scheme", "model"])
        for scheme in SCHEMES:
            self.assertAlmostEqual(baseline[scheme], expected.at[(scheme, MODELS[0]), "r2"], places=12)

    def test_archived_draws_recover_all_models_and_preserve_pairing(self):
        membership = Path(__file__).resolve().parents[3] / "outputs/fig03/predictive_pls_432_FD05_meta_pfc/pfc_definition/pfc_parcel_membership.csv"
        with tempfile.TemporaryDirectory() as folder:
            root = mock(Path(folder)/"mock", membership, n=96, permutations=0, bootstraps=2, docx=False)
            result = recover(root)
            intervals = read_intervals(root)
            self.assertEqual(len(intervals), 8)
            self.assertTrue(intervals.complete.all())
            scores = pd.read_csv(result.parent / "paired_bootstrap_model_r2.csv")
            self.assertEqual(len(scores), 16)
            for replicate, part in scores.groupby("replicate"):
                saved = pd.read_csv(root / "bootstrap/independent" / f"rep_{replicate:06d}.csv")
                p = part.set_index(["scheme", "model"])
                for row in saved.itertuples():
                    pair = {
                        "pfc_minus_covariates": (MODELS[1], MODELS[0]),
                        "distributed_minus_covariates": (MODELS[2], MODELS[0]),
                        "nonpfc_minus_covariates": (MODELS[3], MODELS[0]),
                        "distributed_minus_pfc": (MODELS[2], MODELS[1]),
                        "distributed_minus_nonpfc": (MODELS[2], MODELS[3]),
                    }[row.comparison]
                    self.assertAlmostEqual(p.at[(row.scheme, pair[0]), "r2"] - p.at[(row.scheme, pair[1]), "r2"], row.delta_r2, places=12)
            # Reuse is deterministic; an unrelated edited interval file is rejected.
            original = result.read_bytes()
            recover(root)
            self.assertEqual(original, result.read_bytes())
            table = pd.read_csv(result)
            table.loc[0, "ci_low"] -= .1
            table.to_csv(result, index=False)
            with self.assertRaisesRegex(ValueError, "no longer match"):
                read_intervals(root)


if __name__ == "__main__":
    unittest.main()
