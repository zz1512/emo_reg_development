"""Synthetic checks for leakage, grouping, and the locked signal patterns."""

import inspect
from pathlib import Path
import tempfile
import unittest

import numpy as np
import pandas as pd

from analysis.fig03.predictive_pls.bootstrap import (
    _load_task,
    confirmatory_flags,
    expand_index,
    site_draw,
)
from analysis.fig03.predictive_pls.make_pfc_definition import (
    build_definition,
    coverage_by_parcel,
    membership_from_coverage,
    vertex_areas,
)
from analysis.fig03.predictive_pls.models import (
    PRIMARY_COMPARISON,
    CovariateCoder,
    fit_fold,
    predict_fold,
    select_k,
)
from analysis.fig03.predictive_pls.run_pls import MEMBERSHIP_RULE, MEMBERSHIP_SOURCE, load_membership, nested_loso


def _batch(n, sites, seed):
    rng = np.random.default_rng(seed)
    raw = {
        "age": rng.uniform(7, 18, n),
        "gender": rng.integers(0, 2, n).astype(float),
        "parent_edu": rng.normal(12, 2, n),
        "meanFD": rng.uniform(0.05, 0.4, n),
    }
    return raw, np.asarray(sites)


class PipelineTests(unittest.TestCase):
    def test_vertex_area_and_coverage_rule(self):
        coordinates = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]], dtype=float)
        area = vertex_areas(coordinates, np.array([[0, 1, 2]]))
        self.assertTrue(np.allclose(area, np.full(3, 0.5 / 3)))
        labels = np.array([1, 1, 2, 2])
        inside = np.array([True, False, True, False])
        covered = coverage_by_parcel(labels, inside, np.array([3.0, 1.0, 0.2, 0.2]))
        self.assertAlmostEqual(covered[1][0] / covered[1][1], 0.75)
        info = pd.DataFrame({"parcel": range(1, 433), "name": [f"7Networks_LH_PFC_{i}" if i <= 200 else f"7Networks_RH_PFC_{i}" for i in range(1, 433)]})
        coverage = {parcel: (1.0, 1.0) if parcel <= 10 else (0.0, 1.0) for parcel in range(1, 401)}
        table = membership_from_coverage(coverage, info)
        self.assertEqual(set(table.loc[table["is_pfc"], "parcel"]), set(range(1, 11)))
        self.assertFalse(table.loc[table["parcel"].gt(400), "is_pfc"].any())
        self.assertTrue(table.loc[table["parcel"].eq(1), "boundary"].iloc[0] == False)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "membership.csv"
            table.to_csv(path, index=False)
            loaded = load_membership(path)
            self.assertEqual(int(loaded["is_pfc"].sum()), 10)
            table["definition_source"] = "name_list"
            table.to_csv(path, index=False)
            with self.assertRaises(Exception):
                load_membership(path)

    def test_missing_pfc_source_stops(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            with self.assertRaises(Exception) as caught:
                build_definition(root / "missing.dlabel.nii", root / "info.csv", root / "L.func.gii",
                                 root / "R.func.gii", root / "L.surf.gii", root / "R.surf.gii", root / "out")
        self.assertIn("absent", str(caught.exception))

    def test_knots_and_prediction_ignore_test_target(self):
        raw, _ = _batch(40, ["a"] * 40, 3)
        other, _ = _batch(40, ["b"] * 40, 9)
        first = CovariateCoder("success").fit(raw["age"], raw["gender"], raw["parent_edu"], raw["meanFD"])
        second = CovariateCoder("success").fit(other["age"], other["gender"], other["parent_edu"], other["meanFD"])
        self.assertFalse(np.allclose(first.spline_knots, second.spline_knots))
        self.assertEqual(first.transform(raw["age"], raw["gender"], raw["parent_edu"], raw["meanFD"]).shape[1], 8)
        x_train = np.random.default_rng(1).normal(size=(30, 3))
        x_test = np.random.default_rng(2).normal(size=(10, 3))
        y_train = x_train[:, 0] + 0.1 * raw["parent_edu"][:30]
        state = fit_fold(raw["age"][:30], raw["gender"][:30], raw["parent_edu"][:30], raw["meanFD"][:30],
                         y_train, "success", x_train, np.array([1, 2, 3]), components_through=1)
        prediction = predict_fold(state, raw["age"][30:], raw["gender"][30:], raw["parent_edu"][30:],
                                  raw["meanFD"][30:], x_test, 0)
        self.assertNotIn("y", inspect.signature(predict_fold).parameters)
        changed = x_test.copy()
        changed[0, 0] += 1e6
        moved = predict_fold(state, raw["age"][30:], raw["gender"][30:], raw["parent_edu"][30:],
                             raw["meanFD"][30:], changed, 1)
        original = predict_fold(state, raw["age"][30:], raw["gender"][30:], raw["parent_edu"][30:],
                                raw["meanFD"][30:], x_test, 1)
        self.assertFalse(np.allclose(original, moved))
        refit = fit_fold(raw["age"][:30], raw["gender"][:30], raw["parent_edu"][:30], raw["meanFD"][:30],
                         y_train, "success", x_train, np.array([1, 2, 3]), components_through=1)
        self.assertTrue(np.allclose(state.a, refit.a))
        self.assertEqual(prediction.shape, (10,))

    def test_one_pls_fit_matches_separate_component_fits(self):
        from sklearn.cross_decomposition import PLSRegression

        raw, _ = _batch(60, ["a"] * 60, 5)
        rng = np.random.default_rng(8)
        x_train = rng.normal(size=(40, 6))
        y_train = x_train[:, 0] - 0.4 * x_train[:, 1] + 0.2 * raw["age"][:40]
        state = fit_fold(raw["age"][:40], raw["gender"][:40], raw["parent_edu"][:40], raw["meanFD"][:40],
                         y_train, "success", x_train, np.arange(1, 7), components_through=4)
        design = state.coder.transform(raw["age"][:40], raw["gender"][:40], raw["parent_edu"][:40], raw["meanFD"][:40])
        residual_x = x_train - design @ state.B
        standardized = (residual_x[:, state.keep] - state.x_mean) / state.x_scale
        y_residual = y_train - design @ state.a
        for component in range(1, 5):
            separate = PLSRegression(n_components=component, scale=False).fit(standardized, y_residual)
            direct = np.asarray(separate.predict(standardized)).reshape(-1)
            centered = standardized - state.pls_x_mean
            fast = centered @ state.coef_by_k[component] + state.pls_intercept
            self.assertTrue(np.allclose(direct, fast, rtol=1e-8, atol=1e-8))

    def test_component_tie_prefers_smaller_k(self):
        self.assertEqual(select_k([0.50, 0.20, 0.20 + 1e-12, 0.90]), 1)

    def _run_signal(self, kind):
        raw, sites = _batch(80, np.repeat(["s1", "s2", "s3", "s4"], 20), 11)
        rng = np.random.default_rng(15)
        signal = rng.normal(size=80)
        noise = rng.normal(size=(80, 2))
        pfc = np.column_stack([signal, noise[:, :1]])
        other = np.column_stack([noise[:, 1:], rng.normal(size=80)])
        if kind == "non_pfc":
            y = 3 * other[:, 0]
            other = other.copy()
            other[:, 0] = other[:, 0]
        elif kind == "pfc":
            y = 3 * pfc[:, 0]
        else:
            y = 0.2 * raw["age"] + 0.3 * raw["parent_edu"] + 0.4 * raw["gender"] + 0.1 * raw["meanFD"]
        features = {
            "M1_PFC": (pfc, np.array([1, 2])),
            "M2_distributed432": (np.column_stack([pfc, other]), np.array([1, 2, 3, 4])),
            "M3_nonPFC": (other, np.array([3, 4])),
        }
        predictions, _, _ = nested_loso(raw, features, y, sites, "success", np.array([f"id{i}" for i in range(80)]))
        return predictions

    def test_signal_patterns(self):
        non_pfc = self._run_signal("non_pfc")
        r2 = lambda frame, model: 1 - np.sum((frame["y_true"] - frame[f"pred_{model}"]) ** 2) / np.sum((frame["y_true"] - frame["y_true"].mean()) ** 2)
        self.assertGreater(r2(non_pfc, "M2_distributed432") - r2(non_pfc, "M1_PFC"), 0.2)
        pfc = self._run_signal("pfc")
        self.assertGreater(r2(pfc, "M1_PFC"), 0.5)
        self.assertLess(r2(pfc, "M2_distributed432") - r2(pfc, "M1_PFC"), 0.15)
        null = self._run_signal("null")
        self.assertLess(abs(r2(null, "M2_distributed432") - r2(null, "M1_PFC")), 0.05)
        self.assertEqual(len(null), 80)
        self.assertFalse(null.duplicated("sub_id").any())

    def test_bootstrap_groups_and_confirmation_gate(self):
        sites = np.repeat(["a", "b", "c", "d"], 2)
        draw = site_draw(20261111, 4, np.array(["a", "b", "c", "d"]))
        self.assertTrue(np.array_equal(draw, site_draw(20261111, 4, np.array(["a", "b", "c", "d"]))))
        self.assertEqual(draw.size, 4)
        index = expand_index(sites, np.array(["a", "a", "b", "c"]))
        self.assertEqual(list(sites[index]).count("a"), 4)
        self.assertTrue(confirmatory_flags(1900, 2000, 0.02, "success", PRIMARY_COMPARISON)["supported"])
        self.assertFalse(confirmatory_flags(1900, 2000, 0.02, "age", PRIMARY_COMPARISON)["supported"])
        self.assertFalse(confirmatory_flags(30, 30, 0.02, "success", PRIMARY_COMPARISON)["inference_valid"])
        self.assertFalse(confirmatory_flags(1900, 2000, -0.01, "success", PRIMARY_COMPARISON)["supported"])

    def test_object_id_arrays_reload(self):
        with tempfile.TemporaryDirectory() as directory:
            frozen = Path(directory) / "frozen"
            frozen.mkdir()
            sites = np.array(["site-a", "site-b"], dtype=object)
            np.savez_compressed(frozen / "success_arrays.npz", age=np.array([8., 9.]), gender=np.array([0., 1.]),
                                parent_edu=np.array([12., 14.]), meanFD=np.array([0.2, 0.3]), y=np.array([1., 2.]),
                                sites=sites, sub_id=np.array(["s1", "s2"], dtype=object),
                                X=np.zeros((2, 432)), pfc_mask=np.array([True] + [False] * 431),
                                parcel_ids=np.arange(1, 433))
            loaded = _load_task(Path(directory), "success")
        self.assertEqual(loaded["sites"].dtype.kind, "U")
        self.assertEqual(list(loaded["site_order"]), ["site-a", "site-b"])

    def test_definition_constants(self):
        self.assertEqual(MEMBERSHIP_SOURCE, "Donahue2018_liberal_bilateral")
        self.assertEqual(MEMBERSHIP_RULE, "midthickness_area_fraction>=0.50")


if __name__ == "__main__":
    unittest.main()
