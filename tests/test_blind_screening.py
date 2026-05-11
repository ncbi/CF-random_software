#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Tests for BlindScreening.

Unit tests cover the pure static methods and each private helper in isolation.
Integration tests wire the full __init__ pipeline together with Foldseek and
MDAnalysis mocked out, so no external binaries or real PDB files are required.

Run with:
    pytest tests/test_blind_screening.py -v
"""

import csv
import shutil
import subprocess
from pathlib import Path
from typing import List
from unittest.mock import MagicMock, call, patch

import numpy as np
import pytest

from cf_random.utils.search_foldseek_cluster import (
    FOLDSEEK_BIT_SCORE_BUG_VALUE,
    KMEDOIDS_MIN_CLUSTER_SIZE,
    ZSCORE_OUTLIER_THRESHOLD,
    BlindScreening,
)

RNG = np.random.default_rng(0)


def _make_two_cluster_data(n_per_cluster: int = 30, n_components: int = 4) -> np.ndarray:
    """Return a well-separated two-cluster dataset for clustering tests."""
    c1 = RNG.normal(loc=[-5, 0, 0, 0], scale=0.3, size=(n_per_cluster, n_components))
    c2 = RNG.normal(loc=[5, 0, 0, 0], scale=0.3, size=(n_per_cluster, n_components))
    return np.vstack([c1, c2])


def _make_labels(n_per_cluster: int = 30) -> np.ndarray:
    """Return ground-truth labels matching _make_two_cluster_data."""
    return np.array([0] * n_per_cluster + [1] * n_per_cluster)


@pytest.fixture()
def tmp_blind_path(tmp_path: Path) -> Path:
    """Return a temporary directory that looks like a blind_path."""
    return tmp_path


class TestClusterStructures:
    def test_returns_array(self):
        X = _make_two_cluster_data()
        labels = BlindScreening.cluster_structures(X)
        assert isinstance(labels, np.ndarray)
        assert labels.shape == (len(X),)

    def test_finds_two_clusters_on_separable_data(self):
        X = _make_two_cluster_data(n_per_cluster=40)
        labels = BlindScreening.cluster_structures(X)
        # Noise points (-1) are allowed; meaningful labels should be 0 and 1
        meaningful = labels[labels >= 0]
        assert len(np.unique(meaningful)) == 2

    def test_all_noise_does_not_raise(self):
        # Single tight ball — HDBSCAN may label everything as one cluster or noise
        X = RNG.normal(size=(10, 4))
        labels = BlindScreening.cluster_structures(X)
        assert labels.shape == (10,)

    def test_labels_length_matches_input(self):
        for n in [10, 25, 60]:
            X = _make_two_cluster_data(n_per_cluster=n)
            labels = BlindScreening.cluster_structures(X)
            assert len(labels) == 2 * n


class TestKMedoids:
    def _make_input(self, n_per_cluster=20, n_components=4):
        X = _make_two_cluster_data(n_per_cluster, n_components)
        labels = _make_labels(n_per_cluster)
        return X, labels

    def test_returns_tuple(self):
        X, labels = self._make_input()
        result = BlindScreening.k_medoids(X, 0, labels)
        assert isinstance(result, tuple) and len(result) == 2

    def test_medoid_indices_within_bounds(self):
        X, labels = self._make_input()
        medoids, _ = BlindScreening.k_medoids(X, 0, labels, k=3)
        assert all(0 <= idx < len(X) for idx in medoids)

    def test_medoids_reduce_cost_vs_random(self):
        """Converged medoids should have cost <= initial random assignment.

        The implementation pads non-cluster points with 9999 so they can still
        be selected as medoid indices — checking cluster membership is therefore
        not a valid assertion. Cost reduction is the correct invariant.
        """
        X, labels = self._make_input(n_per_cluster=20)
        medoids_init, cost_init = BlindScreening.k_medoids(X, 0, labels, k=3, max_iter=0)
        medoids_conv, cost_conv = BlindScreening.k_medoids(X, 0, labels, k=3)
        assert cost_conv <= cost_init

    def test_cost_is_finite_for_large_cluster(self):
        X, labels = self._make_input(n_per_cluster=20)
        _, cost = BlindScreening.k_medoids(X, 0, labels, k=3)
        assert np.isfinite(cost)

    def test_small_cluster_returns_all_indices(self):
        """Cluster smaller than KMEDOIDS_MIN_CLUSTER_SIZE returns all members."""
        n_small = KMEDOIDS_MIN_CLUSTER_SIZE - 1
        X = RNG.normal(size=(10, 4))
        labels = np.array([0] * n_small + [1] * (10 - n_small))
        medoids, cost = BlindScreening.k_medoids(X, 0, labels, k=3)
        assert len(medoids) == n_small
        assert np.isnan(cost)

    def test_unknown_cluster_label_returns_empty(self):
        X, labels = self._make_input()
        medoids, cost = BlindScreening.k_medoids(X, 99, labels)
        assert len(medoids) == 0
        assert np.isnan(cost)

    def test_deterministic_with_fixed_seed(self):
        X, labels = self._make_input(n_per_cluster=20)
        m1, c1 = BlindScreening.k_medoids(X, 0, labels, k=3)
        m2, c2 = BlindScreening.k_medoids(X, 0, labels, k=3)
        np.testing.assert_array_equal(m1, m2)
        assert c1 == c2

    def test_k1_returns_single_medoid(self):
        X, labels = self._make_input(n_per_cluster=20)
        medoids, _ = BlindScreening.k_medoids(X, 0, labels, k=1)
        assert len(medoids) == 1


class TestBuildCorrelationMatrix:
    """Test _build_correlation_matrix via a minimal BlindScreening instance
    constructed without running __init__ (object.__new__)."""

    def _make_instance(self):
        obj = object.__new__(BlindScreening)
        return obj

    def _write_foldseek_file(self, path: Path, rows: List[tuple]) -> None:
        """Write a mock .foldseek TSV file."""
        with path.open("w") as fh:
            for row in rows:
                fh.write("\t".join(str(x) for x in row) + "\n")

    def test_basic_matrix_shape(self, tmp_path):
        labels = ["struct_a", "struct_b"]
        f1 = tmp_path / "a-self.foldseek"
        f2 = tmp_path / "b-self.foldseek"
        # query, target, alntmscore, qaln, taln, alnlen, evalue, bits
        self._write_foldseek_file(
            f1,
            [
                ("a", "struct_a", 0.9, "A", "A", 10, 0.001, 100),
                ("a", "struct_b", 0.5, "A", "B", 10, 0.01, 50),
            ],
        )
        self._write_foldseek_file(
            f2,
            [
                ("b", "struct_a", 0.5, "B", "A", 10, 0.01, 50),
                ("b", "struct_b", 0.9, "B", "B", 10, 0.001, 100),
            ],
        )
        obj = self._make_instance()
        mtx = obj._build_correlation_matrix([f1, f2], labels)
        assert mtx.shape == (2, 2)

    def test_known_values(self, tmp_path):
        labels = ["struct_a", "struct_b"]
        f1 = tmp_path / "a-self.foldseek"
        f2 = tmp_path / "b-self.foldseek"
        self._write_foldseek_file(
            f1,
            [
                ("a", "struct_a", 0.9, "A", "A", 10, 0.001, 200),
                ("a", "struct_b", 0.5, "A", "B", 10, 0.01, 75),
            ],
        )
        self._write_foldseek_file(
            f2,
            [
                ("b", "struct_a", 0.5, "B", "A", 10, 0.01, 75),
                ("b", "struct_b", 0.9, "B", "B", 10, 0.001, 300),
            ],
        )
        obj = self._make_instance()
        mtx = obj._build_correlation_matrix([f1, f2], labels)
        assert mtx[0, 0] == 200.0
        assert mtx[0, 1] == 75.0
        assert mtx[1, 0] == 75.0
        assert mtx[1, 1] == 300.0

    def test_foldseek_bug_value_replaced_with_zero(self, tmp_path):
        labels = ["struct_a"]
        f1 = tmp_path / "a-self.foldseek"
        self._write_foldseek_file(
            f1,
            [
                ("a", "struct_a", 0.9, "A", "A", 10, 0.001, FOLDSEEK_BIT_SCORE_BUG_VALUE),
            ],
        )
        obj = self._make_instance()
        mtx = obj._build_correlation_matrix([f1], labels)
        assert mtx[0, 0] == 0.0

    def test_missing_target_defaults_to_zero(self, tmp_path):
        labels = ["struct_a", "struct_b"]
        f1 = tmp_path / "a-self.foldseek"
        # Only struct_a present; struct_b missing → should default to 0
        self._write_foldseek_file(
            f1,
            [
                ("a", "struct_a", 0.9, "A", "A", 10, 0.001, 100),
            ],
        )
        obj = self._make_instance()
        mtx = obj._build_correlation_matrix([f1], labels)
        assert mtx[0, 1] == 0.0

    def test_short_lines_are_skipped(self, tmp_path):
        labels = ["struct_a"]
        f1 = tmp_path / "a-self.foldseek"
        with f1.open("w") as fh:
            fh.write("only\ttwo\tcolumns\n")
            fh.write("a\tstruct_a\t0.9\tA\tA\t10\t0.001\t100\n")
        obj = self._make_instance()
        mtx = obj._build_correlation_matrix([f1], labels)
        assert mtx[0, 0] == 100.0

    def test_non_integer_bit_score_skipped(self, tmp_path):
        labels = ["struct_a"]
        f1 = tmp_path / "a-self.foldseek"
        self._write_foldseek_file(
            f1,
            [
                ("a", "struct_a", 0.9, "A", "A", 10, 0.001, "nan"),
            ],
        )
        obj = self._make_instance()
        mtx = obj._build_correlation_matrix([f1], labels)
        assert mtx[0, 0] == 0.0


class TestFilterUnfolded:
    """Mock MDAnalysis/DSSP so we never need real PDB files."""

    def _make_foldseek_files(self, tmp_path: Path, n: int) -> List[Path]:
        files = []
        for i in range(n):
            ff = tmp_path / f"struct_{i:03d}-self.foldseek"
            ff.touch()
            # Create matching .pdb stub
            pdb = tmp_path / f"struct_{i:03d}.pdb"
            pdb.touch()
            files.append(ff)
        return sorted(files)

    def _dssp_mock(self, loop_counts: List[int]):
        """Return a side_effect list for DSSP that yields the given loop counts."""
        mocks = []
        for lc in loop_counts:
            dssp_array = np.array(["-"] * lc + ["H"] * 10 + ["E"] * 5)
            run_mock = MagicMock()
            run_mock.results.dssp = [dssp_array]
            instance = MagicMock()
            instance.run.return_value = run_mock
            mocks.append(instance)
        return mocks

    @patch("cf_random.utils.search_foldseek_cluster.mda.Universe")
    @patch("cf_random.utils.search_foldseek_cluster.DSSP")
    def test_no_outliers_returns_all(self, mock_dssp_cls, mock_universe, tmp_path):
        files = self._make_foldseek_files(tmp_path, 5)
        # All similar loop counts — no outlier
        loop_counts = [10, 11, 10, 9, 10]
        mock_dssp_cls.side_effect = self._dssp_mock(loop_counts)

        obj = object.__new__(BlindScreening)
        filtered, labels = obj._filter_unfolded(files)

        assert len(filtered) == 5
        assert len(labels) == 5

    @patch("cf_random.utils.search_foldseek_cluster.mda.Universe")
    @patch("cf_random.utils.search_foldseek_cluster.DSSP")
    def test_outlier_removed(self, mock_dssp_cls, mock_universe, tmp_path):
        # 5 samples with value 1000 only reaches z=2.0; use 20 normal + 1 extreme.
        n_normal = 20
        files = self._make_foldseek_files(tmp_path, n_normal + 1)
        loop_counts = [10] * n_normal + [10000]
        mock_dssp_cls.side_effect = self._dssp_mock(loop_counts)

        obj = object.__new__(BlindScreening)
        filtered, labels = obj._filter_unfolded(files)

        assert len(filtered) == n_normal
        assert len(labels) == n_normal

    @patch("cf_random.utils.search_foldseek_cluster.mda.Universe")
    @patch("cf_random.utils.search_foldseek_cluster.DSSP")
    def test_labels_match_slash_dash_rule(self, mock_dssp_cls, mock_universe, tmp_path):
        """Labels must equal str(pdb_path).replace('/', '-')[17:].replace('.pdb', '')."""
        files = self._make_foldseek_files(tmp_path, 3)
        mock_dssp_cls.side_effect = self._dssp_mock([10, 10, 10])

        obj = object.__new__(BlindScreening)
        filtered, labels = obj._filter_unfolded(files)

        for ff, label in zip(filtered, labels):
            pdb_path = str(ff).replace("-self.foldseek", ".pdb")
            expected = pdb_path.replace("/", "-")[17:].replace(".pdb", "")
            assert label == expected

    @patch("cf_random.utils.search_foldseek_cluster.mda.Universe")
    @patch("cf_random.utils.search_foldseek_cluster.DSSP")
    def test_output_is_sorted(self, mock_dssp_cls, mock_universe, tmp_path):
        files = self._make_foldseek_files(tmp_path, 5)
        mock_dssp_cls.side_effect = self._dssp_mock([10] * 5)

        obj = object.__new__(BlindScreening)
        filtered, _ = obj._filter_unfolded(files)

        assert filtered == sorted(filtered)


class TestStagePdbFiles:
    def _make_pdb_tree(self, base: Path, names: List[str]) -> List[Path]:
        paths = []
        for name in names:
            p = base / name
            p.parent.mkdir(parents=True, exist_ok=True)
            p.write_text("ATOM  dummy\n")
            paths.append(p)
        return paths

    def test_staging_directory_created(self, tmp_path):
        self._make_pdb_tree(tmp_path, ["subdir/a.pdb", "subdir/b.pdb"])
        obj = object.__new__(BlindScreening)
        obj.blind_path = tmp_path
        obj._stage_pdb_files()
        assert (tmp_path / "pdbs_for_db").is_dir()

    def test_all_pdbs_staged(self, tmp_path):
        self._make_pdb_tree(tmp_path, ["subdir/a.pdb", "subdir/b.pdb"])
        obj = object.__new__(BlindScreening)
        obj.blind_path = tmp_path
        obj._stage_pdb_files()
        staged = list((tmp_path / "pdbs_for_db").iterdir())
        assert len(staged) == 2

    def test_pdbs_inside_db_directory_excluded(self, tmp_path):
        self._make_pdb_tree(tmp_path, ["subdir/a.pdb"])
        obj = object.__new__(BlindScreening)
        obj.blind_path = tmp_path
        obj._stage_pdb_files()
        # Manually add a .pdb inside pdbs_for_db and re-stage
        (tmp_path / "pdbs_for_db" / "extra.pdb").write_text("ATOM  extra\n")
        obj2 = object.__new__(BlindScreening)
        obj2.blind_path = tmp_path
        obj2._stage_pdb_files()
        staged = [
            f
            for f in (tmp_path / "pdbs_for_db").iterdir()
            if f.suffix == ".pdb" and f.name != "extra.pdb"
        ]
        assert len(staged) == 1

    def test_label_map_populated(self, tmp_path):
        self._make_pdb_tree(tmp_path, ["subdir/a.pdb", "subdir/b.pdb"])
        obj = object.__new__(BlindScreening)
        obj.blind_path = tmp_path
        obj._stage_pdb_files()
        assert len(obj.pdb_label_map) == 2

    def test_label_map_values_match_slash_dash_rule(self, tmp_path):
        self._make_pdb_tree(tmp_path, ["subdir/a.pdb"])
        obj = object.__new__(BlindScreening)
        obj.blind_path = tmp_path
        obj._stage_pdb_files()
        for src, label in obj.pdb_label_map.items():
            expected = str(src).replace("/", "-")[17:].replace(".pdb", "")
            assert label == expected

    def test_raises_if_no_pdbs(self, tmp_path):
        obj = object.__new__(BlindScreening)
        obj.blind_path = tmp_path
        with pytest.raises(FileNotFoundError):
            obj._stage_pdb_files()

    def test_existing_staged_file_not_overwritten(self, tmp_path):
        self._make_pdb_tree(tmp_path, ["subdir/a.pdb"])
        obj = object.__new__(BlindScreening)
        obj.blind_path = tmp_path
        obj._stage_pdb_files()
        # Corrupt the staged file
        staged = list((tmp_path / "pdbs_for_db").glob("*.pdb"))[0]
        staged.write_text("CORRUPTED")
        # Re-stage should not overwrite
        obj2 = object.__new__(BlindScreening)
        obj2.blind_path = tmp_path
        obj2._stage_pdb_files()
        assert staged.read_text() == "CORRUPTED"


class TestSaveClusterPlot:
    def test_file_created(self, tmp_path):
        obj = object.__new__(BlindScreening)
        obj.blind_path = tmp_path
        obj.pdb1_name = "myprotein"
        coords = RNG.normal(size=(20, 4))
        labels = np.array([0] * 10 + [1] * 10)
        obj._save_cluster_plot(coords, labels)
        assert (tmp_path / "myprotein-cluster.png").exists()


class TestSaveStructuresOfInterest:
    def test_csv_created_with_header(self, tmp_path):
        obj = object.__new__(BlindScreening)
        obj.blind_path = tmp_path
        obj.pdb1_name = "myprotein"
        files_of_interest = [(Path("a.foldseek"), 0), (Path("b.foldseek"), 1)]
        pca_of_interest = [np.array([1.0, 2.0, 0, 0]), np.array([3.0, 4.0, 0, 0])]
        obj._save_structures_of_interest(files_of_interest, pca_of_interest)
        out = tmp_path / "myprotein-structures_of_interest.csv"
        assert out.exists()
        with out.open() as fh:
            reader = csv.reader(fh)
            header = next(reader)
        assert header == ["group", "file", "pca_1", "pca_2"]

    def test_csv_row_count(self, tmp_path):
        obj = object.__new__(BlindScreening)
        obj.blind_path = tmp_path
        obj.pdb1_name = "myprotein"
        n = 5
        files_of_interest = [(Path(f"f{i}.foldseek"), i % 2) for i in range(n)]
        pca_of_interest = [np.array([float(i), 0.0, 0, 0]) for i in range(n)]
        obj._save_structures_of_interest(files_of_interest, pca_of_interest)
        out = tmp_path / "myprotein-structures_of_interest.csv"
        with out.open() as fh:
            rows = list(csv.reader(fh))
        assert len(rows) == n + 1  # header + n data rows


class TestSaveAllStructures:
    def test_csv_created(self, tmp_path):
        obj = object.__new__(BlindScreening)
        obj.blind_path = tmp_path
        files = [Path(f"s{i}-self.foldseek") for i in range(4)]
        labels = np.array([0, 0, 1, 1])
        coords = RNG.normal(size=(4, 4))
        obj._save_all_structures(files, labels, coords)
        assert (tmp_path / "structures_all.csv").exists()

    def test_row_count(self, tmp_path):
        obj = object.__new__(BlindScreening)
        obj.blind_path = tmp_path
        n = 6
        files = [Path(f"s{i}-self.foldseek") for i in range(n)]
        labels = np.zeros(n, dtype=int)
        coords = RNG.normal(size=(n, 4))
        obj._save_all_structures(files, labels, coords)
        with (tmp_path / "structures_all.csv").open() as fh:
            rows = list(csv.reader(fh))
        assert len(rows) == n + 1


class TestBlindScreeningIntegration:
    """Run the complete __init__ pipeline without Foldseek, MDAnalysis, or PyMOL."""

    N_STRUCTS = 11  # enough for HDBSCAN to find structure

    @pytest.fixture()
    def pipeline_dir(self, tmp_path: Path) -> Path:
        """Build a minimal directory tree with fake PDB and .foldseek files."""
        for i in range(self.N_STRUCTS):
            sub = tmp_path / f"run_{i:03d}"
            sub.mkdir()
            pdb = sub / f"struct_{i:03d}_rank_001_model.pdb"
            pdb.write_text(f"ATOM  {i}\n")

        # Pre-create .foldseek result files so _run_foldseek_searches skips them
        # Labels must match the slash->dash[17:] rule for this tmp_path
        pdb_files = sorted(tmp_path.rglob("*.pdb"))
        labels = [str(p).replace("/", "-")[17:].replace(".pdb", "") for p in pdb_files]

        for i, pdb in enumerate(pdb_files):
            ff = pdb.with_name(pdb.stem + "-self.foldseek")
            with ff.open("w") as fh:
                for j, label in enumerate(labels):
                    # Give high scores on diagonal, lower elsewhere
                    score = 500 if i == j else (200 if abs(i - j) <= 2 else 50)
                    fh.write(f"query\t{label}\t0.9\tA\tA\t10\t0.001\t{score}\n")

        return tmp_path

    def _make_dssp_side_effect(self, n: int):
        """Return DSSP constructor side effects for n structures (no outliers)."""
        mocks = []
        for _ in range(n):
            dssp_array = np.array(["-"] * 10 + ["H"] * 20 + ["E"] * 10)
            run_mock = MagicMock()
            run_mock.results.dssp = [dssp_array]
            instance = MagicMock()
            instance.run.return_value = run_mock
            mocks.append(instance)
        return mocks

    @patch("cf_random.utils.search_foldseek_cluster.pymol", None)
    @patch("cf_random.utils.search_foldseek_cluster.DSSP")
    @patch("cf_random.utils.search_foldseek_cluster.mda.Universe")
    @patch("cf_random.utils.search_foldseek_cluster.subprocess.run")
    def test_pipeline_completes(self, mock_run, mock_universe, mock_dssp_cls, pipeline_dir):
        mock_run.return_value = MagicMock(returncode=0, stderr="")
        mock_dssp_cls.side_effect = self._make_dssp_side_effect(self.N_STRUCTS)

        bs = BlindScreening("testprot", str(pipeline_dir))

        assert isinstance(bs, BlindScreening)

    @patch("cf_random.utils.search_foldseek_cluster.pymol", None)
    @patch("cf_random.utils.search_foldseek_cluster.DSSP")
    @patch("cf_random.utils.search_foldseek_cluster.mda.Universe")
    @patch("cf_random.utils.search_foldseek_cluster.subprocess.run")
    def test_output_files_created(self, mock_run, mock_universe, mock_dssp_cls, pipeline_dir):
        mock_run.return_value = MagicMock(returncode=0, stderr="")
        mock_dssp_cls.side_effect = self._make_dssp_side_effect(self.N_STRUCTS)

        BlindScreening("testprot", str(pipeline_dir))

        assert (pipeline_dir / "testprot-cluster.png").exists()
        assert (pipeline_dir / "testprot-structures_of_interest.csv").exists()
        assert (pipeline_dir / "structures_all.csv").exists()

    @patch("cf_random.utils.search_foldseek_cluster.pymol", None)
    @patch("cf_random.utils.search_foldseek_cluster.DSSP")
    @patch("cf_random.utils.search_foldseek_cluster.mda.Universe")
    @patch("cf_random.utils.search_foldseek_cluster.subprocess.run")
    def test_structures_of_interest_csv_has_header(
        self, mock_run, mock_universe, mock_dssp_cls, pipeline_dir
    ):
        mock_run.return_value = MagicMock(returncode=0, stderr="")
        mock_dssp_cls.side_effect = self._make_dssp_side_effect(self.N_STRUCTS)

        BlindScreening("testprot", str(pipeline_dir))

        with (pipeline_dir / "testprot-structures_of_interest.csv").open() as fh:
            header = next(csv.reader(fh))
        assert header == ["group", "file", "pca_1", "pca_2"]

    @patch("cf_random.utils.search_foldseek_cluster.pymol", None)
    @patch("cf_random.utils.search_foldseek_cluster.DSSP")
    @patch("cf_random.utils.search_foldseek_cluster.mda.Universe")
    @patch("cf_random.utils.search_foldseek_cluster.subprocess.run")
    def test_all_structures_csv_row_count(
        self, mock_run, mock_universe, mock_dssp_cls, pipeline_dir
    ):
        mock_run.return_value = MagicMock(returncode=0, stderr="")
        mock_dssp_cls.side_effect = self._make_dssp_side_effect(self.N_STRUCTS)

        BlindScreening("testprot", str(pipeline_dir))

        with (pipeline_dir / "structures_all.csv").open() as fh:
            rows = list(csv.reader(fh))
        # header + one row per structure
        assert len(rows) == self.N_STRUCTS + 1

    @patch("cf_random.utils.search_foldseek_cluster.pymol", None)
    @patch("cf_random.utils.search_foldseek_cluster.DSSP")
    @patch("cf_random.utils.search_foldseek_cluster.mda.Universe")
    @patch("cf_random.utils.search_foldseek_cluster.subprocess.run")
    def test_outlier_removed_from_all_structures_csv(
        self, mock_run, mock_universe, mock_dssp_cls, pipeline_dir
    ):
        """One structure with extreme loop count should be excluded from outputs."""
        mock_run.return_value = MagicMock(returncode=0, stderr="")
        # Last structure is an outlier
        normal = [10] * (self.N_STRUCTS - 1)
        outlier_dssp = self._make_dssp_side_effect(self.N_STRUCTS - 1)
        # Build outlier mock with extreme loop count
        # 5000 loops against ~9 normals gives z < 3; use 10000 to reliably exceed threshold.
        extreme_array = np.array(["-"] * 10000 + ["H"] * 5)
        run_mock = MagicMock()
        run_mock.results.dssp = [extreme_array]
        outlier_instance = MagicMock()
        outlier_instance.run.return_value = run_mock
        mock_dssp_cls.side_effect = outlier_dssp + [outlier_instance]

        BlindScreening("testprot", str(pipeline_dir))

        with (pipeline_dir / "structures_all.csv").open() as fh:
            rows = list(csv.reader(fh))
        assert len(rows) == self.N_STRUCTS  # header + (N-1) structures

    def test_raises_on_missing_path(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            BlindScreening("test", str(tmp_path / "does_not_exist"))

    @patch("cf_random.utils.search_foldseek_cluster.pymol", None)
    @patch("cf_random.utils.search_foldseek_cluster.DSSP")
    @patch("cf_random.utils.search_foldseek_cluster.mda.Universe")
    @patch("cf_random.utils.search_foldseek_cluster.subprocess.run")
    def test_foldseek_createdb_called_once(
        self, mock_run, mock_universe, mock_dssp_cls, pipeline_dir
    ):
        mock_run.return_value = MagicMock(returncode=0, stderr="")
        mock_dssp_cls.side_effect = self._make_dssp_side_effect(self.N_STRUCTS)

        BlindScreening("testprot", str(pipeline_dir))

        createdb_calls = [c for c in mock_run.call_args_list if c.args and "createdb" in c.args[0]]
        # .foldseek files already exist so easy-search is skipped;
        # createdb should be called exactly once
        assert len(createdb_calls) == 1

    @patch("cf_random.utils.search_foldseek_cluster.pymol", None)
    @patch("cf_random.utils.search_foldseek_cluster.DSSP")
    @patch("cf_random.utils.search_foldseek_cluster.mda.Universe")
    @patch("cf_random.utils.search_foldseek_cluster.subprocess.run")
    def test_existing_db_skips_createdb(self, mock_run, mock_universe, mock_dssp_cls, pipeline_dir):
        mock_run.return_value = MagicMock(returncode=0, stderr="")
        mock_dssp_cls.side_effect = self._make_dssp_side_effect(self.N_STRUCTS)

        # Pre-create the DB file so the code skips createdb
        db_dir = pipeline_dir / "pdbs_for_db"
        db_dir.mkdir(exist_ok=True)
        (db_dir / "DB").touch()

        BlindScreening("testprot", str(pipeline_dir))

        createdb_calls = [c for c in mock_run.call_args_list if c.args and "createdb" in c.args[0]]
        assert len(createdb_calls) == 0
