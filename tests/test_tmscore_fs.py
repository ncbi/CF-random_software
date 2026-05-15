#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Tests for monomer and multimer fold-switching TM-score pipelines.

Unit tests cover TMScoreFS (monomer), TMScoreFSMulti (multimer),
TMScore for both paths, and TMScoreCalAllVarFS._evaluate_monomer /
_evaluate_multimer in isolation.

Run with:
    pytest tests/test_tmscore_fs.py -v
"""

import glob
from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np
import pytest

from cf_random.utils.convert_multi_single import ConvertM2S
from cf_random.analysis.cal_tmscore_fs_flmsa import TMScoreFS
from cf_random.analysis.cal_tmscore_fs_multimer import TMScoreFSMulti
from cf_random.analysis.tmscore_all_var_fs import TMScoreCalAllVarFS, MSA_MULTIPLIERS
from cf_random.analysis.base import BaseTMScore

PDB1_NAME = "2oug_C"
PDB2_NAME = "6c6s_D"
MONOMER_MODEL_TYPE = "alphafold2_ptm"
MULTIMER_MODEL_TYPE = "alphafold2_multimer_v3"
NUM_SEEDS = 5
NUM_MODELS = 5
NUM_PREDICTIONS = NUM_SEEDS * NUM_MODELS  # 25

RNG = np.random.default_rng(42)

RESNAMES = ["ALA", "GLY", "VAL", "LEU", "ILE", "PRO", "PHE", "TRP", "MET", "SER"]


def _write_minimal_pdb(path: Path, n_residues: int = 10, n_chains: int = 1) -> None:
    """Write a minimal PDB with CA atoms and TER records per chain."""
    with path.open("w") as fh:
        atom_num = 1
        res_num = 1
        for chain_idx in range(n_chains):
            chain_id = chr(ord("A") + chain_idx)
            for i in range(n_residues):
                resname = RESNAMES[i % len(RESNAMES)]
                x, y, z = RNG.uniform(-10, 10, 3)
                fh.write(
                    f"ATOM  {atom_num:5d}  CA  {resname} {chain_id}"
                    f"{res_num:4d}    "
                    f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00 10.00           C\n"
                )
                atom_num += 1
                res_num += 1
            fh.write("TER\n")
        fh.write("END\n")


def _make_monomer_prediction_dir(
    base: Path,
    name: str,
    n_models: int = NUM_PREDICTIONS,
    n_residues: int = 10,
) -> Path:
    """Create a fake ColabFold monomer output directory.

    Monomer filenames: <model>_unrelaxed_rank_NNN_<model_type>.pdb
    No chain prefix, no TER between chains.
    """
    pred_dir = base / name
    pred_dir.mkdir(parents=True, exist_ok=True)
    for i in range(1, n_models + 1):
        pdb = pred_dir / f"model_{i}_unrelaxed_rank_{i:03d}_alphafold2_ptm.pdb"
        _write_minimal_pdb(pdb, n_residues=n_residues, n_chains=1)
    return pred_dir


def _make_multimer_prediction_dir(
    base: Path,
    name: str,
    n_models: int = NUM_PREDICTIONS,
    n_chains: int = 2,
    n_residues: int = 10,
) -> Path:
    """Create a fake ColabFold multimer output directory.

    Multimer filenames: 0_unrelaxed_rank_NNN_alphafold2_multimer_v3_model_1_seed_000.pdb
    """
    pred_dir = base / name
    pred_dir.mkdir(parents=True, exist_ok=True)
    for i in range(1, n_models + 1):
        pdb = pred_dir / (f"0_unrelaxed_rank_{i:03d}_alphafold2_multimer_v3_model_1_seed_000.pdb")
        _write_minimal_pdb(pdb, n_residues=n_residues, n_chains=n_chains)
    return pred_dir


def _make_converted_multimer_dir(
    base: Path,
    name: str,
    n_models: int = NUM_PREDICTIONS,
    n_residues: int = 10,
) -> Path:
    """Create a multimer prediction dir and run ConvertM2S on it."""
    pred_dir = _make_multimer_prediction_dir(base, name, n_models=n_models, n_residues=n_residues)
    pdb2_file = base / f"{PDB2_NAME}.pdb"
    if not pdb2_file.exists():
        _write_minimal_pdb(pdb2_file, n_residues=n_residues, n_chains=1)
    ConvertM2S(str(pred_dir), PDB1_NAME, PDB2_NAME)
    return pred_dir


def _make_range_file(path: Path, fs_range: str = "1-5") -> None:
    """Write a minimal range_fs_pairs_all.txt."""
    with path.open("w") as fh:
        fh.write("# pdb1,pdb2,pred1,pred2,m1,m2\n")
        fh.write(f"{PDB1_NAME},{PDB2_NAME},{fs_range},{fs_range},{fs_range},{fs_range}\n")


def _good_tm_result(value: float = 0.6) -> MagicMock:
    result = MagicMock()
    result.tm_norm_chain1 = value
    return result


class TestTMScoreFSMonomer:
    def test_get_coords_correct_range(self, tmp_path):
        pdb = tmp_path / "test.pdb"
        _write_minimal_pdb(pdb, n_residues=10, n_chains=1)
        scorer = object.__new__(TMScoreFS)
        coords, seq = scorer.get_coords(pdb, "1-5")
        assert len(coords) == 5
        assert len(seq) == 5

    def test_get_coords_full_range(self, tmp_path):
        pdb = tmp_path / "test.pdb"
        _write_minimal_pdb(pdb, n_residues=10, n_chains=1)
        scorer = object.__new__(TMScoreFS)
        coords, seq = scorer.get_coords(pdb, "1-10")
        assert len(coords) == 10
        assert len(seq) == 10

    def test_get_tmscore_returns_zeros_on_empty_dir(self, tmp_path):
        empty = tmp_path / "empty"
        empty.mkdir()
        scorer = object.__new__(TMScoreFS)
        coords = RNG.uniform(size=(5, 3))
        scores = scorer.get_tmscore(coords, "AGVLI", empty, "1-5")
        assert scores == [0.0, 0.0, 0.0, 0.0, 0.0]

    def test_get_tmscore_uses_unrelaxed_glob(self, tmp_path):
        pred_dir = _make_monomer_prediction_dir(tmp_path, "pred", n_models=3)
        matched = glob.glob(str(pred_dir) + "/*_unrelaxed*pdb")
        assert len(matched) == 3
        assert all("rmTER" not in f for f in matched)

    @patch("cf_random.analysis.cal_tmscore_fs_flmsa.tm_align")
    def test_get_tmscore_picks_best_orientation(self, mock_align, tmp_path):
        """Forward scores > reverse → forward selected."""
        pred_dir = _make_monomer_prediction_dir(tmp_path, "pred", n_models=2, n_residues=5)
        scorer = object.__new__(TMScoreFS)

        fwd = MagicMock()
        fwd.tm_norm_chain1 = 0.8
        rev = MagicMock()
        rev.tm_norm_chain1 = 0.3
        mock_align.side_effect = [fwd, rev, fwd, rev]

        coords = RNG.uniform(size=(5, 3))
        scores = scorer.get_tmscore(coords, "AGVLI", pred_dir, "1-5")
        assert all(s == 0.8 for s in scores)

    @patch("cf_random.analysis.cal_tmscore_fs_flmsa.tm_align")
    def test_get_tmscore_picks_reverse_when_higher(self, mock_align, tmp_path):
        pred_dir = _make_monomer_prediction_dir(tmp_path, "pred", n_models=2, n_residues=5)
        scorer = object.__new__(TMScoreFS)

        fwd = MagicMock()
        fwd.tm_norm_chain1 = 0.3
        rev = MagicMock()
        rev.tm_norm_chain1 = 0.8
        mock_align.side_effect = [fwd, rev, fwd, rev]

        coords = RNG.uniform(size=(5, 3))
        scores = scorer.get_tmscore(coords, "AGVLI", pred_dir, "1-5")
        assert all(s == 0.8 for s in scores)

    @patch("cf_random.analysis.cal_tmscore_fs_flmsa.tm_align")
    def test_run_for_models_shape(self, mock_align, tmp_path, monkeypatch):
        """run_for_models does pdb1 loop then pdb2 loop → shape (2, n_models)."""
        monkeypatch.chdir(tmp_path)
        _make_range_file(tmp_path / "range_fs_pairs_all.txt")
        pred_dir = _make_monomer_prediction_dir(tmp_path, "pred", n_models=5, n_residues=10)
        pdb1 = tmp_path / f"{PDB1_NAME}.pdb"
        pdb2 = tmp_path / f"{PDB2_NAME}.pdb"
        _write_minimal_pdb(pdb1, n_residues=10)
        _write_minimal_pdb(pdb2, n_residues=10)

        mock_align.return_value = _good_tm_result(0.6)

        scorer = object.__new__(TMScoreFS)
        scorer.run_for_models(pdb1, pdb2, str(pred_dir), "1-5", "1-5", "1-5")
        assert scorer.tmscores_fs.shape == (2, 5)

    @patch("cf_random.analysis.cal_tmscore_fs_flmsa.tm_align")
    def test_full_init_completes(self, mock_align, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        _make_range_file(tmp_path / "range_fs_pairs_all.txt")
        pred_root = tmp_path / "predictions_all" / PDB1_NAME
        pred_root.mkdir(parents=True)
        _make_monomer_prediction_dir(
            pred_root, f"{PDB1_NAME}_predicted_models_full_rand_0", n_residues=10
        )
        pdb1 = tmp_path / f"{PDB1_NAME}.pdb"
        pdb2 = tmp_path / f"{PDB2_NAME}.pdb"
        _write_minimal_pdb(pdb1, n_residues=10)
        _write_minimal_pdb(pdb2, n_residues=10)
        mock_align.return_value = _good_tm_result(0.6)

        scorer = TMScoreFS(pdb1, PDB1_NAME, pdb2, PDB2_NAME)
        assert scorer.tmscores_fs is not None
        assert scorer.tmscores_fs.shape[0] == 2


class TestTMScoreFSMulti:
    def test_get_coords_correct_range(self, tmp_path):
        pdb = tmp_path / "test.pdb"
        _write_minimal_pdb(pdb, n_residues=10, n_chains=1)
        scorer = object.__new__(TMScoreFSMulti)
        coords, seq = scorer.get_coords(pdb, "1-5")
        assert len(coords) == 5
        assert len(seq) == 5

    def test_get_tmscore_uses_single_files(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        pred_dir = _make_converted_multimer_dir(tmp_path, "pred", n_models=3)
        matched = glob.glob(str(pred_dir) + "/single_0_unrelaxed*pdb")
        assert len(matched) == 3
        assert all("rmTER" not in f for f in matched)

    def test_get_tmscore_does_not_use_rmter_files(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        pred_dir = _make_converted_multimer_dir(tmp_path, "pred", n_models=3)
        matched = glob.glob(str(pred_dir) + "/single_0_unrelaxed*pdb")
        # rmTER files exist but must not be in single_ glob
        rmter = glob.glob(str(pred_dir) + "/rmTER*pdb")
        assert rmter  # they exist
        assert not any("rmTER" in f for f in matched)

    def test_get_tmscore_returns_zeros_on_empty_dir(self, tmp_path):
        empty = tmp_path / "empty"
        empty.mkdir()
        scorer = object.__new__(TMScoreFSMulti)
        coords = RNG.uniform(size=(5, 3))
        scores = scorer.get_tmscore(coords, "AGVLI", empty, "1-5")
        assert scores == [0.0, 0.0, 0.0, 0.0, 0.0]

    @patch("cf_random.analysis.cal_tmscore_fs_multimer.tm_align")
    def test_get_tmscore_forward_only(self, mock_align, tmp_path, monkeypatch):
        """TMScoreFSMulti does not do reverse alignment unlike TMScoreFS."""
        monkeypatch.chdir(tmp_path)
        pred_dir = _make_converted_multimer_dir(tmp_path, "pred", n_models=2, n_residues=5)
        mock_align.return_value = _good_tm_result(0.7)
        scorer = object.__new__(TMScoreFSMulti)
        coords = RNG.uniform(size=(5, 3))
        scores = scorer.get_tmscore(coords, "AGVLI", pred_dir, "1-5")
        # One call per model, no reverse
        assert mock_align.call_count == 2

    @patch("cf_random.analysis.cal_tmscore_fs_multimer.tm_align")
    def test_run_for_models_interleaves_folds(self, mock_align, tmp_path, monkeypatch):
        """run_for_models appends pdb1 row then pdb2 row per directory → (2, n_models)."""
        monkeypatch.chdir(tmp_path)
        _make_range_file(tmp_path / "range_fs_pairs_all.txt")
        pred_dir = _make_converted_multimer_dir(tmp_path, "pred", n_models=5, n_residues=10)
        pdb1 = tmp_path / f"{PDB1_NAME}.pdb"
        pdb2 = tmp_path / f"{PDB2_NAME}.pdb"
        _write_minimal_pdb(pdb1, n_residues=10)
        _write_minimal_pdb(pdb2, n_residues=10)
        mock_align.return_value = _good_tm_result(0.6)

        scorer = object.__new__(TMScoreFSMulti)
        scorer.run_for_models(pdb1, pdb2, str(pred_dir), "1-5", "1-5", "1-5")
        assert scorer.tmscores_fs.shape == (2, 5)

    @patch("cf_random.analysis.cal_tmscore_fs_multimer.tm_align")
    def test_full_init_completes(self, mock_align, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        _make_range_file(tmp_path / "range_fs_pairs_all.txt")
        pred_dir = _make_converted_multimer_dir(tmp_path, "pred", n_residues=10)
        pdb1 = tmp_path / f"{PDB1_NAME}.pdb"
        pdb2 = tmp_path / f"{PDB2_NAME}.pdb"
        _write_minimal_pdb(pdb1, n_residues=10)
        _write_minimal_pdb(pdb2, n_residues=10)
        mock_align.return_value = _good_tm_result(0.6)

        scorer = TMScoreFSMulti(str(pred_dir), pdb1, PDB1_NAME, pdb2, PDB2_NAME)
        assert scorer.tmscores_fs.shape == (2, NUM_PREDICTIONS)


class TestTMScoreGetPredictedFiles:
    def _make_scorer(self, pred_dir, model_type):
        scorer = object.__new__(BaseTMScore)
        scorer.pred_dir = str(pred_dir)
        scorer.pdb1_name = PDB1_NAME
        scorer.pdb2_name = PDB2_NAME
        scorer.model_type = model_type
        return scorer

    def test_monomer_returns_unrelaxed_files(self, tmp_path):
        pred_dir = _make_monomer_prediction_dir(tmp_path, "pred", n_models=5)
        scorer = self._make_scorer(pred_dir, MONOMER_MODEL_TYPE)
        files = scorer._resolve_models()
        assert len(files) == 5
        assert all("unrelaxed" in f for f in files)
        assert all("rmTER" not in f for f in files)

    def test_multimer_returns_rmter_files(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        pred_dir = _make_converted_multimer_dir(tmp_path, "pred")
        scorer = self._make_scorer(pred_dir, MULTIMER_MODEL_TYPE)
        files = scorer._resolve_models()
        assert len(files) == NUM_PREDICTIONS
        assert all("rmTER" in f for f in files)

    def test_multimer_triggers_conversion_when_missing(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        pred_dir = _make_multimer_prediction_dir(tmp_path, "pred")
        pdb2 = tmp_path / f"{PDB2_NAME}.pdb"
        _write_minimal_pdb(pdb2)
        scorer = self._make_scorer(pred_dir, MULTIMER_MODEL_TYPE)
        files = scorer._resolve_models()
        assert len(files) == NUM_PREDICTIONS
        assert list(pred_dir.glob("rmTER_0_unrelaxed*pdb"))

    def test_multimer_does_not_reconvert_if_rmter_exists(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        pred_dir = _make_converted_multimer_dir(tmp_path, "pred", n_models=3)
        rmter = list(pred_dir.glob("rmTER_0_unrelaxed*pdb"))
        rmter[0].write_text("SENTINEL")
        scorer = self._make_scorer(pred_dir, MULTIMER_MODEL_TYPE)
        scorer._resolve_models()
        assert rmter[0].read_text() == "SENTINEL"

    def test_monomer_does_not_create_rmter_files(self, tmp_path):
        pred_dir = _make_monomer_prediction_dir(tmp_path, "pred")
        scorer = self._make_scorer(pred_dir, MONOMER_MODEL_TYPE)
        scorer._resolve_models()
        assert not list(pred_dir.glob("rmTER*"))

    def test_resolve_models_strips_extension_and_prepends_pwd(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        pred_dir = _make_monomer_prediction_dir(tmp_path, "pred", n_models=2)
        scorer = self._make_scorer(pred_dir, MONOMER_MODEL_TYPE)
        files = scorer._resolve_models()
        assert all(not f.endswith(".pdb") for f in files)
        assert all(f.startswith(str(tmp_path)) for f in files)


class TestEvaluateMonomerIntegration:
    @pytest.fixture()
    def monomer_pipeline_dir(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        _make_range_file(tmp_path / "range_fs_pairs_all.txt", fs_range="1-5")
        pdb1 = tmp_path / f"{PDB1_NAME}.pdb"
        pdb2 = tmp_path / f"{PDB2_NAME}.pdb"
        _write_minimal_pdb(pdb1, n_residues=10)
        _write_minimal_pdb(pdb2, n_residues=10)

        pred_root = tmp_path / "predictions_all" / PDB1_NAME
        pred_root.mkdir(parents=True)

        # Full MSA
        _make_monomer_prediction_dir(
            pred_root,
            f"{PDB1_NAME}_predicted_models_full_rand_0",
            n_residues=10,
        )
        # Shallow random MSA dirs
        max_msa, ext_msa = 1, 2
        for mult in MSA_MULTIPLIERS:
            max_msa *= mult
            ext_msa *= mult
            _make_monomer_prediction_dir(
                pred_root,
                f"{PDB1_NAME}_predicted_models_rand_0_max_{max_msa}_ext_{ext_msa}",
                n_residues=10,
            )
        return tmp_path

    @patch("cf_random.analysis.cal_tmscore_fs_flmsa.tm_align")
    @patch("cf_random.analysis.base.tm_align")
    def test_completes(self, mock_whole, mock_fs, monomer_pipeline_dir):
        mock_whole.return_value = _good_tm_result(0.6)
        mock_fs.return_value = _good_tm_result(0.6)
        scorer = TMScoreCalAllVarFS(
            pdb1=f"{PDB1_NAME}.pdb",
            pdb1_name=PDB1_NAME,
            pdb2=f"{PDB2_NAME}.pdb",
            pdb2_name=PDB2_NAME,
            num_msa=0,
            option="FS",
            model_type=MONOMER_MODEL_TYPE,
        )
        assert len(scorer.size_selection) == 1

    @patch("cf_random.analysis.cal_tmscore_fs_flmsa.tm_align")
    @patch("cf_random.analysis.base.tm_align")
    def test_csv_files_written(self, mock_whole, mock_fs, monomer_pipeline_dir):
        mock_whole.return_value = _good_tm_result(0.6)
        mock_fs.return_value = _good_tm_result(0.6)
        TMScoreCalAllVarFS(
            pdb1=f"{PDB1_NAME}.pdb",
            pdb1_name=PDB1_NAME,
            pdb2=f"{PDB2_NAME}.pdb",
            pdb2_name=PDB2_NAME,
            num_msa=0,
            option="FS",
            model_type=MONOMER_MODEL_TYPE,
        )
        base = monomer_pipeline_dir
        assert (base / f"TMScore_full-MSA_{PDB1_NAME}.csv").exists()
        assert (base / f"TMScore_fs_full-MSA_{PDB1_NAME}.csv").exists()
        assert (base / f"TMScore_random-MSA_{PDB1_NAME}.csv").exists()
        assert (base / f"TMScore_fs_random-MSA_{PDB1_NAME}.csv").exists()

    @patch("cf_random.analysis.cal_tmscore_fs_flmsa.tm_align")
    @patch("cf_random.analysis.base.tm_align")
    def test_bad_predictions_raise(self, mock_whole, mock_fs, monomer_pipeline_dir):
        mock_whole.return_value = _good_tm_result(0.1)
        mock_fs.return_value = _good_tm_result(0.1)
        with pytest.raises(RuntimeError):
            TMScoreCalAllVarFS(
                pdb1=f"{PDB1_NAME}.pdb",
                pdb1_name=PDB1_NAME,
                pdb2=f"{PDB2_NAME}.pdb",
                pdb2_name=PDB2_NAME,
                num_msa=0,
                option="FS",
                model_type=MONOMER_MODEL_TYPE,
            )

    @patch("cf_random.analysis.cal_tmscore_fs_flmsa.tm_align")
    @patch("cf_random.analysis.base.tm_align")
    def test_ref_alt_determination(self, mock_whole, mock_fs, monomer_pipeline_dir):
        """Row 0 scores higher → pdb1 is reference, pdb2 is alternative."""
        high = _good_tm_result(0.8)
        low = _good_tm_result(0.2)

        def whole_side_effect(*args, **kwargs):
            # tm_align called alternately for pdb1 then pdb2 rows
            whole_side_effect.count = getattr(whole_side_effect, "count", 0) + 1
            return high if whole_side_effect.count % 2 == 1 else low

        mock_whole.side_effect = whole_side_effect
        mock_fs.return_value = _good_tm_result(0.6)

        scorer = TMScoreCalAllVarFS(
            pdb1=f"{PDB1_NAME}.pdb",
            pdb1_name=PDB1_NAME,
            pdb2=f"{PDB2_NAME}.pdb",
            pdb2_name=PDB2_NAME,
            num_msa=0,
            option="FS",
            model_type=MONOMER_MODEL_TYPE,
        )
        assert len(scorer.size_selection) == 1


class TestEvaluateMultimerIntegration:
    @pytest.fixture()
    def multimer_pipeline_dir(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        _make_range_file(tmp_path / "range_fs_pairs_all.txt", fs_range="1-5")
        pdb1 = tmp_path / f"{PDB1_NAME}.pdb"
        pdb2 = tmp_path / f"{PDB2_NAME}.pdb"
        _write_minimal_pdb(pdb1, n_residues=10)
        _write_minimal_pdb(pdb2, n_residues=10)

        pred_root = tmp_path / "predictions_all" / PDB1_NAME
        pred_root.mkdir(parents=True)

        # Full MSA
        full_name = f"{PDB1_NAME}_predicted_models_full_rand_0"
        _make_multimer_prediction_dir(pred_root, full_name, n_residues=10)
        ConvertM2S(str(pred_root / full_name), PDB1_NAME, PDB2_NAME)

        # Shallow random MSA dirs
        max_msa, ext_msa = 1, 2
        for mult in MSA_MULTIPLIERS:
            max_msa *= mult
            ext_msa *= mult
            rand_name = f"{PDB1_NAME}_predicted_models_rand_0_max_{max_msa}_ext_{ext_msa}"
            _make_multimer_prediction_dir(pred_root, rand_name, n_residues=10)
            ConvertM2S(str(pred_root / rand_name), PDB1_NAME, PDB2_NAME)

        return tmp_path

    @patch("cf_random.analysis.cal_tmscore_fs_multimer.tm_align")
    @patch("cf_random.analysis.base.tm_align")
    def test_completes(self, mock_whole, mock_fs, multimer_pipeline_dir):
        mock_whole.return_value = _good_tm_result(0.6)
        mock_fs.return_value = _good_tm_result(0.6)
        scorer = TMScoreCalAllVarFS(
            pdb1=f"{PDB1_NAME}.pdb",
            pdb1_name=PDB1_NAME,
            pdb2=f"{PDB2_NAME}.pdb",
            pdb2_name=PDB2_NAME,
            num_msa=0,
            option="FS",
            model_type=MULTIMER_MODEL_TYPE,
        )
        assert len(scorer.size_selection) == 1

    @patch("cf_random.analysis.cal_tmscore_fs_multimer.tm_align")
    @patch("cf_random.analysis.base.tm_align")
    def test_csv_files_written(self, mock_whole, mock_fs, multimer_pipeline_dir):
        mock_whole.return_value = _good_tm_result(0.6)
        mock_fs.return_value = _good_tm_result(0.6)
        TMScoreCalAllVarFS(
            pdb1=f"{PDB1_NAME}.pdb",
            pdb1_name=PDB1_NAME,
            pdb2=f"{PDB2_NAME}.pdb",
            pdb2_name=PDB2_NAME,
            num_msa=0,
            option="FS",
            model_type=MULTIMER_MODEL_TYPE,
        )
        base = multimer_pipeline_dir
        assert (base / f"TMScore_full-MSA_{PDB1_NAME}.csv").exists()
        assert (base / f"TMScore_fs_full-MSA_{PDB1_NAME}.csv").exists()
        assert (base / f"TMScore_random-MSA_{PDB1_NAME}.csv").exists()
        assert (base / f"TMScore_fs_random-MSA_{PDB1_NAME}.csv").exists()

    @patch("cf_random.analysis.cal_tmscore_fs_multimer.tm_align")
    @patch("cf_random.analysis.base.tm_align")
    def test_size_selection_is_int(self, mock_whole, mock_fs, multimer_pipeline_dir):
        mock_whole.return_value = _good_tm_result(0.6)
        mock_fs.return_value = _good_tm_result(0.6)
        scorer = TMScoreCalAllVarFS(
            pdb1=f"{PDB1_NAME}.pdb",
            pdb1_name=PDB1_NAME,
            pdb2=f"{PDB2_NAME}.pdb",
            pdb2_name=PDB2_NAME,
            num_msa=0,
            option="FS",
            model_type=MULTIMER_MODEL_TYPE,
        )
        assert isinstance(scorer.size_selection[0], int)

    @patch("cf_random.analysis.cal_tmscore_fs_multimer.tm_align")
    @patch("cf_random.analysis.base.tm_align")
    def test_bad_predictions_raise(self, mock_whole, mock_fs, multimer_pipeline_dir):
        mock_whole.return_value = _good_tm_result(0.1)
        mock_fs.return_value = _good_tm_result(0.1)
        with pytest.raises(RuntimeError):
            TMScoreCalAllVarFS(
                pdb1=f"{PDB1_NAME}.pdb",
                pdb1_name=PDB1_NAME,
                pdb2=f"{PDB2_NAME}.pdb",
                pdb2_name=PDB2_NAME,
                num_msa=0,
                option="FS",
                model_type=MULTIMER_MODEL_TYPE,
            )

    @patch("cf_random.analysis.cal_tmscore_fs_multimer.tm_align")
    @patch("cf_random.analysis.base.tm_align")
    def test_multimer_uses_rmter_for_whole_structure(
        self, mock_whole, mock_fs, multimer_pipeline_dir
    ):
        seen_files = []
        original_resolve = BaseTMScore._resolve_models

        def capturing_resolve(self_inner):
            files = original_resolve(self_inner)
            seen_files.extend(files)
            return files

        mock_whole.return_value = _good_tm_result(0.6)
        mock_fs.return_value = _good_tm_result(0.6)

        with patch.object(BaseTMScore, "_resolve_models", capturing_resolve):
            TMScoreCalAllVarFS(
                pdb1=f"{PDB1_NAME}.pdb",
                pdb1_name=PDB1_NAME,
                pdb2=f"{PDB2_NAME}.pdb",
                pdb2_name=PDB2_NAME,
                num_msa=0,
                option="FS",
                model_type=MULTIMER_MODEL_TYPE,
            )

        assert seen_files
        assert all("rmTER" in f for f in seen_files)

    @patch("cf_random.analysis.cal_tmscore_fs_multimer.tm_align")
    @patch("cf_random.analysis.base.tm_align")
    def test_multimer_uses_single_for_fs_region(self, mock_whole, mock_fs, multimer_pipeline_dir):
        """FS scorer must receive single_ files, not rmTER or originals."""
        seen_files = []

        original_get = TMScoreFSMulti.get_tmscore

        def capturing_get(self_inner, coords1, seq1, predfilepath, res_range):
            matched = glob.glob(str(predfilepath) + "/single_0_unrelaxed*pdb")
            seen_files.extend(matched)
            return original_get(self_inner, coords1, seq1, predfilepath, res_range)

        mock_whole.return_value = _good_tm_result(0.6)
        mock_fs.return_value = _good_tm_result(0.6)

        with patch.object(TMScoreFSMulti, "get_tmscore", capturing_get):
            TMScoreCalAllVarFS(
                pdb1=f"{PDB1_NAME}.pdb",
                pdb1_name=PDB1_NAME,
                pdb2=f"{PDB2_NAME}.pdb",
                pdb2_name=PDB2_NAME,
                num_msa=0,
                option="FS",
                model_type=MULTIMER_MODEL_TYPE,
            )

        assert seen_files
        assert all("single_" in f for f in seen_files)
        assert all("rmTER" not in f for f in seen_files)
