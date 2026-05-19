#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Tests for alternative conformation TM-score pipelines.

Unit tests cover TMScore, TMScoreCalAllVar._evaluate_monomer /
_evaluate_multimer, and select_size logic.

Run with:
    pytest tests/test_alternative_conformation.py -v
"""

from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np
import pytest

from cf_random.utils.convert_multi_single import ConvertM2S
from cf_random.analysis.tmscore_all_var import TMScore, TMScoreCalAllVar, MSA_MULTIPLIERS
from cf_random.analysis.base import BaseTMScore

PDB1_NAME = "5olw_A"
PDB2_NAME = "5olx_A"
MONOMER_MODEL_TYPE = "alphafold2_ptm"
MULTIMER_MODEL_TYPE = "alphafold2_multimer_v3"
NUM_SEEDS = 5
NUM_MODELS = 5
NUM_PREDICTIONS = NUM_SEEDS * NUM_MODELS  # 25

RNG = np.random.default_rng(42)

RESNAMES = ["ALA", "GLY", "VAL", "LEU", "ILE", "PRO", "PHE", "TRP", "MET", "SER"]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _write_minimal_pdb(path: Path, n_residues: int = 10, n_chains: int = 1) -> None:
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
    pred_dir = _make_multimer_prediction_dir(base, name, n_models=n_models, n_residues=n_residues)
    pdb2_file = base / f"{PDB2_NAME}.pdb"
    if not pdb2_file.exists():
        _write_minimal_pdb(pdb2_file, n_residues=n_residues, n_chains=1)
    ConvertM2S(str(pred_dir), PDB1_NAME, PDB2_NAME)
    return pred_dir


def _good_tm_result(value: float = 0.6) -> MagicMock:
    result = MagicMock()
    result.tm_norm_chain1 = value
    return result


def _make_pipeline_dirs(tmp_path, multimer: bool = False) -> Path:
    """Create full prediction directory structure for integration tests."""
    pdb1 = tmp_path / f"{PDB1_NAME}.pdb"
    pdb2 = tmp_path / f"{PDB2_NAME}.pdb"
    _write_minimal_pdb(pdb1, n_residues=10)
    _write_minimal_pdb(pdb2, n_residues=10)

    pred_root = tmp_path / "predictions_all" / PDB1_NAME
    pred_root.mkdir(parents=True)

    if multimer:
        full_name = f"{PDB1_NAME}_predicted_models_full_rand_0"
        _make_multimer_prediction_dir(pred_root, full_name, n_residues=10)
        ConvertM2S(str(pred_root / full_name), PDB1_NAME, PDB2_NAME)
        max_msa, ext_msa = 1, 2
        for mult in MSA_MULTIPLIERS:
            max_msa *= mult
            ext_msa *= mult
            rand_name = f"{PDB1_NAME}_predicted_models_rand_0_max_{max_msa}_ext_{ext_msa}"
            _make_multimer_prediction_dir(pred_root, rand_name, n_residues=10)
            ConvertM2S(str(pred_root / rand_name), PDB1_NAME, PDB2_NAME)
    else:
        _make_monomer_prediction_dir(
            pred_root, f"{PDB1_NAME}_predicted_models_full_rand_0", n_residues=10
        )
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


# ---------------------------------------------------------------------------
# TMScore.select_size
# ---------------------------------------------------------------------------


class TestSelectSize:
    def _make_scorer(self) -> TMScore:
        scorer = object.__new__(TMScore)
        scorer.pdb1_name = PDB1_NAME
        scorer.pdb2_name = PDB2_NAME
        return scorer

    def test_selects_pdb2_as_alternative(self):
        scorer = self._make_scorer()
        num_seeds = 5
        # Row 1 (pdb2) has high scores at MSA depth 0 → selection should be 0
        data = np.zeros((14, num_seeds * 5))
        data[1, :] = 0.8
        scorer.select_size(data.flatten(), PDB1_NAME, PDB2_NAME, PDB2_NAME, num_seeds)
        assert scorer.selection == 0

    def test_selects_pdb1_as_alternative(self):
        scorer = self._make_scorer()
        num_seeds = 5
        # Row 0 (pdb1) has high scores at MSA depth 0 → selection should be 0
        data = np.zeros((14, num_seeds * 5))
        data[0, :] = 0.8
        scorer.select_size(data.flatten(), PDB1_NAME, PDB2_NAME, PDB1_NAME, num_seeds)
        assert scorer.selection == 0

    def test_raises_when_no_scores_above_threshold(self):
        scorer = self._make_scorer()
        num_seeds = 5
        data = np.full((14, num_seeds * 5), 0.3)
        with pytest.raises(RuntimeError):
            scorer.select_size(data.flatten(), PDB1_NAME, PDB2_NAME, PDB2_NAME, num_seeds)

    def test_selection_is_int(self):
        scorer = self._make_scorer()
        num_seeds = 5
        data = np.zeros((14, num_seeds * 5))
        data[1, :] = 0.8
        scorer.select_size(data.flatten(), PDB1_NAME, PDB2_NAME, PDB2_NAME, num_seeds)
        assert isinstance(scorer.selection, int)

    def test_picks_highest_sum_depth(self):
        """Row with highest sum across seeds should be selected."""
        scorer = self._make_scorer()
        num_seeds = 5
        data = np.zeros((14, num_seeds * 5))
        # MSA depth 2 (row index 5 for pdb2) has highest scores
        data[5, :] = 0.9
        data[1, :] = 0.6
        scorer.select_size(data.flatten(), PDB1_NAME, PDB2_NAME, PDB2_NAME, num_seeds)
        assert scorer.selection == 2


# ---------------------------------------------------------------------------
# TMScoreCalAllVar._determine_alternative / _extract_alternative_rows
# ---------------------------------------------------------------------------


class TestHelperMethods:
    def _make_cal(self) -> TMScoreCalAllVar:
        cal = object.__new__(TMScoreCalAllVar)
        cal.pdb1_name = PDB1_NAME
        cal.pdb2_name = PDB2_NAME
        return cal

    def test_determine_alternative_row0_higher(self):
        cal = self._make_cal()
        scores = np.array([0.8, 0.3])
        assert cal._determine_alternative(scores) == PDB2_NAME

    def test_determine_alternative_row1_higher(self):
        cal = self._make_cal()
        scores = np.array([0.3, 0.8])
        assert cal._determine_alternative(scores) == PDB1_NAME

    def test_determine_alternative_equal_scores(self):
        """Equal scores → pdb2 is alternative (row0 >= row1 condition)."""
        cal = self._make_cal()
        scores = np.array([0.5, 0.5])
        assert cal._determine_alternative(scores) == PDB2_NAME

    def test_extract_alternative_rows_pdb2(self):
        matrix = np.arange(14 * 5).reshape(14, 5).astype(float)
        result = TMScoreCalAllVar._extract_alternative_rows(matrix, PDB2_NAME, PDB1_NAME, PDB2_NAME)
        expected = matrix[1::2, :]
        np.testing.assert_array_equal(result, expected)

    def test_extract_alternative_rows_pdb1(self):
        matrix = np.arange(14 * 5).reshape(14, 5).astype(float)
        result = TMScoreCalAllVar._extract_alternative_rows(matrix, PDB1_NAME, PDB1_NAME, PDB2_NAME)
        expected = matrix[0::2, :]
        np.testing.assert_array_equal(result, expected)


# ---------------------------------------------------------------------------
# Monomer integration
# ---------------------------------------------------------------------------


class TestEvaluateMonomerIntegration:
    @pytest.fixture()
    def monomer_dir(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        return _make_pipeline_dirs(tmp_path, multimer=False)

    @patch("cf_random.analysis.base.tm_align")
    def test_completes(self, mock_align, monomer_dir):
        mock_align.return_value = _good_tm_result(0.6)
        scorer = TMScoreCalAllVar(
            pdb1=f"{PDB1_NAME}.pdb",
            pdb1_name=PDB1_NAME,
            pdb2=f"{PDB2_NAME}.pdb",
            pdb2_name=PDB2_NAME,
            num_msa=0,
            option="AC",
            model_type=MONOMER_MODEL_TYPE,
        )
        assert len(scorer.size_selection) == 1

    @patch("cf_random.analysis.base.tm_align")
    def test_csv_files_written(self, mock_align, monomer_dir):
        mock_align.return_value = _good_tm_result(0.6)
        TMScoreCalAllVar(
            pdb1=f"{PDB1_NAME}.pdb",
            pdb1_name=PDB1_NAME,
            pdb2=f"{PDB2_NAME}.pdb",
            pdb2_name=PDB2_NAME,
            num_msa=0,
            option="AC",
            model_type=MONOMER_MODEL_TYPE,
        )
        assert (monomer_dir / f"TMScore_full-MSA_{PDB1_NAME}.csv").exists()
        assert (monomer_dir / f"TMScore_random-MSA_{PDB1_NAME}.csv").exists()

    @patch("cf_random.analysis.base.tm_align")
    def test_no_fs_csv_written(self, mock_align, monomer_dir):
        """AC mode must not write FS-region CSV files."""
        mock_align.return_value = _good_tm_result(0.6)
        TMScoreCalAllVar(
            pdb1=f"{PDB1_NAME}.pdb",
            pdb1_name=PDB1_NAME,
            pdb2=f"{PDB2_NAME}.pdb",
            pdb2_name=PDB2_NAME,
            num_msa=0,
            option="AC",
            model_type=MONOMER_MODEL_TYPE,
        )
        assert not (monomer_dir / f"TMScore_fs_full-MSA_{PDB1_NAME}.csv").exists()
        assert not (monomer_dir / f"TMScore_fs_random-MSA_{PDB1_NAME}.csv").exists()

    @patch("cf_random.analysis.base.tm_align")
    def test_bad_full_msa_raises(self, mock_align, monomer_dir):
        mock_align.return_value = _good_tm_result(0.1)
        with pytest.raises(RuntimeError):
            TMScoreCalAllVar(
                pdb1=f"{PDB1_NAME}.pdb",
                pdb1_name=PDB1_NAME,
                pdb2=f"{PDB2_NAME}.pdb",
                pdb2_name=PDB2_NAME,
                num_msa=0,
                option="AC",
                model_type=MONOMER_MODEL_TYPE,
            )

    @patch("cf_random.analysis.base.tm_align")
    def test_failed_predictions_moved(self, mock_align, monomer_dir):
        mock_align.return_value = _good_tm_result(0.1)
        with pytest.raises(RuntimeError):
            TMScoreCalAllVar(
                pdb1=f"{PDB1_NAME}.pdb",
                pdb1_name=PDB1_NAME,
                pdb2=f"{PDB2_NAME}.pdb",
                pdb2_name=PDB2_NAME,
                num_msa=0,
                option="AC",
                model_type=MONOMER_MODEL_TYPE,
            )
        failed_root = monomer_dir / "failed_predictions" / PDB1_NAME
        assert failed_root.exists()
        assert any(failed_root.iterdir())

    @patch("cf_random.analysis.base.tm_align")
    def test_size_selection_is_int(self, mock_align, monomer_dir):
        mock_align.return_value = _good_tm_result(0.6)
        scorer = TMScoreCalAllVar(
            pdb1=f"{PDB1_NAME}.pdb",
            pdb1_name=PDB1_NAME,
            pdb2=f"{PDB2_NAME}.pdb",
            pdb2_name=PDB2_NAME,
            num_msa=0,
            option="AC",
            model_type=MONOMER_MODEL_TYPE,
        )
        assert isinstance(scorer.size_selection[0], int)

    @patch("cf_random.analysis.base.tm_align")
    def test_ref_alt_determination(self, mock_align, monomer_dir):
        """Row 0 averaging higher → pdb2 is alternative."""
        high = _good_tm_result(0.8)
        low = _good_tm_result(0.2)

        def side_effect(*args, **kwargs):
            side_effect.count = getattr(side_effect, "count", 0) + 1
            return high if side_effect.count % 2 == 1 else low

        mock_align.side_effect = side_effect
        scorer = TMScoreCalAllVar(
            pdb1=f"{PDB1_NAME}.pdb",
            pdb1_name=PDB1_NAME,
            pdb2=f"{PDB2_NAME}.pdb",
            pdb2_name=PDB2_NAME,
            num_msa=0,
            option="AC",
            model_type=MONOMER_MODEL_TYPE,
        )
        assert len(scorer.size_selection) == 1


# ---------------------------------------------------------------------------
# Multimer integration
# ---------------------------------------------------------------------------


class TestEvaluateMultimerIntegration:
    @pytest.fixture()
    def multimer_dir(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        return _make_pipeline_dirs(tmp_path, multimer=True)

    @patch("cf_random.analysis.base.tm_align")
    def test_completes(self, mock_align, multimer_dir):
        mock_align.return_value = _good_tm_result(0.6)
        scorer = TMScoreCalAllVar(
            pdb1=f"{PDB1_NAME}.pdb",
            pdb1_name=PDB1_NAME,
            pdb2=f"{PDB2_NAME}.pdb",
            pdb2_name=PDB2_NAME,
            num_msa=0,
            option="AC",
            model_type=MULTIMER_MODEL_TYPE,
        )
        assert len(scorer.size_selection) == 1

    @patch("cf_random.analysis.base.tm_align")
    def test_csv_files_written(self, mock_align, multimer_dir):
        mock_align.return_value = _good_tm_result(0.6)
        TMScoreCalAllVar(
            pdb1=f"{PDB1_NAME}.pdb",
            pdb1_name=PDB1_NAME,
            pdb2=f"{PDB2_NAME}.pdb",
            pdb2_name=PDB2_NAME,
            num_msa=0,
            option="AC",
            model_type=MULTIMER_MODEL_TYPE,
        )
        assert (multimer_dir / f"TMScore_full-MSA_{PDB1_NAME}.csv").exists()
        assert (multimer_dir / f"TMScore_random-MSA_{PDB1_NAME}.csv").exists()

    @patch("cf_random.analysis.base.tm_align")
    def test_bad_full_msa_raises(self, mock_align, multimer_dir):
        mock_align.return_value = _good_tm_result(0.1)
        with pytest.raises(RuntimeError):
            TMScoreCalAllVar(
                pdb1=f"{PDB1_NAME}.pdb",
                pdb1_name=PDB1_NAME,
                pdb2=f"{PDB2_NAME}.pdb",
                pdb2_name=PDB2_NAME,
                num_msa=0,
                option="AC",
                model_type=MULTIMER_MODEL_TYPE,
            )

    @patch("cf_random.analysis.base.tm_align")
    def test_size_selection_is_int(self, mock_align, multimer_dir):
        mock_align.return_value = _good_tm_result(0.6)
        scorer = TMScoreCalAllVar(
            pdb1=f"{PDB1_NAME}.pdb",
            pdb1_name=PDB1_NAME,
            pdb2=f"{PDB2_NAME}.pdb",
            pdb2_name=PDB2_NAME,
            num_msa=0,
            option="AC",
            model_type=MULTIMER_MODEL_TYPE,
        )
        assert isinstance(scorer.size_selection[0], int)

    @patch("cf_random.analysis.base.tm_align")
    def test_multimer_uses_rmter_files(self, mock_align, multimer_dir):
        """BaseTMScore must resolve rmTER files for multimer whole-structure scoring."""
        seen_files = []
        original_resolve = BaseTMScore._resolve_models

        def capturing_resolve(self_inner):
            files = original_resolve(self_inner)
            seen_files.extend(files)
            return files

        mock_align.return_value = _good_tm_result(0.6)

        with patch.object(BaseTMScore, "_resolve_models", capturing_resolve):
            TMScoreCalAllVar(
                pdb1=f"{PDB1_NAME}.pdb",
                pdb1_name=PDB1_NAME,
                pdb2=f"{PDB2_NAME}.pdb",
                pdb2_name=PDB2_NAME,
                num_msa=0,
                option="AC",
                model_type=MULTIMER_MODEL_TYPE,
            )

        assert seen_files
        assert all("rmTER" in f for f in seen_files)
