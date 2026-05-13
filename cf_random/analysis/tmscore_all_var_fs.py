#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""TM-score analysis for fold-switching prediction workflows."""

import glob
import shutil
from pathlib import (
    Path,
)
from typing import (
    List,
    Optional,
)

import numpy as np

from ..prediction.pred_cal_tmscore_multimer import (
    PredictionAllMultimerFS,
)
from ..analysis.cal_tmscore_fs_multimer import (
    TMScoreFSMulti,
)
from .base import (
    MULTIMER_MODEL_TYPE,
    BaseTMScore,
    logger,
)
from .cal_tmscore_fs_flmsa import (
    TMScoreFS,
)

PREDICTIONS_ROOT = Path("predictions_all")
FAILED_ROOT = Path("failed_prediction")

# Must match the multiplier sequence used in prediction_all_var.py
MSA_MULTIPLIERS = (1, 2, 2, 2, 2, 2, 2)


class TMScore(BaseTMScore):
    """Compute whole-structure TM-scores for a prediction set.

    Inherits forward/reverse scoring, glob resolution, and multimer
    conversion from BaseTMScore. Adds fold-switching select_size logic.
    """

    def select_size(
        self,
        tmscores_random: np.ndarray,
        tmscores_fs_random: np.ndarray,
        pdb1_name: str,
        pdb2_name: str,
        alt_name: str,
        num_seeds: int,
    ) -> None:
        """Select optimal shallow MSA depth using whole-structure and FS scores.

            1. Extract alternative rows from both whole and FS matrices.
            2. Sum whole-structure group scores and pick argmax location.
            3. Primary check: whole >= 0.5 AND fs >= 0.5 at that location.
            4. Fallback: if fs < 0.5 at primary location, scan all pairs
               for any combination where whole >= 0.4 and fs >= 0.5.
            5. Raise RuntimeError if no valid selection found.

        Args:
            tmscores_random: Flat array of shallow whole-structure TM-scores.
            tmscores_fs_random: Flat array of shallow FS-region TM-scores.
            pdb1_name: Name of first reference structure.
            pdb2_name: Name of second reference structure.
            alt_name: Name of the alternative conformation structure.
            num_seeds: Number of prediction seeds.

        Raises:
            RuntimeError: If no viable shallow MSA selection can be found.
        """
        tmscores_random_reshape = tmscores_random.reshape(14, num_seeds * 5)
        tmscores_fs_random_reshape = tmscores_fs_random.reshape(14, num_seeds * 5)
        tmscores_random_locat = np.zeros((7, num_seeds * 5))
        tmscores_fs_random_locat = np.zeros((7, num_seeds * 5))

        # Extract every-other row for the alternative structure (both matrices)
        if alt_name == pdb2_name:
            tmp_cnt = 0
            for i in range(1, 14, 2):
                tmscores_random_locat[tmp_cnt, :] = tmscores_random_reshape[i, :]
                tmscores_fs_random_locat[tmp_cnt, :] = tmscores_fs_random_reshape[i, :]
                tmp_cnt += 1
        else:
            tmp_cnt = 0
            for i in range(0, 13, 2):
                tmscores_random_locat[tmp_cnt, :] = tmscores_random_reshape[i, :]
                tmscores_fs_random_locat[tmp_cnt, :] = tmscores_fs_random_reshape[i, :]
                tmp_cnt += 1

        # Sum whole-structure groups and find best location
        tmscore_data_sum = np.zeros((7, 1))
        for i in range(tmscores_random_locat.shape[0]):
            tmscore_data_sum[i] = np.sum(tmscores_random_locat[i])

        location = int(np.argmax(np.max(tmscore_data_sum, axis=1)))

        # Map group index back to full-matrix row
        if alt_name == pdb2_name:
            location_full = (location * 2) + 1
        else:
            location_full = location * 2

        tmscore_data = tmscores_random_reshape
        tmscore_fs_data = tmscores_fs_random_reshape

        # Primary check: both whole and FS must have >= 0.5 at location
        if alt_name == pdb2_name and (
            np.any(tmscore_data[location_full, :] >= 0.5)
            and np.any(tmscore_fs_data[location_full, :] >= 0.5)
        ):
            self.selection = int((location_full - 1) / 2)
            logger.info("Selected shallow MSA index %s for %s", self.selection, pdb1_name)

        elif alt_name == pdb1_name and (
            np.any(tmscore_data[location_full, :] >= 0.5)
            and np.any(tmscore_fs_data[location_full, :] >= 0.5)
        ):
            self.selection = int(location_full / 2)
            logger.info("Selected shallow MSA index %s for %s", self.selection, pdb1_name)

        # Scan all pairs for whole >= 0.4 AND fs >= 0.5 in any combination
        elif np.any(tmscore_fs_data[location_full, :] < 0.5):
            found = False
            for jj in range(int(tmscore_data.shape[0] / 2)):
                cross1 = np.any(tmscore_data[jj * 2, :] >= 0.4) and np.any(
                    tmscore_fs_data[(jj * 2) + 1, :] >= 0.5
                )
                cross2 = np.any(tmscore_data[(jj * 2) + 1, :] >= 0.4) and np.any(
                    tmscore_fs_data[jj * 2, :] >= 0.5
                )
                same1 = np.any(tmscore_data[jj * 2, :] >= 0.4) and np.any(
                    tmscore_fs_data[jj * 2, :] >= 0.5
                )
                same2 = np.any(tmscore_data[(jj * 2) + 1, :] >= 0.4) and np.any(
                    tmscore_fs_data[(jj * 2) + 1, :] >= 0.5
                )
                if cross1 or cross2 or same1 or same2:
                    self.selection = jj
                    found = True
                    logger.info("Fallback selection: shallow MSA index %s for %s", jj, pdb1_name)
                    break
                elif jj == (int(tmscore_data.shape[0]) - 1) and np.all(tmscore_data[jj, :] < 0.5):
                    raise RuntimeError("Predictions are bad: no viable shallow MSA selection found")
            if not found:
                raise RuntimeError("Predictions are bad: no viable shallow MSA selection found")

        else:
            raise RuntimeError("Predictions are bad: whole-structure predictions are bad")


class TMScoreCalAllVarFS:
    """Evaluate TM-scores for fold-switching prediction results."""

    def __init__(
        self,
        pdb1: str,
        pdb1_name: str,
        pdb2: str,
        pdb2_name: str,
        num_msa: int,
        option: str,
        model_type: str,
        search_dir: Optional[str] = None,
        search_multi_dir: Optional[str] = None,
    ) -> None:
        self.pdb1 = pdb1
        self.pdb2 = pdb2
        self.pdb1_name = pdb1_name
        self.pdb2_name = pdb2_name
        self.num_msa = num_msa
        self.option = option
        self.model_type = model_type
        self.search_dir = search_dir
        self.search_multi_dir = search_multi_dir
        self.size_selection: List[int] = []

        if self.model_type != MULTIMER_MODEL_TYPE:
            self._evaluate_monomer()
        else:
            self._evaluate_multimer()

    def _move_failed_full_outputs(self) -> None:
        """Move failed full-MSA prediction folders to failed_prediction/."""
        FAILED_ROOT.mkdir(parents=True, exist_ok=True)
        pattern = (
            str(PREDICTIONS_ROOT / self.pdb1_name)
            + f"/{self.pdb1_name}_predicted_models_full_rand_*"
        )
        for candidate in glob.glob(pattern):
            destination = FAILED_ROOT / self.pdb1_name / Path(candidate).name
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.move(candidate, str(destination))
            logger.info("Moved failed prediction %s -> %s", candidate, destination)

    @staticmethod
    def _msa_pairs() -> List[tuple]:
        """Generate (max_msa, ext_msa) pairs matching the prediction multiplier sequence."""
        pairs: List[tuple] = []
        max_msa = 1
        ext_msa = 2
        for multiplier in MSA_MULTIPLIERS:
            max_msa *= multiplier
            ext_msa *= multiplier
            pairs.append((max_msa, ext_msa))
        return pairs

    def _evaluate_monomer(self) -> None:
        """Run the full monomer FS TM-score evaluation pipeline."""
        num_seeds = 5 + self.num_msa

        # Full MSA whole-structure TM-scores
        # Passed as glob pattern string so BaseTMScore._resolve_models expands it
        pdb1_basename = self.pdb1_name.split("/")[-1]

        full_pred_dir = (
            str(PREDICTIONS_ROOT / self.pdb1_name)
            + f"/{pdb1_basename}_predicted_models_full_rand_*"
        )
        msa_full_tmscore = TMScore(
            full_pred_dir,
            self.pdb1,
            self.pdb1_name,
            self.pdb2,
            self.pdb2_name,
            self.model_type,
        )
        full_tmscore = np.asarray(msa_full_tmscore.tmscores, dtype=float).reshape(2, num_seeds * 5)

        # Full MSA fold-switching region TM-scores
        msa_fs_tmscore = TMScoreFS(
            self.pdb1,
            self.pdb1_name,
            self.pdb2,
            self.pdb2_name,
        )
        fs_tmscore = np.asarray(msa_fs_tmscore.tmscores_fs, dtype=float).reshape(2, num_seeds * 5)

        # Two-branch quality check using both whole and FS scores
        if np.average(full_tmscore[0, :]) > np.average(full_tmscore[1, :]):
            if np.any(fs_tmscore[0, :] >= 0.5) and np.any(full_tmscore[0, :] >= 0.5):
                ref_name = self.pdb1_name
                alt_name = self.pdb2_name
            elif np.any(fs_tmscore[1, :] >= 0.5) and np.any(full_tmscore[1, :] >= 0.5):
                ref_name = self.pdb2_name
                alt_name = self.pdb1_name
            else:
                self._move_failed_full_outputs()
                raise RuntimeError("Prediction with deep MSA was failed")
        else:
            if np.any(fs_tmscore[1, :] >= 0.5) and np.any(full_tmscore[1, :] >= 0.5):
                ref_name = self.pdb2_name
                alt_name = self.pdb1_name
            elif np.any(fs_tmscore[0, :] >= 0.5) and np.any(full_tmscore[0, :] >= 0.5):
                ref_name = self.pdb1_name
                alt_name = self.pdb2_name
            else:
                self._move_failed_full_outputs()
                raise RuntimeError("Prediction with deep MSA was failed")

        logger.info("Reference structure: %s", ref_name)
        logger.info("Alternative structure: %s", alt_name)

        # Save full MSA scores for both whole structure and FS region
        np.savetxt(f"TMScore_full-MSA_{self.pdb1_name}.csv", full_tmscore, fmt="%2.3f")
        np.savetxt(f"TMScore_fs_full-MSA_{self.pdb1_name}.csv", fs_tmscore, fmt="%2.3f")

        # Shallow random MSA TM-scores (whole + FS)
        tmscores_random: List[float] = []
        tmscores_fs_random: List[float] = []
        last_shallow: Optional[TMScore] = None

        for max_msa, ext_msa in self._msa_pairs():
            pred_dir = (
                str(PREDICTIONS_ROOT / self.pdb1_name)
                + f"/{self.pdb1_name}_predicted_models_rand_*"
                + f"_max_{max_msa}_ext_{ext_msa}"
            )
            logger.debug("Shallow MSA dir pattern: %s", pred_dir)

            shallow = TMScore(
                pred_dir,
                self.pdb1,
                self.pdb1_name,
                self.pdb2,
                self.pdb2_name,
                self.model_type,
            )
            tmscores_random = list(np.append(tmscores_random, shallow.tmscores))
            last_shallow = shallow

            shallow_fs = TMScoreFS(
                self.pdb1,
                self.pdb1_name,
                self.pdb2,
                self.pdb2_name,
                pred_dir_override=pred_dir,
            )
            tmscores_fs_random = list(np.append(tmscores_fs_random, shallow_fs.tmscores_fs))

        random_array = np.asarray(tmscores_random, dtype=float)
        fs_random_array = np.asarray(tmscores_fs_random, dtype=float)

        tmscores_random_reshape = random_array.reshape(14, num_seeds * 5)
        tmscores_fs_random_reshape = fs_random_array.reshape(14, num_seeds * 5)

        # Extract alternative rows for quality check (both matrices)
        if alt_name == self.pdb2_name:
            tmscores_random_alter = tmscores_random_reshape[1::2, :]
            tmscores_fs_random_alter = tmscores_fs_random_reshape[1::2, :]
        else:
            tmscores_random_alter = tmscores_random_reshape[0::2, :]
            tmscores_fs_random_alter = tmscores_fs_random_reshape[0::2, :]

        # Quality gate: both whole and FS must have at least one value > 0.5
        if np.any(tmscores_random_alter > 0.5) and np.any(tmscores_fs_random_alter > 0.5):
            np.savetxt(
                f"TMScore_random-MSA_{self.pdb1_name}.csv",
                tmscores_random_reshape,
                fmt="%2.3f",
            )
            np.savetxt(
                f"TMScore_fs_random-MSA_{self.pdb1_name}.csv",
                tmscores_fs_random_reshape,
                fmt="%2.3f",
            )

            logger.info("Finding optimal size of random MSA...")
            last_shallow.select_size(
                tmscores_random_reshape,
                tmscores_fs_random_reshape,
                self.pdb1_name,
                self.pdb2_name,
                alt_name,
                num_seeds,
            )
            self.size_selection = [last_shallow.selection]
            logger.info("Selected MSA size index: %s", self.size_selection)
        else:
            raise RuntimeError(
                "Full-MSA prediction is not tightly aligned to crystal structure "
                "with additional seeds"
            )

    def _evaluate_multimer(self) -> None:
        """Run the full multimer FS TM-score evaluation pipeline."""
            num_seeds = 5 + self.num_msa
            pdb1_basename = self.pdb1_name.split("/")[-1]

            full_pred_dir = (
                str(PREDICTIONS_ROOT / self.pdb1_name)
                + f"/{pdb1_basename}_predicted_models_full_rand_*"
            )
            msa_full_tmscore = TMScore(
                full_pred_dir,
                self.pdb1,
                self.pdb1_name,
                self.pdb2,
                self.pdb2_name,
                self.model_type,
            )
            full_tmscore = np.asarray(msa_full_tmscore.tmscores, dtype=float).reshape(
                2, num_seeds * 5
            )

            msa_fs_tmscore = TMScoreFSMulti(
                full_pred_dir,
                self.pdb1,
                self.pdb1_name,
                self.pdb2,
                self.pdb2_name,
            )
            fs_tmscore = np.asarray(msa_fs_tmscore.tmscores_fs, dtype=float).reshape(
                2, num_seeds * 5
            )

            if np.average(full_tmscore[0, :]) > np.average(full_tmscore[1, :]):
                if np.any(fs_tmscore[0, :] >= 0.5) and np.any(full_tmscore[0, :] >= 0.5):
                    ref_name = self.pdb1_name
                    alt_name = self.pdb2_name
                elif np.any(fs_tmscore[1, :] >= 0.5) and np.any(full_tmscore[1, :] >= 0.5):
                    ref_name = self.pdb2_name
                    alt_name = self.pdb1_name
                else:
                    self._move_failed_full_outputs()
                    raise RuntimeError("Prediction with deep MSA was failed")
            else:
                if np.any(fs_tmscore[1, :] >= 0.5) and np.any(full_tmscore[1, :] >= 0.5):
                    ref_name = self.pdb2_name
                    alt_name = self.pdb1_name
                elif np.any(fs_tmscore[0, :] >= 0.5) and np.any(full_tmscore[0, :] >= 0.5):
                    ref_name = self.pdb1_name
                    alt_name = self.pdb2_name
                else:
                    self._move_failed_full_outputs()
                    raise RuntimeError("Prediction with deep MSA was failed")

            logger.info("Reference structure: %s", ref_name)
            logger.info("Alternative structure: %s", alt_name)

            np.savetxt(f"TMScore_full-MSA_{self.pdb1_name}.csv", full_tmscore, fmt="%2.3f")
            np.savetxt(f"TMScore_fs_full-MSA_{self.pdb1_name}.csv", fs_tmscore, fmt="%2.3f")

            tmscores_random: List[float] = []
            tmscores_fs_random: List[float] = []
            last_shallow: Optional[TMScore] = None

            for max_msa, ext_msa in self._msa_pairs():
                pred_dir = (
                    str(PREDICTIONS_ROOT / self.pdb1_name)
                    + f"/{self.pdb1_name}_predicted_models_rand_*"
                    + f"_max_{max_msa}_ext_{ext_msa}"
                )
                logger.debug("Shallow MSA dir pattern: %s", pred_dir)

                shallow = TMScore(
                    pred_dir,
                    self.pdb1,
                    self.pdb1_name,
                    self.pdb2,
                    self.pdb2_name,
                    self.model_type,
                )
                tmscores_random = list(np.append(tmscores_random, shallow.tmscores))
                last_shallow = shallow

                shallow_fs = TMScoreFSMulti(
                    pred_dir,
                    self.pdb1,
                    self.pdb1_name,
                    self.pdb2,
                    self.pdb2_name,
                )
                tmscores_fs_random = list(np.append(tmscores_fs_random, shallow_fs.tmscores_fs))

            random_array = np.asarray(tmscores_random, dtype=float)
            fs_random_array = np.asarray(tmscores_fs_random, dtype=float)

            tmscores_random_reshape = random_array.reshape(14, num_seeds * 5)
            tmscores_fs_random_reshape = fs_random_array.reshape(14, num_seeds * 5)

            if alt_name == self.pdb2_name:
                tmscores_random_alter = tmscores_random_reshape[1::2, :]
                tmscores_fs_random_alter = tmscores_fs_random_reshape[1::2, :]
            else:
                tmscores_random_alter = tmscores_random_reshape[0::2, :]
                tmscores_fs_random_alter = tmscores_fs_random_reshape[0::2, :]

            if np.any(tmscores_random_alter > 0.5) and np.any(tmscores_fs_random_alter > 0.5):
                np.savetxt(
                    f"TMScore_random-MSA_{self.pdb1_name}.csv",
                    tmscores_random_reshape,
                    fmt="%2.3f",
                )
                np.savetxt(
                    f"TMScore_fs_random-MSA_{self.pdb1_name}.csv",
                    tmscores_fs_random_reshape,
                    fmt="%2.3f",
                )

                logger.info("Finding optimal size of random MSA...")
                last_shallow.select_size(
                    tmscores_random_reshape,
                    tmscores_fs_random_reshape,
                    self.pdb1_name,
                    self.pdb2_name,
                    alt_name,
                    num_seeds,
                )
                self.size_selection = [last_shallow.selection]
                logger.info("Selected MSA size index: %s", self.size_selection)
            else:
                raise RuntimeError(
                    "Full-MSA prediction is not tightly aligned to crystal structure "
                    "with additional seeds"
                )
