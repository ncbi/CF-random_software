#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""TM-score analysis for alternative conformation prediction workflows."""

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

from .base import (
    MULTIMER_MODEL_TYPE,
    BaseTMScore,
    logger,
)

PREDICTIONS_ROOT = Path("predictions_all")
FAILED_ROOT = Path("failed_prediction")
MSA_MULTIPLIERS = (1, 2, 2, 2, 2, 2, 2)


class TMScore(BaseTMScore):
    """Compute TM-scores for one prediction directory.

    Inherits all scoring logic from BaseTMScore (forward/reverse comparison,
    glob resolution, multimer conversion). Only adds select_size, which
    identifies the optimal shallow MSA depth index from clustered TM-scores.
    """

    def select_size(
        self,
        tmscores_random: np.ndarray,
        pdb1_name: str,
        pdb2_name: str,
        alt_name: str,
        num_seeds: int,
    ) -> None:
        """Select the optimal shallow MSA depth index.

        Mirrors the original select_size logic exactly:
            1. Reshape flat array to (14, num_seeds*5).
            2. Extract every-other row for the alternative structure.
            3. Sum each group row and find the argmax.
            4. Map back to the full-matrix row via alt_name offset.
            5. Verify at least one score >= 0.5 threshold.

        Args:
            tmscores_random: Flat array of all shallow-MSA TM-scores.
            pdb1_name: Name of first reference structure.
            pdb2_name: Name of second reference structure.
            alt_name: Name of the alternative conformation structure.
            num_seeds: Number of prediction seeds.

        Raises:
            RuntimeError: If no model exceeds the 0.5 TM-score threshold.
        """
        tmscores_random_reshape = tmscores_random.reshape(14, num_seeds * 5)
        tmscores_random_locat = np.zeros((7, num_seeds * 5))

        # Extract rows for the alternative structure
        if alt_name == pdb2_name:
            tmp_cnt = 0
            for i in range(1, 14, 2):
                tmscores_random_locat[tmp_cnt, :] = tmscores_random_reshape[i, :]
                tmp_cnt += 1
        else:
            tmp_cnt = 0
            for i in range(0, 13, 2):
                tmscores_random_locat[tmp_cnt, :] = tmscores_random_reshape[i, :]
                tmp_cnt += 1

        # Sum each MSA-depth group and pick the best via argmax of max
        tmscore_data_sum = np.zeros((7, 1))
        for i in range(tmscores_random_locat.shape[0]):
            tmscore_data_sum[i] = np.sum(tmscores_random_locat[i])

        location = int(np.argmax(np.max(tmscore_data_sum, axis=1)))

        # Map group index back to full-matrix row using alt_name offset
        if alt_name == pdb2_name:
            location_full = (location * 2) + 1
        else:
            location_full = location * 2

        tmscore_check = tmscores_random_reshape

        if alt_name == pdb2_name and np.any(tmscore_check[location_full, :] >= 0.5):
            self.selection = int((location_full - 1) / 2)
        elif alt_name == pdb1_name and np.any(tmscore_check[location_full, :] >= 0.5):
            self.selection = int(location_full / 2)
        else:
            raise RuntimeError(
                "Predictions are bad: no alternative-conformation models exceed "
                "the 0.5 TM-score threshold"
            )

        logger.info("Selected shallow MSA index %s for %s", self.selection, pdb1_name)


class TMScoreCalAllVar:
    """Evaluate TM-scores for alternative-conformation prediction results."""

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

    def _evaluate_monomer(self) -> None:
        """Run the full monomer TM-score evaluation pipeline."""
        num_seeds = 5 + self.num_msa

        # Pass as glob pattern string — BaseTMScore._resolve_models expands it
        full_pred_dir = (
            str(PREDICTIONS_ROOT / self.pdb1_name)
            + f"/{self.pdb1_name}_predicted_models_full_rand_*"
        )
        msa_full = TMScore(
            full_pred_dir,
            self.pdb1,
            self.pdb1_name,
            self.pdb2,
            self.pdb2_name,
            self.model_type,
        )
        full_scores_array = np.asarray(msa_full.tmscores, dtype=float).reshape(2, num_seeds * 5)

        # Three-branch quality check
        if np.any(full_scores_array[0, :] > 0.5) or np.any(full_scores_array[1, :] > 0.5):
            alt_name = self._determine_alternative(np.average(full_scores_array, axis=1))
        elif np.all(full_scores_array[0, :] < 0.5) and np.all(full_scores_array[1, :] < 0.5):
            self._move_failed_full_outputs()
            raise RuntimeError("All predictions with deep MSA are failed")
        else:
            alt_name = self._determine_alternative(np.average(full_scores_array, axis=1))

        logger.info(
            "Reference: %s  Alternative: %s",
            self.pdb1_name if alt_name == self.pdb2_name else self.pdb2_name,
            alt_name,
        )
        np.savetxt(f"TMScore_full-MSA_{self.pdb1_name}.csv", full_scores_array, fmt="%2.3f")

        # Shallow random MSA TM-scores
        max_msa = 1
        ext_msa = 2
        tmscores_random: List[float] = []
        last_shallow: Optional[TMScore] = None

        for multi in MSA_MULTIPLIERS:
            max_msa *= multi
            ext_msa *= multi

            pred_dir = (
                str(PREDICTIONS_ROOT / self.pdb1_name)
                + f"/{self.pdb1_name}_predicted_models_rand_*"
                + f"_max_{max_msa}_ext_{ext_msa}/"
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

        random_array = np.asarray(tmscores_random, dtype=float)
        tmscores_random_reshape = random_array.reshape(14, num_seeds * 5)

        # Extract alternative rows for quality check
        tmscores_random_alter = self._extract_alternative_rows(
            tmscores_random_reshape, alt_name, self.pdb1_name, self.pdb2_name
        )
        if np.all(tmscores_random_alter < 0.5):
            raise RuntimeError("All shallow predictions are failed")

        logger.info("Finding optimal size of random MSA...")
        last_shallow.select_size(random_array, self.pdb1_name, self.pdb2_name, alt_name, num_seeds)
        self.size_selection = [last_shallow.selection]
        logger.info("Selected MSA size index: %s", self.size_selection)

        np.savetxt(f"TMScore_random-MSA_{self.pdb1_name}.csv", tmscores_random_reshape, fmt="%2.3f")

    def _determine_alternative(self, reference_scores: np.ndarray) -> str:
        """Return the name of the alternative conformation structure."""
        if reference_scores[0] >= reference_scores[1]:
            return self.pdb2_name
        return self.pdb1_name

    @staticmethod
    def _extract_alternative_rows(
        matrix: np.ndarray,
        alt_name: str,
        pdb1_name: str,
        pdb2_name: str,
    ) -> np.ndarray:
        """Extract rows corresponding to the alternative structure."""
        if alt_name == pdb2_name:
            return matrix[1::2, :]
        return matrix[0::2, :]

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

    def _evaluate_multimer(self) -> None:
        """Run multimer TM-score evaluation pipeline."""
