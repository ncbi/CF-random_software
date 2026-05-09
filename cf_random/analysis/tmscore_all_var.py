#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""TM-score analysis for alternative conformation prediction workflows."""

import shutil
from pathlib import (
    Path,
)
from typing import (
    List,
    Optional,
)

import numpy as np
from tmtools import (
    tm_align,
)
from tmtools.io import (
    get_residue_data,
    get_structure,
)
from tmtools.testing import (
    get_pdb_path,
)

from ..prediction.pred_cal_tmscore_multimer import (
    PredictionAllMultimerFS,
)
from .base import (
    MULTIMER_MODEL_TYPE,
    BaseTMScore,
    logger,
)

PREDICTIONS_ROOT = Path("predictions_all")
FAILED_ROOT = Path("failed_prediction")


class TMScore(BaseTMScore):
    """Compute TM-scores for one prediction directory."""

    def select_size(
        self,
        tmscores_random: np.ndarray,
        pdb1_name: str,
        pdb2_name: str,
        alt_name: str,
        num_seeds: int,
    ) -> None:
        if tmscores_random.size != 14 * num_seeds * 5:
            raise ValueError("Unexpected TM-score matrix size for selection")

        reshaped = tmscores_random.reshape(14, num_seeds * 5)
        if alt_name == pdb2_name:
            alternative_rows = list(range(1, 14, 2))
            offset = 1
        else:
            alternative_rows = list(range(0, 14, 2))
            offset = 0

        grouped = reshaped[alternative_rows, :]
        best_group = int(np.argmax(grouped.sum(axis=1)))
        location = best_group * 2 + offset
        self.selection = int((location - offset) / 2)

        if not np.any(reshaped[location, :] >= 0.5):
            logger.error("No acceptable TM-score found for selection at location %s", location)
            raise RuntimeError(
                "Predictions are bad: no alternative-conformation models exceed TM-score threshold"
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
        nMSA: int,
        option: str,
        model_type: str,
        search_dir: Optional[str] = None,
        search_multi_dir: Optional[str] = None,
    ) -> None:
        self.pdb1 = Path(pdb1)
        self.pdb2 = Path(pdb2)
        self.pdb1_name = pdb1_name
        self.pdb2_name = pdb2_name
        self.nMSA = nMSA
        self.option = option
        self.model_type = model_type
        self.search_dir = search_dir
        self.search_multi_dir = search_multi_dir
        self.size_selection: List[int] = []

        if self.model_type != MULTIMER_MODEL_TYPE:
            self._evaluate_monomer()
        else:
            self._evaluate_multimer()

    def _find_single_output(self, pattern: str) -> Optional[Path]:
        candidates = sorted((PREDICTIONS_ROOT / self.pdb1_name).glob(pattern))
        return candidates[0] if candidates else None

    def _move_failed_full_outputs(self) -> None:
        FAILED_ROOT.mkdir(parents=True, exist_ok=True)
        for candidate in sorted(
            (PREDICTIONS_ROOT / self.pdb1_name).glob(
                f"*{self.pdb1_name}_predicted_models_full_rand_*"
            )
        ):
            destination = FAILED_ROOT / self.pdb1_name / candidate.name
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.move(str(candidate), str(destination))

    def _evaluate_monomer(self) -> None:
        num_seeds = 5 + self.nMSA
        full_pattern = f"{self.pdb1_name}_predicted_models_full_rand_*"
        full_dir = self._find_single_output(full_pattern)
        if full_dir is None:
            raise FileNotFoundError(f"Full-MSA prediction directory not found for {self.pdb1_name}")

        full_scores = TMScore(
            str(full_dir),
            str(self.pdb1),
            self.pdb1_name,
            str(self.pdb2),
            self.pdb2_name,
            self.model_type,
        ).tmscores
        full_scores_array = np.asarray(full_scores, dtype=float).reshape(2, num_seeds * 5)

        reference_scores = np.average(full_scores_array, axis=1)
        if np.all(full_scores_array < 0.5):
            self._move_failed_full_outputs()
            raise RuntimeError("All predictions with deep MSA are failed")

        alt_name = self._determine_alternative(reference_scores)
        logger.info(
            "Reference structure: %s, alternative structure: %s",
            self.pdb1_name if alt_name == self.pdb2_name else self.pdb2_name,
            alt_name,
        )

        np.savetxt(f"TMScore_full-MSA_{self.pdb1_name}.csv", full_scores_array, fmt="%2.3f")

        tmscores_random: List[float] = []
        for max_msa, ext_msa in self._msa_pairs():
            pattern = f"{self.pdb1_name}_predicted_models_rand_*_max_{max_msa}_ext_{ext_msa}"
            shallow_dir = PREDICTIONS_ROOT / self.pdb1_name / pattern
            tmscores_random.extend(
                TMScore(
                    str(shallow_dir),
                    str(self.pdb1),
                    self.pdb1_name,
                    str(self.pdb2),
                    self.pdb2_name,
                    self.model_type,
                ).tmscores
            )

        random_array = np.asarray(tmscores_random, dtype=float)
        if random_array.size != 14 * num_seeds * 5:
            raise RuntimeError("Unexpected shallow MSA TM-score count")

        random_matrix = random_array.reshape(14, num_seeds * 5)
        random_alter = self._extract_alternative_rows(
            random_matrix, alt_name, self.pdb1_name, self.pdb2_name
        )
        if np.all(random_alter < 0.5):
            raise RuntimeError("All shallow predictions are failed")

        selector = TMScore(
            str(shallow_dir),
            str(self.pdb1),
            self.pdb1_name,
            str(self.pdb2),
            self.pdb2_name,
            self.model_type,
        )
        selector.select_size(random_array, self.pdb1_name, self.pdb2_name, alt_name, num_seeds)
        self.size_selection = [selector.selection]
        np.savetxt(f"TMScore_random-MSA_{self.pdb1_name}.csv", random_matrix, fmt="%2.3f")

    def _determine_alternative(self, reference_scores: np.ndarray) -> str:
        if reference_scores[0] >= reference_scores[1]:
            return self.pdb2_name
        return self.pdb1_name

    @staticmethod
    def _msa_pairs() -> List[tuple]:
        pairs = []
        max_msa = 1
        ext_msa = 2
        for multiplier in (1, 2, 2, 2, 2, 2, 2):
            max_msa *= multiplier
            ext_msa *= multiplier
            pairs.append((max_msa, ext_msa))
        return pairs

    @staticmethod
    def _extract_alternative_rows(
        matrix: np.ndarray, alt_name: str, pdb1_name: str, pdb2_name: str
    ) -> np.ndarray:
        if alt_name == pdb2_name:
            return matrix[1::2, :]
        return matrix[0::2, :]

    def _evaluate_multimer(self) -> None:
        if not self.search_dir and not self.search_multi_dir:
            raise ValueError("Multimer evaluation requires search_dir or search_multi_dir")

        PredictionAllMultimerFS(
            self.pdb1_name,
            self.pdb2_name,
            self.search_dir,
            self.nMSA,
            self.model_type,
            self.search_multi_dir,
        )
        self.size_selection = []
