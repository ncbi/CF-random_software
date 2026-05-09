#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""TM-score analysis for fold-switching prediction workflows."""

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
from .cal_tmscore_fs_flmsa import (
    TMScoreFS,
)

PREDICTIONS_ROOT = Path("predictions_all")
FAILED_ROOT = Path("failed_prediction")


class TMScore(BaseTMScore):
    """Compute whole-structure TM-scores for a prediction set."""

    def select_size(
        self,
        tmscores_random: np.ndarray,
        tmscores_fs_random: np.ndarray,
        pdb1_name: str,
        pdb2_name: str,
        alt_name: str,
        num_seeds: int,
    ) -> None:
        if tmscores_random.size != 14 * num_seeds * 5:
            raise ValueError("Unexpected random TM-score matrix size")

        reshaped_random = tmscores_random.reshape(14, num_seeds * 5)
        reshaped_fs = tmscores_fs_random.reshape(14, num_seeds * 5)

        if alt_name == pdb2_name:
            alternative_indices = list(range(1, 14, 2))
            location_offset = 1
        else:
            alternative_indices = list(range(0, 14, 2))
            location_offset = 0

        alternative_matrix = reshaped_random[alternative_indices, :]
        best_group = int(np.argmax(alternative_matrix.sum(axis=1)))
        location = best_group * 2 + location_offset

        if not (
            np.any(reshaped_random[location, :] >= 0.5) and np.any(reshaped_fs[location, :] >= 0.5)
        ):
            raise RuntimeError("Predictions are bad: no viable shallow MSA selection")

        self.selection = int((location - location_offset) / 2)
        logger.info("Selected shallow MSA index %s for %s", self.selection, pdb1_name)


class TMScoreCalAllVarFS:
    """Evaluate TM-scores for fold-switching prediction results."""

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

    def _find_prediction_directory(self, pattern: str) -> Optional[Path]:
        candidates = sorted((PREDICTIONS_ROOT / self.pdb1_name).glob(pattern))
        return candidates[0] if candidates else None

    def _move_failed_full_outputs(self) -> None:
        target = FAILED_ROOT / self.pdb1_name
        target.mkdir(parents=True, exist_ok=True)
        for candidate in sorted(
            (PREDICTIONS_ROOT / self.pdb1_name).glob(
                f"{self.pdb1_name}_predicted_models_full_rand_*"
            )
        ):
            shutil.move(str(candidate), str(target / candidate.name))

    @staticmethod
    def _msa_pairs() -> List[tuple]:
        pairs: List[tuple] = []
        max_msa = 1
        ext_msa = 2
        for multiplier in (1, 2, 2, 2, 2, 2, 2):
            max_msa *= multiplier
            ext_msa *= multiplier
            pairs.append((max_msa, ext_msa))
        return pairs

    def _evaluate_monomer(self) -> None:
        num_seeds = 5 + self.nMSA
        full_dir = self._find_prediction_directory(f"{self.pdb1_name}_predicted_models_full_rand_*")
        if full_dir is None:
            raise FileNotFoundError(f"Full-MSA prediction directory not found for {self.pdb1_name}")

        full_scores = np.asarray(
            TMScore(
                str(full_dir),
                str(self.pdb1),
                self.pdb1_name,
                str(self.pdb2),
                self.pdb2_name,
                self.model_type,
            ).tmscores,
            dtype=float,
        ).reshape(2, num_seeds * 5)

        if np.all(full_scores < 0.5):
            self._move_failed_full_outputs()
            raise RuntimeError("All predictions with deep MSA are failed")

        if np.average(full_scores[0, :]) >= np.average(full_scores[1, :]):
            alt_name = self.pdb2_name
        else:
            alt_name = self.pdb1_name

        np.savetxt(f"TMScore_full-MSA_{self.pdb1_name}.csv", full_scores, fmt="%2.3f")

        random_scores: List[float] = []
        fs_random_scores: List[float] = []
        for max_msa, ext_msa in self._msa_pairs():
            pattern = f"{self.pdb1_name}_predicted_models_rand_*_max_{max_msa}_ext_{ext_msa}"
            pred_path = PREDICTIONS_ROOT / self.pdb1_name / pattern
            random_scores.extend(
                TMScore(
                    str(pred_path),
                    str(self.pdb1),
                    self.pdb1_name,
                    str(self.pdb2),
                    self.pdb2_name,
                    self.model_type,
                ).tmscores
            )
            fs_random_scores.extend(
                TMScoreFS(
                    str(pred_path), str(self.pdb1), self.pdb1_name, self.pdb2_name, self.pdb2_name
                ).tmscores_fs
            )

        random_array = np.asarray(random_scores, dtype=float)
        fs_random_array = np.asarray(fs_random_scores, dtype=float)
        if random_array.size != 14 * num_seeds * 5 or fs_random_array.size != 14 * num_seeds * 5:
            raise RuntimeError("Unexpected shallow MSA TM-score matrix size")

        selector = TMScore(
            str(pred_path),
            str(self.pdb1),
            self.pdb1_name,
            str(self.pdb2),
            self.pdb2_name,
            self.model_type,
        )
        selector.select_size(
            random_array, fs_random_array, self.pdb1_name, self.pdb2_name, alt_name, num_seeds
        )
        self.size_selection = [selector.selection]

        np.savetxt(
            f"TMScore_random-MSA_{self.pdb1_name}.csv",
            random_array.reshape(14, num_seeds * 5),
            fmt="%2.3f",
        )
        np.savetxt(
            f"TMScore_fs_random-MSA_{self.pdb1_name}.csv",
            fs_random_array.reshape(14, num_seeds * 5),
            fmt="%2.3f",
        )

    def _evaluate_multimer(self) -> None:
        if not self.search_dir and not self.search_multi_dir:
            raise ValueError(
                "Multimer fold-switching evaluation requires search_dir or search_multi_dir"
            )

        PredictionAllMultimerFS(
            self.pdb1_name,
            self.pdb2_name,
            self.search_dir,
            self.nMSA,
            self.model_type,
            self.search_multi_dir,
            str(self.pdb1),
            str(self.pdb2),
        )
        self.size_selection = []
