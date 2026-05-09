# -*- coding: utf-8 -*-
"""Base classes and utilities for structural analysis."""

import logging
from pathlib import (
    Path,
)
from typing import (
    List,
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

logger = logging.getLogger(__name__)

ZERO_TM_SCORES = [0.0, 0.0, 0.0, 0.0, 0.0]
MULTIMER_MODEL_TYPE = "alphafold2_multimer_v3"


class BaseTMScore:
    """Base class for computing TM-scores for predicted protein models.

    This class provides common functionality for resolving model files
    and computing TM-scores against reference structures.
    """

    def __init__(
        self,
        pred_dir: str,
        pdb1: str,
        pdb1_name: str,
        pdb2: str,
        pdb2_name: str,
        model_type: str,
    ):
        self.pred_dir = Path(pred_dir)
        self.pdb1 = Path(pdb1)
        self.pdb2 = Path(pdb2)
        self.pdb1_name = pdb1_name
        self.pdb2_name = pdb2_name
        self.model_type = model_type
        self.tmscores: List[float] = self._calculate_scores()

    def _resolve_models(self) -> List[Path]:
        """Resolve paths to predicted model files."""
        if any(char in str(self.pred_dir) for char in "*?["):
            candidate_paths = sorted(Path().glob(str(self.pred_dir)))
            if len(candidate_paths) == 1 and candidate_paths[0].is_dir():
                return sorted(candidate_paths[0].glob("*_unrelaxed*pdb"))
            return [p for p in candidate_paths if p.suffix.lower() == ".pdb"]

        if self.model_type != MULTIMER_MODEL_TYPE:
            return sorted(self.pred_dir.glob("*_unrelaxed*pdb"))

        output_files = sorted(self.pred_dir.glob("rmTER*_unrelaxed*pdb"))
        if not output_files:
            self._convert_multimer_to_single()
            output_files = sorted(self.pred_dir.glob("rmTER*_unrelaxed*pdb"))
        return output_files

    def _convert_multimer_to_single(self):
        """Convert multimer predictions to single chains if needed."""
        from ..utils.convert_multi_single import (
            convert_m2s,
        )

        convert_m2s(str(self.pred_dir), self.pdb1_name, self.pdb2_name)

    def _calculate_scores(self) -> List[float]:
        """Calculate TM-scores against both reference structures."""
        predicted_models = self._resolve_models()
        if not predicted_models:
            logger.warning("No predicted models found for %s", self.pred_dir)
            return ZERO_TM_SCORES.copy()

        reference1 = get_structure(get_pdb_path(str(self.pdb1)))
        ref1_coords, ref1_seq = get_residue_data(reference1)
        reference2 = get_structure(get_pdb_path(str(self.pdb2)))
        ref2_coords, ref2_seq = get_residue_data(reference2)

        scores_with_ref1: List[float] = []
        scores_with_ref2: List[float] = []

        for model_path in predicted_models:
            model_path = model_path.with_suffix("")
            model_structure = get_structure(get_pdb_path(str(model_path)))
            model_coords, model_seq = get_residue_data(model_structure)

            alignment1 = tm_align(model_coords, ref1_coords, model_seq, ref1_seq)
            scores_with_ref1.append(round(alignment1.tm_norm_chain1, 5))

            alignment2 = tm_align(model_coords, ref2_coords, model_seq, ref2_seq)
            scores_with_ref2.append(round(alignment2.tm_norm_chain1, 5))

        if max(scores_with_ref1, default=0.0) >= max(scores_with_ref2, default=0.0):
            logger.debug("Using reference1 orientation for %s", self.pred_dir)
            return scores_with_ref1

        logger.debug("Using reference2 orientation for %s", self.pred_dir)
        return scores_with_ref2
