#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Base classes and utilities for structural analysis."""

import glob
import logging
import os
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

    For each predicted model, TM-scores are computed in both the forward
    (model→reference) and reverse (reference→model) directions against
    both reference structures. Whichever direction yields the higher maximum
    score is kept.
    """

    def __init__(
        self,
        pred_dir: str,
        pdb1: str,
        pdb1_name: str,
        pdb2: str,
        pdb2_name: str,
        model_type: str,
    ) -> None:
        """Initialise and compute TM-scores.

        Args:
            pred_dir: Glob pattern or path to prediction directory.
            pdb1: Path to first reference PDB (without extension).
            pdb1_name: Name of first reference structure.
            pdb2: Path to second reference PDB (without extension).
            pdb2_name: Name of second reference structure.
            model_type: ColabFold model type string.
        """
        self.pred_dir = pred_dir
        self.pdb1 = pdb1
        self.pdb1_name = pdb1_name
        self.pdb2 = pdb2
        self.pdb2_name = pdb2_name
        self.model_type = model_type
        self.tmscores: List[float] = self._calculate_scores()

    def _resolve_models(self) -> List[str]:
        """Resolve paths to predicted model PDB files.

        Handles both plain directory paths and glob patterns.
        For multimer models, converts to single-chain first if needed.

        Returns:
            List of absolute PDB file path strings (without extension).
        """
        pwd = os.getcwd() + "/"

        if self.model_type != MULTIMER_MODEL_TYPE:
            pattern = str(self.pred_dir) + "/*_unrelaxed*pdb"
            files_list = glob.glob(pattern)
        else:
            check_pattern = str(self.pred_dir) + "/rmTER*_unrelaxed*pdb"
            files_list = glob.glob(check_pattern)
            if not files_list:
                self._convert_multimer_to_single()
                files_list = glob.glob(check_pattern)

        # Prepend pwd and strip extension
        return [pwd + f.replace(".pdb", "") for f in files_list]

    def _convert_multimer_to_single(self) -> None:
        """Convert multimer predictions to single chains."""
        from ..utils.convert_multi_single import (
            convert_m2s,
        )

        convert_m2s(str(self.pred_dir), self.pdb1_name, self.pdb2_name)

    def _calculate_scores(self) -> List[float]:
        """Calculate TM-scores against both reference structures.

        Computes scores in both forward (model→ref) and reverse (ref→model)
        directions. Whichever direction yields the higher max is returned.

        Returns:
            List of TM-scores, one per predicted model, for both references
            concatenated (pdb1 scores first, then pdb2 scores).
        """
        pwd = os.getcwd() + "/"
        predicted_models = self._resolve_models()

        if not predicted_models:
            logger.warning("No predicted models found for %s", self.pred_dir)
            return ZERO_TM_SCORES.copy()

        # Load reference structures
        ref1 = get_structure(get_pdb_path(pwd + self.pdb1_name))
        ref1_coords, ref1_seq = get_residue_data(ref1)

        ref2 = get_structure(get_pdb_path(pwd + self.pdb2_name))
        ref2_coords, ref2_seq = get_residue_data(ref2)

        tmscores_ord: List[float] = []
        tmscores_rev: List[float] = []

        for model_path in predicted_models:
            s = get_structure(get_pdb_path(model_path))
            coords, seq = get_residue_data(s)

            res = tm_align(coords, ref1_coords, seq, ref1_seq)
            tmscores_ord.append(round(res.tm_norm_chain1, 5))

            res = tm_align(ref1_coords, coords, ref1_seq, seq)
            tmscores_rev.append(round(res.tm_norm_chain1, 5))

        for model_path in predicted_models:
            s = get_structure(get_pdb_path(model_path))
            coords, seq = get_residue_data(s)

            res = tm_align(coords, ref2_coords, seq, ref2_seq)
            tmscores_ord.append(round(res.tm_norm_chain1, 5))

            res = tm_align(ref2_coords, coords, ref2_seq, seq)
            tmscores_rev.append(round(res.tm_norm_chain1, 5))

        logger.debug("TM-scores forward: %s", tmscores_ord)
        logger.debug("TM-scores reverse: %s", tmscores_rev)

        # Return whichever direction gives the higher maximum
        if np.max(tmscores_ord) > np.max(tmscores_rev):
            return tmscores_ord
        return tmscores_rev
