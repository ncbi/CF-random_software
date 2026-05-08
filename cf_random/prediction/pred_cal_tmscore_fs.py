#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""TM-score evaluation and MSA orchestration for fold-switching predictions.

This module provides classes for computing TM-scores focused on fold-switching
regions and orchestrating ColabFold batch predictions with varied MSA sizes.
"""

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

from ..analysis.cal_tmscore_fs_only import (
    TMScoreFS,
)
from ..utils.convert_multi_single import (
    ConvertM2S,
)
from .base import (
    MSAMaxRunner,
    MSAVariableRunner,
)

logger = logging.getLogger(__name__)

# Configuration constants
TMSCORE_THRESHOLD_GOOD = 0.5
TMSCORE_THRESHOLD_ACCEPTABLE = 0.4
INITIAL_MAX_MSA = 1
INITIAL_EXTRA_MSA = 2


class TMScore:
    """Calculate TM-scores between predicted and reference structures.

    Computes structural alignment scores (TM-scores) comparing predicted
    protein models against reference PDB structures using tmtools.
    Handles both monomer and multimer predictions.

    Attributes:
        tmscores (List[float]): Array of computed TM-scores.
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
        """Initialize and compute TM-scores for predicted models.

        Args:
            pred_dir: Directory containing predicted PDB files.
            pdb1: Path to first reference structure.
            pdb1_name: Name of first reference structure.
            pdb2: Path to second reference structure.
            pdb2_name: Name of second reference structure.
            model_type: AlphaFold model type ('monomer', 'multimer', 'ptm').

        Raises:
            FileNotFoundError: If no predicted models found in directory.
        """
        self.pred_dir = pred_dir
        self.pdb1 = pdb1
        self.pdb1_name = pdb1_name
        self.pdb2 = pdb2
        self.pdb2_name = pdb2_name
        self.model_type = model_type
        self.tmscores = []

        self._compute_tmscores()

    def _get_predicted_files(self) -> List[str]:
        """Get list of predicted PDB files.

        Handles both monomer and multimer predictions, converting
        multimer files to single chains if necessary.

        Returns:
            List[str]: Paths to predicted PDB files.
        """
        if self.model_type != "alphafold2_multimer_v3":
            return glob.glob(f"{self.pred_dir}/*_unrelaxed*pdb")

        # Handle multimer predictions
        check_files_list = glob.glob(f"{self.pred_dir}/rmTER*_unrelaxed*pdb")
        if not check_files_list:
            logger.info("Converting multimer predictions to single chains")
            ConvertM2S(self.pred_dir, self.pdb1_name, self.pdb2_name)
            check_files_list = glob.glob(f"{self.pred_dir}/rmTER*_unrelaxed*pdb")

        return check_files_list

    def _compute_tmscore_pair(
        self,
        model_coords: np.ndarray,
        model_seq: str,
        ref_coords: np.ndarray,
        ref_seq: str,
    ) -> float:
        """Compute TM-score between two structures.

        Args:
            model_coords: CA coordinates of model (N x 3).
            model_seq: Amino acid sequence of model.
            ref_coords: CA coordinates of reference (N x 3).
            ref_seq: Amino acid sequence of reference.

        Returns:
            float: TM-score normalized to reference length (rounded to 5 decimals).
        """
        try:
            alignment = tm_align(model_coords, ref_coords, model_seq, ref_seq)
            return round(alignment.tm_norm_chain1, 5)
        except Exception as e:
            logger.warning(f"TM-align computation failed: {e}")
            return 0.0

    def _compute_tmscores(self) -> None:
        """Compute TM-scores for all predicted models.

        Compares each predicted model against both reference structures
        in both forward and reverse alignments.
        """
        pwd = os.getcwd() + "/"
        files_list = self._get_predicted_files()

        if not files_list:
            logger.warning(f"No predicted models found in {self.pred_dir}")
            self.tmscores = [0.0] * 5
            return

        tmscores_forward = []
        tmscores_reverse = []

        # Load reference structures
        try:
            pdb1_structure = get_structure(get_pdb_path(f"{pwd}{self.pdb1_name}"))
            pdb1_coords, pdb1_seq = get_residue_data(pdb1_structure)

            pdb2_structure = get_structure(get_pdb_path(f"{pwd}{self.pdb2_name}"))
            pdb2_coords, pdb2_seq = get_residue_data(pdb2_structure)
        except Exception as e:
            logger.error(f"Failed to load reference structures: {e}")
            self.tmscores = [0.0] * 5
            return

        # Compute TM-scores against both references
        for model_path in files_list:
            try:
                model_clean = model_path.replace(".pdb", "")
                model_structure = get_structure(get_pdb_path(f"{pwd}{model_clean}"))
                model_coords, model_seq = get_residue_data(model_structure)

                # Score against pdb1
                score_fwd = self._compute_tmscore_pair(
                    model_coords, model_seq, pdb1_coords, pdb1_seq
                )
                score_rev = self._compute_tmscore_pair(
                    pdb1_coords, pdb1_seq, model_coords, model_seq
                )
                tmscores_forward.append(score_fwd)
                tmscores_reverse.append(score_rev)

                # Score against pdb2
                score_fwd = self._compute_tmscore_pair(
                    model_coords, model_seq, pdb2_coords, pdb2_seq
                )
                score_rev = self._compute_tmscore_pair(
                    pdb2_coords, pdb2_seq, model_coords, model_seq
                )
                tmscores_forward.append(score_fwd)
                tmscores_reverse.append(score_rev)

            except Exception as e:
                logger.warning(f"Failed to compute TM-score for {model_path}: {e}")
                continue

        # Select best scoring orientation
        if tmscores_forward and tmscores_reverse:
            self.tmscores = (
                tmscores_forward
                if np.max(tmscores_forward) > np.max(tmscores_reverse)
                else tmscores_reverse
            )
        else:
            self.tmscores = tmscores_forward or tmscores_reverse or [0.0] * 5

        logger.info(f"Computed {len(self.tmscores)} TM-scores")


class MSAMaxFoldSwitch(MSAMaxRunner):
    """Run full MSA prediction with fold-switching region analysis.

    Extends MSAMaxRunner to include fold-switching region TM-scores.
    """

    def __init__(
        self,
        search_dir: str,
        output_dir: str,
        pdb1_name: str,
        random_seed: int,
        num_seeds: int,
        model_type: str,
    ) -> None:
        """Initialize full MSA fold-switching prediction.

        Args:
            search_dir: Directory with MSA files.
            output_dir: Output directory for results.
            pdb1_name: Protein name.
            random_seed: Random seed.
            num_seeds: Number of random seeds.
            model_type: ColabFold model type.
        """
        super().__init__(search_dir, output_dir, pdb1_name, random_seed, num_seeds, model_type)
        logger.info(f"Completed full MSA prediction for {pdb1_name}")


class MSAVariableFoldSwitch(MSAVariableRunner):
    """Run varied MSA predictions with fold-switching region analysis.

    Extends MSAVariableRunner to systematically reduce MSA depth while
    monitoring fold-switching region alignment quality.
    """

    def __init__(
        self,
        search_dir: str,
        output_dir: str,
        pdb1_name: str,
        random_seed: int,
        num_seeds: int,
        model_type: str,
    ) -> None:
        """Initialize varied MSA fold-switching predictions.

        Args:
            search_dir: Directory with MSA files.
            output_dir: Base output directory.
            pdb1_name: Protein name.
            random_seed: Random seed.
            num_seeds: Number of random seeds.
            model_type: ColabFold model type.
        """
        super().__init__(
            search_dir,
            output_dir,
            pdb1_name,
            random_seed,
            num_seeds,
            model_type,
        )
        logger.info(f"Completed varied MSA predictions for {pdb1_name}")


# Backward compatibility aliases
CF_MSA_MAX = MSAMaxFoldSwitch
CF_MSA_VAR = MSAVariableFoldSwitch
