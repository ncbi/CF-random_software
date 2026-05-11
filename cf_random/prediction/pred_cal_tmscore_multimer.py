#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Multimer-aware TM-score utilities and ColabFold orchestration.

Provides TM-score computations tailored for multimer predictions and
helpers to run ColabFold batches for monomer/multimer cases.
"""

import glob
import logging
import os
from pathlib import (
    Path,
)
from typing import (
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

from .base import (
    MSAMaxRunner,
    MSAVariableRunner,
)

logger = logging.getLogger(__name__)

# Configuration constants
TMSCORE_THRESHOLD = 0.5
INITIAL_MAX_MSA = 1
INITIAL_EXTRA_MSA = 2


class TMScoreMonomer:
    """Calculate TM-scores for monomer chains in multimer predictions.

    Computes TM-align scores for individual chains extracted from
    multimer prediction results.

    Attributes:
        tmscores_monomer (List[float]): Array of computed TM-scores.
    """

    def __init__(
        self,
        pred_dir: str,
        pdb1_name: str,
        pdb2_name: str,
    ) -> None:
        """Initialize monomer TM-score calculation.

        Args:
            pred_dir: Directory containing predicted PDB files.
            pdb1_name: Name of first reference structure.
            pdb2_name: Name of second reference structure.
        """
        self.pred_dir = pred_dir
        self.pdb1_name = pdb1_name
        self.pdb2_name = pdb2_name
        self.tmscores_monomer = []
        self._compute_scores()

    def _compute_scores(self) -> None:
        """Compute TM-scores for monomer chains."""
        pwd = os.getcwd() + "/"
        files_list = glob.glob(f"{self.pred_dir}/*_unrelaxed*pdb")

        if not files_list:
            logger.warning(f"No predicted models found in {self.pred_dir}")
            self.tmscores_monomer = [0.0] * 5
            return

        try:
            pdb1_structure = get_structure(get_pdb_path(f"{pwd}{self.pdb1_name}"))
            pdb1_coords, pdb1_seq = get_residue_data(pdb1_structure)
        except Exception as e:
            logger.error(f"Failed to load reference structure {self.pdb1_name}: {e}")
            self.tmscores_monomer = [0.0] * 5
            return

        for model_path in files_list:
            try:
                model_clean = model_path.replace(".pdb", "")
                model_structure = get_structure(get_pdb_path(f"{pwd}{model_clean}"))
                model_coords, model_seq = get_residue_data(model_structure)

                alignment = tm_align(model_coords, pdb1_coords, model_seq, pdb1_seq)
                tmscore = round(alignment.tm_norm_chain1, 5)
                self.tmscores_monomer.append(tmscore)
            except Exception as e:
                logger.warning(f"Failed to compute TM-score for {model_path}: {e}")
                continue

        logger.info(f"Computed {len(self.tmscores_monomer)} monomer TM-scores")


class MSAMaxMultimer(MSAMaxRunner):
    """Run full MSA prediction for multimer structures.

    Extends MSAMaxRunner for multimer-specific predictions.
    """

    def __init__(
        self,
        search_dir: str,
        output_dir: str,
        pdb_name: str,
        random_seed: int,
        num_seeds: int,
        model_type: str,
    ) -> None:
        """Initialize full MSA multimer prediction.

        Args:
            search_dir: Directory with MSA files.
            output_dir: Output directory for results.
            pdb_name: Protein name.
            random_seed: Random seed.
            num_seeds: Number of random seeds.
            model_type: ColabFold model type.
        """
        super().__init__(search_dir, output_dir, pdb_name, random_seed, num_seeds, model_type)
        logger.info(f"Completed full MSA multimer prediction for {pdb_name}")


class MSAVariableMultimer(MSAVariableRunner):
    """Run varied MSA predictions for multimer structures.

    Systematically reduces MSA depth while monitoring prediction quality
    for multimer complexes.
    """

    def __init__(
        self,
        search_dir: str,
        output_dir: str,
        pdb1_name: str,
        pdb2_name: str,
        random_seed: int,
        num_seeds: int,
        model_type: str,
    ) -> None:
        """Initialize varied MSA multimer predictions.

        Args:
            search_dir: Directory with MSA files.
            output_dir: Base output directory.
            pdb1_name: First protein name.
            pdb2_name: Second protein name.
            random_seed: Random seed.
            num_seeds: Number of random seeds.
            model_type: ColabFold model type.
        """
        self.pdb1_name = pdb1_name
        self.pdb2_name = pdb2_name
        super().__init__(
            search_dir,
            output_dir,
            pdb1_name,
            random_seed,
            num_seeds,
            model_type,
        )
        logger.info(f"Completed varied MSA multimer predictions for {pdb1_name}")


class PredictionAllMultimerFS:
    """Compatibility wrapper for multimer prediction workflows.

    This helper exists to preserve the legacy analysis entry point used by
    `cf_random.analysis` while delegating to the new multimer prediction runners.
    """

    def __init__(
        self,
        pdb1_name: str,
        pdb2_name: str,
        search_dir: str,
        num_msa: int,
        model_type: str,
        search_multi_dir: Optional[str] = None,
        pdb1: Optional[str] = None,
        pdb2: Optional[str] = None,
    ) -> None:
        if model_type != "alphafold2_multimer_v3":
            raise ValueError("PredictionAllMultimerFS only supports multimer model_type")

        self.pdb1_name = pdb1_name
        self.pdb2_name = pdb2_name
        self.search_dir = search_multi_dir or search_dir
        self.num_seeds = 5 + num_msa
        self.random_seed = int(np.random.randint(0, 100, 1))
        self.size_selection = []

        output_dir = Path("predictions_all") / pdb1_name
        output_dir.mkdir(parents=True, exist_ok=True)

        full_output_dir = output_dir / f"{pdb1_name}_predicted_models_full_rand_{self.random_seed}"
        MSAMaxMultimer(
            self.search_dir,
            str(full_output_dir),
            pdb1_name,
            self.random_seed,
            self.num_seeds,
            model_type,
        )

        variable_output_dir = output_dir / f"{pdb1_name}_predicted_models_rand_"
        MSAVariableMultimer(
            self.search_dir,
            str(variable_output_dir),
            pdb1_name,
            pdb2_name,
            self.random_seed,
            self.num_seeds,
            model_type,
        )
