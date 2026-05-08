#!/bin/env python3
# -*- coding: utf-8 -*-
"""Helpers to run ColabFold prediction batches across varied MSAs.

Small utilities used by higher-level workflows to orchestrate
multiple ColabFold runs with different MSA settings.
"""

import logging
from pathlib import (
    Path,
)

import numpy as np

from .base import (
    MSAMaxRunner,
    MSAVariableRunner,
)


class PredictionAll:
    """High-level orchestration for full and variable MSA predictions."""

    def __init__(
        self,
        pdb1_name: str,
        search_dir: str,
        search_multi_dir: str,
        nMSA: int,
        model_type: str,
    ) -> None:
        """Run the full and varied MSA prediction pipeline.

        Args:
            pdb1_name: Name of the target protein.
            search_dir: Path to the single-chain MSA folder.
            search_multi_dir: Path to the multimer MSA folder.
            nMSA: Number of additional MSA seeds to add to the default 5.
            model_type: ColabFold model type.
        """
        self.pdb1_name = pdb1_name
        self.search_dir = search_dir
        self.search_multi_dir = search_multi_dir
        self.model_type = model_type

        self.base_output_dir = Path("predictions_all") / pdb1_name
        self.base_output_dir.mkdir(parents=True, exist_ok=True)

        num_seeds = nMSA + 5
        random_seed = int(np.random.randint(0, 100, 1))

        full_output_dir = (
            self.base_output_dir / f"{pdb1_name}_predicted_models_full_rand_{random_seed}"
        )
        logging.info("Running full MSA prediction: %s", full_output_dir)
        MSAMaxRunner(
            self.search_dir,
            str(full_output_dir),
            self.pdb1_name,
            random_seed,
            num_seeds,
            self.model_type,
        )

        if self.model_type == "alphafold2_multimer_v3":
            variable_search_dir = self.search_multi_dir
        else:
            variable_search_dir = self.search_dir

        variable_output_dir = self.base_output_dir / f"{pdb1_name}_predicted_models_rand_"
        logging.info("Running variable MSA predictions under: %s", variable_output_dir)
        MSAVariableRunner(
            variable_search_dir,
            str(variable_output_dir),
            self.pdb1_name,
            random_seed,
            num_seeds,
            self.model_type,
        )
