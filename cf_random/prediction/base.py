#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Base classes and utilities for ColabFold prediction orchestration.

This module provides the foundation for running ColabFold predictions with
various MSA configurations. It consolidates common patterns for managing
prediction workflows and result organization.
"""

import glob
import logging
import os
import shutil
from abc import (
    ABC,
)
from pathlib import (
    Path,
)
from typing import (
    Optional,
)

import numpy as np
from colabfold.batch import (
    get_queries,
    run,
)
from colabfold.utils import (
    setup_logging,
)

logger = logging.getLogger(__name__)

# Configuration constants
DEFAULT_NUM_MODELS = 5
MSA_DEPTH_MULTIPLIERS = (1, 2, 2, 2, 2, 2, 2)
INITIAL_MAX_MSA = 1
INITIAL_EXTRA_MSA = 2


class ColabFoldRunner(ABC):
    """Base class for orchestrating ColabFold predictions.

    Provides shared functionality for running predictions with ColabFold,
    including setup, query preparation, and result handling.
    """

    def __init__(
        self,
        search_dir: str,
        output_dir: str,
        pdb_name: str,
        num_seeds: int,
        model_type: str,
    ) -> None:
        """Initialize ColabFold runner.

        Args:
            search_dir: Directory containing MSA files or fasta sequences.
            output_dir: Directory where prediction results will be saved.
            pdb_name: Name of the target protein for file naming.
            num_seeds: Number of random seeds for predictions.
            model_type: ColabFold model type (e.g., 'ptm', 'monomer', 'multimer').
        """
        self.search_dir = search_dir
        self.output_dir = output_dir
        self.pdb_name = pdb_name
        self.num_seeds = num_seeds
        self.model_type = model_type
        self._setup_logging()

    def _setup_logging(self) -> None:
        """Configure logging for the prediction run."""
        output_path = Path(self.output_dir)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        setup_logging(output_path / "log.txt")

    def _run_colabfold(
        self,
        random_seed: int,
        max_seq: Optional[int] = None,
        max_extra_seq: Optional[int] = None,
    ) -> None:
        """Execute ColabFold prediction.

        Args:
            random_seed: Random seed for reproducibility.
            max_seq: Maximum sequence depth (optional).
            max_extra_seq: Maximum extra sequence depth (optional).

        Raises:
            RuntimeError: If query extraction or prediction fails.
        """
        try:
            queries, is_complex = get_queries(self.search_dir)
        except Exception as e:
            raise RuntimeError(f"Query extraction failed: {e}") from e

        run_kwargs = {
            "queries": queries,
            "result_dir": self.output_dir,
            "num_models": DEFAULT_NUM_MODELS,
            "is_complex": is_complex,
            "model_type": self.model_type,
            "num_seeds": int(self.num_seeds),
            "random_seed": int(random_seed),
            "data_dir": Path("."),
        }

        if max_seq is not None:
            run_kwargs["max_seq"] = int(max_seq)
        if max_extra_seq is not None:
            run_kwargs["max_extra_seq"] = int(max_extra_seq)

        try:
            run(**run_kwargs)
        except Exception as e:
            raise RuntimeError(f"Prediction failed: {e}") from e

    @staticmethod
    def _ensure_dir(dir_path: str) -> None:
        """Create directory if it doesn't exist."""
        os.makedirs(dir_path, exist_ok=True)

    @staticmethod
    def _move_results(source_pattern: str, dest_dir: str) -> None:
        """Move prediction results matching a glob pattern to destination directory.

        Args:
            source_pattern: Glob pattern for source directories/files.
            dest_dir: Destination directory path.

        Raises:
            RuntimeError: If no files matched or move operation fails.
        """
        ColabFoldRunner._ensure_dir(dest_dir)
        matched = glob.glob(source_pattern)
        if not matched:
            raise RuntimeError(f"No files matched pattern: {source_pattern}")
        for source_path in matched:
            shutil.move(source_path, dest_dir)
            logger.info("Moved %s -> %s", source_path, dest_dir)


class MSAMaxRunner(ColabFoldRunner):
    """Run ColabFold prediction with maximum (full) MSA depth."""

    def __init__(
        self,
        search_dir: str,
        output_dir: str,
        pdb_name: str,
        random_seed: int,
        num_seeds: int,
        model_type: str,
    ) -> None:
        """Initialize and execute full MSA prediction.

        Args:
            search_dir: Directory containing MSA files.
            output_dir: Local output directory for ColabFold results.
            pdb_name: Protein name for file naming.
            random_seed: Random seed for prediction (string or int).
            num_seeds: Number of random seeds.
            model_type: ColabFold model type.
        """
        # ColabFold writes locally; use a local output dir then move
        local_output_dir = f"{pdb_name}_predicted_models_full_rand_{random_seed}"
        super().__init__(search_dir, local_output_dir, pdb_name, num_seeds, model_type)

        logger.info("Running full MSA prediction for %s", pdb_name)
        self._run_colabfold(int(random_seed) if not isinstance(random_seed, int) else random_seed)

        # Move completed folder into predictions_all/<pdb_name>/ (matches original)
        dest_dir = str(Path("predictions_all") / pdb_name)
        self._move_results(local_output_dir + "/", dest_dir)


class MSAVariableRunner(ColabFoldRunner):
    """Run ColabFold predictions across varied MSA depths."""

    def __init__(
        self,
        search_dir: str,
        output_dir: str,
        pdb_name: str,
        random_seed,
        num_seeds: int,
        model_type: str,
        multipliers: tuple = MSA_DEPTH_MULTIPLIERS,
        initial_max_msa: int = INITIAL_MAX_MSA,
        initial_extra_msa: int = INITIAL_EXTRA_MSA,
    ) -> None:
        """Initialize and execute varied MSA predictions.

        Args:
            search_dir: Directory containing MSA files.
            output_dir: Base name for local output directories.
            pdb_name: Protein name for file naming.
            random_seed: Random seed for predictions (list, ndarray, or int).
            num_seeds: Number of random seeds.
            model_type: ColabFold model type.
            multipliers: Tuple of multipliers for MSA depth variation.
            initial_max_msa: Starting maximum sequence depth.
            initial_extra_msa: Starting extra sequence depth.
        """
        # Normalise random seed to string, matching original join logic
        if isinstance(random_seed, (list, np.ndarray)):
            random_seed_str = "".join(map(str, random_seed))
        else:
            random_seed_str = str(random_seed)

        # Base output_dir is just a prefix used to name local folders
        super().__init__(search_dir, output_dir, pdb_name, num_seeds, model_type)

        self.random_seed_str = random_seed_str
        self.multipliers = multipliers
        self.pdb_name = pdb_name

        logger.info(
            "Running variable MSA predictions for %s with %d depth variations",
            pdb_name,
            len(multipliers),
        )
        self._run_varied_predictions(initial_max_msa, initial_extra_msa)

    def _run_varied_predictions(self, initial_max_msa: int, initial_extra_msa: int) -> None:
        """Execute predictions across varied MSA depths and move results.

        Args:
            initial_max_msa: Starting maximum sequence depth.
            initial_extra_msa: Starting extra sequence depth.
        """
        max_msa = initial_max_msa
        extra_msa = initial_extra_msa
        dest_dir = str(Path("predictions_all") / self.pdb_name)

        for multiplier in self.multipliers:
            max_msa = int(max_msa * multiplier)
            extra_msa = int(extra_msa * multiplier)

            # Local folder name matches original pattern exactly
            local_output_dir = (
                f"{self.pdb_name}_predicted_models_rand_"
                f"{self.random_seed_str}_max_{max_msa}_ext_{extra_msa}"
            )

            logger.info(
                "Running prediction: max_seq=%d, max_extra_seq=%d -> %s",
                max_msa,
                extra_msa,
                local_output_dir,
            )

            # Write to local dir, then move (matches original mv logic)
            self.output_dir = local_output_dir
            self._setup_logging()

            self._run_colabfold(
                int(self.random_seed_str),
                max_seq=max_msa,
                max_extra_seq=extra_msa,
            )

            # Move completed folder into predictions_all/<pdb_name>/
            self._move_results(local_output_dir + "/", dest_dir)
