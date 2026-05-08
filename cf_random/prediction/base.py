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
            logger.error(f"Failed to extract queries from {self.search_dir}: {e}")
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
            logger.error(f"ColabFold prediction failed: {e}")
            raise RuntimeError(f"Prediction failed: {e}") from e

    @staticmethod
    def _ensure_dir(dir_path: str) -> None:
        """Create directory if it doesn't exist.

        Args:
            dir_path: Path to directory to create.
        """
        os.makedirs(dir_path, exist_ok=True)

    @staticmethod
    def _move_results(source_pattern: str, dest_dir: str) -> None:
        """Move prediction results to destination directory.

        Args:
            source_pattern: Glob pattern for source files.
            dest_dir: Destination directory path.

        Raises:
            RuntimeError: If move operation fails.
        """
        try:
            ColabFoldRunner._ensure_dir(dest_dir)
            moved = False
            for source_path in glob.glob(source_pattern):
                shutil.move(source_path, dest_dir)
                moved = True
            if not moved:
                raise RuntimeError(f"No files matched {source_pattern}")
            logger.info(f"Moved results to {dest_dir}")
        except Exception as e:
            logger.error(f"Failed to move results: {e}")
            raise RuntimeError(f"Move operation failed: {e}") from e


class MSAMaxRunner(ColabFoldRunner):
    """Run ColabFold prediction with maximum (full) MSA depth.

    This runner uses full-length MSAs for structure prediction,
    providing the most comprehensive predictions.
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
        """Initialize and execute full MSA prediction.

        Args:
            search_dir: Directory containing MSA files.
            output_dir: Output directory for results.
            pdb_name: Protein name for file naming.
            random_seed: Random seed for prediction.
            num_seeds: Number of random seeds.
            model_type: ColabFold model type.
        """
        super().__init__(search_dir, output_dir, pdb_name, num_seeds, model_type)
        logger.info(f"Running full MSA prediction for {pdb_name}")
        self._run_colabfold(random_seed)


class MSAVariableRunner(ColabFoldRunner):
    """Run ColabFold predictions across varied MSA depths.

    This runner systematically reduces MSA depth by specified multipliers,
    allowing analysis of prediction robustness across different MSA sizes.
    """

    def __init__(
        self,
        search_dir: str,
        output_dir: str,
        pdb_name: str,
        random_seed: int,
        num_seeds: int,
        model_type: str,
        multipliers: tuple = MSA_DEPTH_MULTIPLIERS,
        initial_max_msa: int = INITIAL_MAX_MSA,
        initial_extra_msa: int = INITIAL_EXTRA_MSA,
    ) -> None:
        """Initialize and execute varied MSA predictions.

        Args:
            search_dir: Directory containing MSA files.
            output_dir: Base output directory for results.
            pdb_name: Protein name for file naming.
            random_seed: Random seed for predictions.
            num_seeds: Number of random seeds.
            model_type: ColabFold model type.
            multipliers: Tuple of multipliers for MSA depth variation.
            initial_max_msa: Starting maximum sequence depth.
            initial_extra_msa: Starting extra sequence depth.
        """
        super().__init__(search_dir, output_dir, pdb_name, num_seeds, model_type)

        # Convert random seed to string if needed
        if isinstance(random_seed, (list, np.ndarray)):
            random_seed_str = "".join(map(str, random_seed))
        else:
            random_seed_str = str(random_seed)

        self.random_seed_str = random_seed_str
        self.multipliers = multipliers
        logger.info(
            f"Running variable MSA prediction for {pdb_name} "
            f"with {len(multipliers)} depth variations"
        )
        self._run_varied_predictions(initial_max_msa, initial_extra_msa)

    def _run_varied_predictions(self, initial_max_msa: int, initial_extra_msa: int) -> None:
        """Execute predictions across varied MSA depths.

        Args:
            initial_max_msa: Starting maximum sequence depth.
            initial_extra_msa: Starting extra sequence depth.
        """
        max_msa = initial_max_msa
        extra_msa = initial_extra_msa

        for multiplier in self.multipliers:
            max_msa = int(max_msa * multiplier)
            extra_msa = int(extra_msa * multiplier)

            output_dir_var = (
                f"{self.output_dir}{self.random_seed_str}_max_{max_msa}_ext_{extra_msa}"
            )

            logger.info(f"Running prediction: max_seq={max_msa}, max_extra_seq={extra_msa}")

            # Override output directory temporarily
            original_output_dir = self.output_dir
            self.output_dir = output_dir_var
            self._setup_logging()

            try:
                self._run_colabfold(
                    int(self.random_seed_str),
                    max_seq=max_msa,
                    max_extra_seq=extra_msa,
                )
            finally:
                # Restore original output directory
                self.output_dir = original_output_dir
