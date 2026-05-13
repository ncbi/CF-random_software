#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Utilities for converting multimer PDB files to single-chain structures.

This module provides functionality to convert multimer prediction outputs
into single-chain PDB files, removing TER records and extracting specific chains.
"""

import glob
import logging
from pathlib import (
    Path,
)

logger = logging.getLogger(__name__)


class ConvertM2S:
    """Convert multimer PDB structures to single-chain PDB files.

    Processes multimer prediction outputs by removing TER records and
    extracting individual chains for separate analysis.
    """

    def __init__(self, pred_path: str, pdb1_name: str, pdb2_name: str) -> None:
        """Initialize and execute multimer to single-chain conversion.

        Args:
            pred_path: Path to directory containing multimer predictions.
            pdb1_name: Name of first reference structure (used for naming).
            pdb2_name: Name of second reference structure (used for conversion).

        Raises:
            FileNotFoundError: If prediction directory or reference PDB not found.
            RuntimeError: If conversion commands fail.
        """
        self.pred_path = Path(pred_path)
        self.pdb1_name = pdb1_name
        self.pdb2_name = pdb2_name

        if not self.pred_path.exists():
            raise FileNotFoundError(f"Prediction directory not found: {pred_path}")

        try:
            self._remove_ter_records()
            self._extract_single_chains()
            logger.info("Successfully converted multimer predictions to single chains")
        except Exception as e:
            logger.error(f"Conversion failed: {e}")
            raise

    def _remove_ter_records(self) -> None:
        """Remove TER records from predicted multimer PDB files.

        Also creates cleaned versions of reference structures.
        """
        # Process predicted files
        files_list = glob.glob(str(self.pred_path / "0_unrelaxed*pdb"))

        for pred_file in files_list:
            try:
                output_file = pred_file.replace(".pdb", "").split("/")[-1]
                output_path = self.pred_path / f"rmTER_{output_file}.pdb"

                with open(pred_file, "r", encoding="utf-8") as infile:
                    with open(output_path, "w", encoding="utf-8") as outfile:
                        for line in infile:
                            if "TER" not in line:
                                outfile.write(line)

                logger.debug(f"Removed TER records: {output_path}")
            except Exception as e:
                logger.warning(f"Failed to process {pred_file}: {e}")
                continue

        # Process reference structure
        try:
            ref_file = Path(f"{self.pdb2_name}.pdb")
            if ref_file.exists():
                output_path = Path(f"{self.pdb2_name}_rmTER.pdb")
                with open(ref_file, "r", encoding="utf-8") as infile:
                    with open(output_path, "w", encoding="utf-8") as outfile:
                        for line in infile:
                            if "TER" not in line:
                                outfile.write(line)
                logger.debug(f"Removed TER records: {output_path}")
        except Exception as e:
            logger.warning(f"Failed to process reference {self.pdb2_name}: {e}")

    def _extract_single_chains(self) -> None:
        """Extract individual chains from multimer PDB files.

        Creates single-chain PDB files for the first chain found in each prediction.
        """
        files_list = glob.glob(str(self.pred_path / "0_unrelaxed*pdb"))

        for pred_file in files_list:
            try:
                output_basename = pred_file.replace(".pdb", "").split("/")[-1]
                output_path = self.pred_path / f"single_{output_basename}.pdb"

                # Extract first chain (up to first TER record)
                with open(pred_file, "r", encoding="utf-8") as infile:
                    with open(output_path, "w", encoding="utf-8") as outfile:
                        for line in infile:
                            outfile.write(line)
                            if "TER" in line:
                                break

                logger.debug(f"Extracted single chain: {output_path}")
            except Exception as e:
                logger.warning(f"Failed to extract chain from {pred_file}: {e}")
                continue
