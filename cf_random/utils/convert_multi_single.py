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

    def _find_unrelaxed_files(self) -> list:
        """Find unrelaxed PDB files in the prediction directory.

        Tries the ColabFold default prefix first, then falls back to a
        wildcard match to handle sequence-ID-prefixed filenames.

        Returns:
            List of matched file path strings.
        """
        files = glob.glob(str(self.pred_path / "0_unrelaxed*pdb"))
        if not files:
            files = glob.glob(str(self.pred_path / "*_unrelaxed*pdb"))
            if files:
                logger.debug(
                    "Default prefix not found; matched %d file(s) with wildcard in %s",
                    len(files),
                    self.pred_path,
                )
        return files

    def _remove_ter_records(self) -> None:
        """Remove TER records from predicted multimer PDB files.

        Also creates cleaned versions of reference structures.
        """
        for pred_file in self._find_unrelaxed_files():
            try:
                output_file = pred_file.replace(".pdb", "").split("/")[-1]
                output_path = self.pred_path / f"rmTER_{output_file}.pdb"

                with open(pred_file, "r", encoding="utf-8") as infile:
                    with open(output_path, "w", encoding="utf-8") as outfile:
                        for line in infile:
                            if "TER" not in line:
                                outfile.write(line)

                logger.debug("Removed TER records: %s", output_path)
            except Exception as e:
                logger.warning("Failed to process %s: %s", pred_file, e)
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
                logger.debug("Removed TER records: %s", output_path)
        except Exception as e:
            logger.warning("Failed to process reference %s: %s", self.pdb2_name, e)

    def _extract_single_chains(self) -> None:
        """Extract individual chains from multimer PDB files.

        Creates single-chain PDB files for the first chain found in each prediction.
        """
        for pred_file in self._find_unrelaxed_files():
            try:
                output_basename = pred_file.replace(".pdb", "").split("/")[-1]
                output_path = self.pred_path / f"single_{output_basename}.pdb"

                with open(pred_file, "r", encoding="utf-8") as infile:
                    with open(output_path, "w", encoding="utf-8") as outfile:
                        for line in infile:
                            outfile.write(line)
                            if "TER" in line:
                                break

                logger.debug("Extracted single chain: %s", output_path)
            except Exception as e:
                logger.warning("Failed to extract chain from %s: %s", pred_file, e)
                continue
