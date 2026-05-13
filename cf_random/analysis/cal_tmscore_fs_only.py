#!/bin/env python3
# -*- coding: utf-8 -*-
"""Utilities to compute TM-scores focused on fold-switching regions.

This module provides tools to extract region coordinates and compute
TM-align-based comparisons for predicted models.
"""

import glob
import logging
import os
import sys
from pathlib import (
    Path,
)
from typing import (
    Dict,
    List,
    Tuple,
    Union,
)

import numpy as np
from Bio.PDB import (
    PDBParser,
)
from tmtools import (
    tm_align,
)

from cf_random.utils.constants import (
    AA3TO1,
)

logger = logging.getLogger(__name__)


class TMScoreFS:
    """Calculates TM-scores for fold-switching regions between PDB structures.

    Compares predicted protein models against original PDB structures,
    focusing on fold-switching regions. Computes TM-align scores for structural
    alignments.

    Attributes:
        tmscores_fs (numpy.ndarray): Array of TM-scores for fold-switching comparisons.
    """

    def __init__(
        self,
        pred_path: Union[str, Path],
        pdb1: Union[str, Path],
        pdb1_name: str,
        pdb2: Union[str, Path],
        pdb2_name: str,
    ) -> None:
        """Initializes TM-score calculation for fold-switching analysis.

        Args:
            pred_path (str or Path): Path to directory containing predicted model subdirectories.
            pdb1 (str or Path): Path to first PDB file.
            pdb1_name (str): Name/ID of first PDB structure.
            pdb2 (str or Path): Path to second PDB file.
            pdb2_name (str): Name/ID of second PDB structure.

        Raises:
            SystemExit: If PDB names are not found in range file.
        """
        current_dir = os.getcwd() + "/"

        range_file = current_dir + "range_fs_pairs_all.txt"
        fs_res: Dict[str, Tuple[str, str]] = {}

        with open(range_file, "r", encoding="utf-8") as infile:
            next(infile)  # Skip header
            for line in infile:
                line = line.strip()
                n1, n2, p1, p2, m1, m2 = line.split(",")
                if n1 not in fs_res:
                    fs_res[n1] = (p1, m1)
                if n2 not in fs_res:
                    fs_res[n2] = (p2, m2)

        logger.info("Running fold-switching TM-score for pair %s / %s", pdb1_name, pdb2_name)

        try:
            range_pdb1 = fs_res[pdb1_name]
            range_pdb2 = fs_res[pdb2_name]
        except KeyError:
            logger.error(
                "PDB ID(s) not found in range file — check identifiers: %s, %s",
                pdb1_name,
                pdb2_name,
            )
            sys.exit(1)

        range_pred = range_pdb1[1]
        self.run_for_models(pdb1, pdb2, pred_path, range_pred, range_pdb1[0], range_pdb2[0])

    def get_coords(self, pdbfile: Union[str, Path], fs_range: str) -> Tuple[np.ndarray, str]:
        """Extracts CA coordinates and sequence for fold-switching region from PDB.

        Args:
            pdbfile (str or Path): Path to the PDB file.
            fs_range (str): Residue range for fold-switching region, e.g., "112-162".

        Returns:
            tuple: (coords_np, seq) where coords_np is numpy array of CA coordinates
                   (N x 3) and seq is the one-letter amino acid sequence string.
        """
        pdb_parser = PDBParser(QUIET=True)
        struct = pdb_parser.get_structure("x", str(pdbfile))
        coords: List[List[float]] = []
        seq_dict: Dict[int, str] = {}

        start, stop = fs_range.split("-")
        res_range = range(int(start), int(stop) + 1)

        for atom in struct.get_atoms():
            residue = atom.get_parent()
            res_id = residue.get_id()[1]
            resname = residue.get_resname()

            if res_id in res_range and atom.get_name() == "CA":
                x, y, z = atom.get_coord()
                coords.append([x, y, z])
                if res_id not in seq_dict:
                    seq_dict[res_id] = AA3TO1[resname]

        coords_np = np.array(coords)
        seq = "".join(item[1] for item in sorted(seq_dict.items()))

        logger.debug("Extracted %d CA atoms from %s (range %s)", len(coords_np), pdbfile, fs_range)
        return coords_np, seq

    def get_tmscore(
        self, coords1: np.ndarray, seq1: str, predfilepath: Union[str, Path], res_range: str
    ) -> List[float]:
        """Calculates TM-scores between reference and predicted structures.

        Args:
            coords1 (numpy.ndarray): Coordinates of reference structure.
            seq1 (str): Sequence of reference structure.
            predfilepath (str or Path): Path to directory containing predicted models.
            res_range (str): Residue range for fold-switching in predicted models.

        Returns:
            list: TM-scores for each predicted model (rounded to 2 decimals).
                  Returns [0.0, 0.0, 0.0, 0.0, 0.0] if no models found.
        """
        modelfiles = glob.glob(str(predfilepath) + "/*_unrelaxed*pdb")

        if not modelfiles:
            logger.warning("No unrelaxed model files found in %s", predfilepath)
            return [0.0, 0.0, 0.0, 0.0, 0.0]

        tmscores_ord: List[float] = []
        tmscores_rev: List[float] = []

        for model in modelfiles:
            coords2, seq2 = self.get_coords(Path(model), res_range)

            res = tm_align(coords1, coords2, seq1, seq2)
            tmscores_ord.append(round(res.tm_norm_chain1, 2))

            # Note: both directions rounded to 2 for consistent comparison
            res = tm_align(coords2, coords1, seq2, seq1)
            tmscores_rev.append(round(res.tm_norm_chain1, 2))

        if np.max(tmscores_ord) > np.max(tmscores_rev):
            tmscores = tmscores_ord
        else:
            tmscores = tmscores_rev

        logger.debug("TM-scores for %s: %s", predfilepath, tmscores)
        return tmscores

    def run_for_models(
        self,
        pdbfile1: Union[str, Path],
        pdbfile2: Union[str, Path],
        data_dir: Union[str, Path],
        pred_range: str,
        res_range1: str,
        res_range2: str,
    ) -> None:
        """Compares predicted models against both original PDB structures.

        Calculates TM-scores for fold-switching regions by comparing predicted
        models against both fold states (pdbfile1 and pdbfile2).

        Args:
            pdbfile1 (str or Path): Path to first PDB structure (Fold1).
            pdbfile2 (str or Path): Path to second PDB structure (Fold2).
            data_dir (str or Path): Path to directory containing predicted model subdirectories.
            pred_range (str): Residue range for fold-switching in predicted models.
            res_range1 (str): Residue range for fold-switching in pdbfile1.
            res_range2 (str): Residue range for fold-switching in pdbfile2.

        Returns:
            None: Stores results in self.tmscores_fs attribute.
        """
        all_sub_dir_paths = glob.glob(str(data_dir))

        if not all_sub_dir_paths:
            logger.warning("No prediction subdirectories matched pattern: %s", data_dir)
            return

        logger.info(
            "Found %d prediction director%s under %s",
            len(all_sub_dir_paths),
            "y" if len(all_sub_dir_paths) == 1 else "ies",
            data_dir,
        )

        coords1, seq1 = self.get_coords(pdbfile1, res_range1)
        coords2, seq2 = self.get_coords(pdbfile2, res_range2)

        tmscores_fs: List[List[float]] = []

        for subdir in all_sub_dir_paths:
            preddir = Path(subdir)
            if not preddir.exists():
                logger.debug("Skipping missing directory: %s", preddir)
                continue

            tmscores_fs.append(self.get_tmscore(coords1, seq1, preddir, pred_range))
            tmscores_fs.append(self.get_tmscore(coords2, seq2, preddir, pred_range))

        self.tmscores_fs = np.array(tmscores_fs)
        logger.info("TM-score array shape: %s", self.tmscores_fs.shape)
