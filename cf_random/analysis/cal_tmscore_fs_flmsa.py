#!/bin/env python3
# -*- coding: utf-8 -*-
"""Fold-switching TM-score computation.

Extracts coordinates for FS regions and computes TM-align scores
against both reference structures.
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
    Optional,
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

    Compares predicted protein models against reference PDB structures,
    focusing on fold-switching regions. Supports both monomer and multimer
    predicted models.

    Attributes:
        tmscores_fs (numpy.ndarray): Array of TM-scores for fold-switching comparisons.
    """

    def __init__(
        self,
        pred_path: str,
        pdb1: Union[str, Path],
        pdb1_name: str,
        pdb2: Union[str, Path],
        pdb2_name: str,
        model_glob: str = "*_unrelaxed*pdb",
    ) -> None:
        """Initializes TM-score calculation for fold-switching analysis.

        Args:
            pred_path (str): Path or glob pattern to directory containing predicted models.
            pdb1 (str or Path): Path to first PDB file.
            pdb1_name (str): Name/ID of first PDB structure.
            pdb2 (str or Path): Path to second PDB file.
            pdb2_name (str): Name/ID of second PDB structure.
            model_glob (str): Glob pattern for model files within each prediction directory.
                Defaults to "*_unrelaxed*pdb". Use "single_0_unrelaxed*pdb" for multimer.

        Raises:
            SystemExit: If PDB names are not found in range file.
        """
        self.tmscores_fs = None
        self.model_glob = model_glob

        range_file = os.path.join(os.getcwd(), "range_fs_pairs_all.txt")
        fs_res: Dict[str, Tuple[str, str]] = {}

        with open(range_file, "r", encoding="utf-8") as f:
            next(f)  # Skip header
            for line in f:
                n1, n2, p1, p2, m1, m2 = line.strip().split(",")
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

        self._run_for_models(pdb1, pdb2, pred_path, range_pdb1[1], range_pdb1[0], range_pdb2[0])

    def _get_coords(self, pdb_file: Union[str, Path], fs_range: str) -> Tuple[np.ndarray, str]:
        """Extracts CA coordinates and sequence for a fold-switching region from a PDB file.

        Args:
            pdb_file (str or Path): Path to the PDB file.
            fs_range (str): Residue range for fold-switching region, e.g., "112-162".

        Returns:
            tuple: (coords_np, seq) where coords_np is a numpy array of CA coordinates
                   (N x 3) and seq is the one-letter amino acid sequence string.
        """
        struct = PDBParser(QUIET=True).get_structure("x", str(pdb_file))
        coords: List[List[float]] = []
        seq_dict: Dict[int, str] = {}

        start, stop = fs_range.split("-")
        res_range = range(int(start), int(stop) + 1)

        for atom in struct.get_atoms():
            residue = atom.get_parent()
            res_id = residue.get_id()[1]
            if res_id in res_range and atom.get_name() == "CA":
                x, y, z = atom.get_coord()
                coords.append([x, y, z])
                if res_id not in seq_dict:
                    seq_dict[res_id] = AA3TO1[residue.get_resname()]

        coords_np = np.array(coords)
        seq = "".join(v for _, v in sorted(seq_dict.items()))

        logger.debug("Extracted %d CA atoms from %s (range %s)", len(coords_np), pdb_file, fs_range)
        return coords_np, seq

    def _get_tmscore(
        self,
        coords1: np.ndarray,
        seq1: str,
        pred_dir: Union[str, Path],
        res_range: str,
    ) -> List[float]:
        """Calculates TM-scores between a reference structure and predicted models.

        Scores are computed in both directions and the higher-scoring direction is returned.

        Args:
            coords1 (numpy.ndarray): CA coordinates of the reference structure.
            seq1 (str): Sequence of the reference structure.
            pred_dir (str or Path): Path to directory containing predicted models.
            res_range (str): Residue range for the fold-switching region in predicted models.

        Returns:
            list: TM-scores for each predicted model (rounded to 2 decimals).
                  Returns [0.0, 0.0, 0.0, 0.0, 0.0] if no models are found.
        """
        model_files = sorted(glob.glob(str(pred_dir) + f"/{self.model_glob}"))

        if not model_files:
            logger.warning("No unrelaxed model files found in %s", pred_dir)
            return [0.0, 0.0, 0.0, 0.0, 0.0]

        tmscores_fwd: List[float] = []
        tmscores_rev: List[float] = []

        for model in model_files:
            coords2, seq2 = self._get_coords(model, res_range)
            tmscores_fwd.append(round(tm_align(coords1, coords2, seq1, seq2).tm_norm_chain1, 2))
            tmscores_rev.append(round(tm_align(coords2, coords1, seq2, seq1).tm_norm_chain1, 2))

        tmscores = tmscores_fwd if np.max(tmscores_fwd) > np.max(tmscores_rev) else tmscores_rev
        logger.debug("TM-scores for %s: %s", pred_dir, tmscores)
        return tmscores

    def _run_for_models(
        self,
        pdb_file1: Union[str, Path],
        pdb_file2: Union[str, Path],
        pred_path: Union[str, Path],
        pred_range: str,
        res_range1: str,
        res_range2: str,
    ) -> None:
        """Compares predicted models against both reference PDB structures.

        Args:
            pdb_file1 (str or Path): Path to the first reference PDB (Fold1).
            pdb_file2 (str or Path): Path to the second reference PDB (Fold2).
            pred_path (str or Path): Path or glob pattern to predicted model directories.
            pred_range (str): Residue range for fold-switching in predicted models.
            res_range1 (str): Residue range for fold-switching in pdb_file1.
            res_range2 (str): Residue range for fold-switching in pdb_file2.
        """
        all_subdirs = glob.glob(str(pred_path))

        if not all_subdirs:
            logger.warning("No prediction subdirectories matched pattern: %s", pred_path)
            return

        logger.info(
            "Found %d prediction director%s under %s",
            len(all_subdirs),
            "y" if len(all_subdirs) == 1 else "ies",
            pred_path,
        )

        coords1, seq1 = self._get_coords(pdb_file1, res_range1)
        coords2, seq2 = self._get_coords(pdb_file2, res_range2)

        tmscores_fs: List[List[float]] = []

        for subdir in all_subdirs:
            preddir = Path(subdir)
            if not preddir.exists():
                logger.debug("Skipping missing directory: %s", preddir)
                continue
            tmscores_fs.append(self._get_tmscore(coords1, seq1, preddir, pred_range))

        for subdir in all_subdirs:
            preddir = Path(subdir)
            if not preddir.exists():
                logger.debug("Skipping missing directory: %s", preddir)
                continue
            tmscores_fs.append(self._get_tmscore(coords2, seq2, preddir, pred_range))

        self.tmscores_fs = np.array(tmscores_fs)
        logger.info("TM-score array shape: %s", self.tmscores_fs.shape)
