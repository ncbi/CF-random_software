#!/bin/env python3
# -*- coding: utf-8 -*-
"""Utilities to compute TM-scores focused on fold-switching regions.

This module provides tools to extract region coordinates and compute
TM-align-based comparisons for predicted models.
"""

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

# Amino acid three-letter to one-letter code mapping
AA3TO1: Dict[str, str] = {
    "CYS": "C",
    "ASP": "D",
    "SER": "S",
    "GLN": "Q",
    "LYS": "K",
    "ILE": "I",
    "PRO": "P",
    "THR": "T",
    "PHE": "F",
    "ASN": "N",
    "GLY": "G",
    "HIS": "H",
    "LEU": "L",
    "ARG": "R",
    "TRP": "W",
    "ALA": "A",
    "VAL": "V",
    "GLU": "E",
    "TYR": "Y",
    "MET": "M",
}


class TMScoreFS:
    """Calculates TM-scores for fold-switching regions between PDB structures.

    Compares predicted protein models against original PDB structures,
    focusing on fold-switching regions. Computes TM-align scores for structural
    alignments.

    Attributes:
        tmscores_fs (numpy.ndarray): Array of TM-scores for fold-switching comparisons.
    """

    def get_coords(self, pdbfile: Union[str, Path], FSRange: str) -> Tuple[np.ndarray, str]:
        """Extracts CA coordinates and sequence for fold-switching region from PDB.

        Args:
            pdbfile (str or Path): Path to the PDB file.
            FSRange (str): Residue range for fold-switching region, e.g., "112-162".

        Returns:
            tuple: (coords_np, seq) where coords_np is numpy array of CA coordinates
                   (N x 3) and seq is the one-letter amino acid sequence string.
        """
        from Bio.PDB import (
            PDBParser,
        )

        pdb_parser = PDBParser(QUIET=True)
        struct = pdb_parser.get_structure("x", str(pdbfile))
        coords: List[List[float]] = []
        seq_dict: Dict[int, str] = {}

        # Convert string range to residue range
        start, stop = FSRange.split("-")
        res_range = range(int(start), int(stop) + 1)

        # Extract CA coordinates and sequence for residues in range
        for atom in struct.get_atoms():
            residue = atom.get_parent()
            res_id = residue.get_id()[1]
            resname = residue.get_resname()

            if res_id in res_range and atom.get_name() == "CA":
                x, y, z = atom.get_coord()
                coords.append([x, y, z])
                if res_id not in seq_dict:
                    seq_dict[res_id] = AA3TO1[resname]

        # Convert to numpy array and build sequence string
        coords_np = np.array(coords)
        sorted_data = sorted(seq_dict.items())
        seq = "".join(item[1] for item in sorted_data)

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
        import glob

        from tmtools import (
            tm_align,
        )

        tmscores: List[float] = []
        tmscores_ord: List[float] = []
        tmscores_rev: List[float] = []

        modelfiles = glob.glob(str(predfilepath) + "/*_unrelaxed*pdb")

        if len(modelfiles) == 0:
            return [0.0, 0.0, 0.0, 0.0, 0.0]

        for model in modelfiles:
            modelpath = Path(model)
            coords2, seq2 = self.get_coords(modelpath, res_range)

            # Calculate TM-score: coords1 vs coords2
            res = tm_align(coords1, coords2, seq1, seq2)
            tmscore_ord = round(res.tm_norm_chain1, 2)
            tmscores_ord.append(tmscore_ord)

            # Calculate reverse TM-score: coords2 vs coords1
            res = tm_align(coords2, coords1, seq2, seq1)
            tmscore_rev = round(res.tm_norm_chain1, 5)
            tmscores_rev.append(tmscore_rev)

        # Select the set with higher max score
        if np.max(tmscores_ord) > np.max(tmscores_rev):
            tmscores = tmscores_ord
        else:
            tmscores = tmscores_rev

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
        import glob

        # Get list of subdirectories
        all_sub_dir_paths = glob.glob(str(data_dir))
        tmscores_fs: List[List[float]] = []

        if len(all_sub_dir_paths) == 0:
            return

        # Compute coordinates and sequences outside loops for efficiency
        coords1, seq1 = self.get_coords(pdbfile1, res_range1)
        coords2, seq2 = self.get_coords(pdbfile2, res_range2)

        for subdir in all_sub_dir_paths:
            preddir = Path(subdir)
            if not preddir.exists():
                continue

            # Compare against both fold states
            tmscore_lst1 = self.get_tmscore(coords1, seq1, preddir, pred_range)
            tmscore_lst2 = self.get_tmscore(coords2, seq2, preddir, pred_range)

            tmscores_fs.extend([tmscore_lst1, tmscore_lst2])

        print("         ")
        tmscores_fs_array = np.array(tmscores_fs)
        print("tmscores_fs")
        self.tmscores_fs = tmscores_fs_array

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
        import os

        current_dir = os.getcwd() + "/"
        data_dir = Path(pred_path)

        # Read fold-switching ranges from file
        range_file = current_dir + "range_fs_pairs_all.txt"
        fs_res: Dict[str, Tuple[str, str]] = {}

        with open(range_file, "r") as infile:
            next(infile)  # Skip header
            for line in infile:
                line = line.strip()
                n1, n2, p1, p2, m1, m2 = line.split(",")
                if n1 not in fs_res:
                    fs_res[n1] = (p1, m1)
                if n2 not in fs_res:
                    fs_res[n2] = (p2, m2)

        print("Running for pair", pdb1_name, pdb2_name, end="..\n")
        print("comparing predictions of", pdb1_name, end="...\n")
