#!/bin/env python3
# -*- coding: utf-8 -*-
"""Multimer-aware TM-score calculations for fold-switching regions.

Provides utilities to compute TM-scores for multimeric predicted models
and compare fold-switching regions to reference structures.
"""

from cf_random.utils.constants import (
    AA3TO1,
)


class TMScoreFSMulti:
    """Calculates TM-scores for fold-switching regions in multimeric protein structures.

    This class compares predicted protein models against original PDB structures,
    focusing on fold-switching regions. It computes TM-align scores for structural
    alignments between predicted and reference structures.

    Attributes:
        tmscores_fs (numpy.ndarray): Array of TM-scores for fold-switching comparisons.
    """

    def __init__(self, pred_path, pdb1, pdb1_name, pdb2, pdb2_name):
        """Initializes TM-score calculation for fold-switching multimer analysis.

        Args:
            pred_path (str): Path to directory containing predicted model subdirectories.
            pdb1 (str or Path): Path to first PDB file.
            pdb1_name (str): Name/ID of first PDB structure.
            pdb2 (str or Path): Path to second PDB file.
            pdb2_name (str): Name/ID of second PDB structure.

        Raises:
            SystemExit: If PDB names are not found in range file.
        """
        import os
        import sys
        from pathlib import (
            Path,
        )

        current_dir = os.getcwd() + "/"
        data_dir = Path(pred_path)

        # the range of the fold-switching region
        range_file = current_dir + "range_fs_pairs_all.txt"

        # convert this file into a dictionary for reference later
        fs_res = {}

        # The range_file file has the fold-switching residue ranges
        # for the original PDB/PDB1, alternate PDB/PDB2
        # Predicted model for PDB1, predicted model for PDB2
        with open(range_file, "r") as Infile:
            next(Infile)  # skip header line "# pdb1,pdb2,pred1,pred2"
            for line in Infile:
                line = line.strip()
                n1, n2, p1, p2, m1, m2 = line.split(",")
                # the value of the dictionary is a tuple
                # the first element of tuple is the fs range in the original PDB
                # followed by the range in the predicted model
                if n1 not in fs_res:
                    fs_res[n1] = (p1, m1)
                if n2 not in fs_res:
                    fs_res[n2] = (p2, m2)

        print("Running for pair ", pdb1_name, pdb2_name, end="..")
        print("                ")
        print("comparing predictions of ", pdb1_name, end="...")
        print("                ")

        try:
            range_pdb1 = fs_res[
                pdb1_name
            ]  # so if pdb1 is '1nqd_A', fs_res['1nqd_A']=('895-919', '1-33')
            range_pdb2 = fs_res[
                pdb2_name
            ]  # and if pdb2 is '1nqj_B', fs_res['1nqj_B']=('894-919', '1-33')
        except:
            print("check PDBIDs ", pdb1_name, pdb2_name)
            sys.exit(1)

        range_pred = range_pdb1[1]
        self.run_for_models(pdb1, pdb2, data_dir, range_pred, range_pdb1[0], range_pdb2[0])

    def get_coords(self, pdbfile, FSRange):
        """Extracts coordinates and sequence for fold-switching region from PDB file.

        Args:
            pdbfile (str or Path): Path to the PDB file.
            FSRange (str): Residue range for fold-switching region, e.g., "112-162".

        Returns:
            tuple: (coords_np, seq) where coords_np is numpy array of CA coordinates
                   and seq is the one-letter amino acid sequence.
        """
        import numpy as np
        from Bio.PDB import (
            PDBParser,
        )

        pdb_parser = PDBParser(QUIET=True)
        seq = ""
        struct = pdb_parser.get_structure("x", str(pdbfile))
        coords = []
        seq_dict = {}

        # for residues within a certain range, using numpy to save the coords
        # and save the sequence as a dict and then sorted list of tuples
        # return the coords and the seq

        # convert str to residue range for the fs region
        start, stop = FSRange.split("-")
        res_range = range(int(start), int(stop) + 1)

        for atom in struct.get_atoms():
            residue = atom.get_parent()  # from atom we can get the parent residue
            res_id = residue.get_id()[1]
            resname = residue.get_resname()
            if res_id in res_range and atom.get_name() == "CA":
                x, y, z = atom.get_coord()
                coords.append([x, y, z])
                if res_id not in seq_dict:
                    seq_dict[res_id] = AA3TO1[resname]

        # convert to np array
        coords_np = np.array(coords)
        # sort the seq_dict by keys a.k.a res_ids
        sorted_data = sorted(seq_dict.items())
        for i in sorted_data:
            seq += i[1]

        return coords_np, seq

    def get_tmscore(self, coords1, seq1, predfilepath, res_range):
        """Calculates TM-scores between reference structure and predicted models.

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
        from pathlib import (
            Path,
        )

        from tmtools import (
            tm_align,
        )

        tmscores = []
        # modelfiles = sorted(glob.glob(str(predfilepath) + "/*_unrelaxed*pdb"))
        modelfiles = glob.glob(str(predfilepath) + "/single*_unrelaxed*pdb")

        if len(modelfiles) == 0:
            tmscores = [0.0, 0.0, 0.0, 0.0, 0.0]
            return tmscores

        for model in modelfiles:
            modelpath = Path(model)
            coords2, seq2 = self.get_coords(modelpath, res_range)
            res = tm_align(coords1, coords2, seq1, seq2)
            tmscore = round(res.tm_norm_chain1, 2)  # wrt to model
            tmscores.append(tmscore)

        return tmscores

    def run_for_models(self, pdbfile1, pdbfile2, data_dir, pred_range, res_range1, res_range2):
        """Compares predicted models against both original PDB structures.

        Calculates TM-scores for fold-switching regions by comparing predicted models
        against both fold states (pdbfile1 and pdbfile2).

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
        from pathlib import (
            Path,
        )

        import numpy as np

        # print(res_range1,res_range2)

        # get list of subdirectories
        all_sub_dir_paths = glob.glob(str(data_dir))
        tmscores_fs = []

        if len(all_sub_dir_paths) == 0:
            return

        # Compute coordinates and sequences outside the loop for efficiency
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
        tmscores_fs = np.array(tmscores_fs)
        print("tmscores_fs")
        self.tmscores_fs = tmscores_fs
