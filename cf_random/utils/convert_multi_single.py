#!/bin/env python3
# -*- coding: utf-8 -*-
"""Utilities for converting multimer PDBs to single-chain PDBs.

This module provides a small utility class used by the TM-score
calculation routines to generate single-chain PDBs from multimer
predictions. The implementation intentionally performs filesystem
operations in the constructor for compatibility with the original
behaviour.

Classes:
    ConvertM2S: Convert multimer prediction outputs into single-chain
        PDB files.
"""

import glob
import os


class ConvertM2S:
    """Convert multimer prediction outputs to single-chain PDB files.

    The original code executed the conversion in the class
    constructor. To remain backward compatible we preserve that
    behaviour.

    Args:
        pred_path: Path to prediction output directory (globbed for *_unrelaxed*pdb).
        pdb1_name: Reference PDB name for conversion of the second structure.
        pdb2_name: Second PDB name used when producing *_rmTER.pdb.
    """

    def __init__(self, pred_path: str, pdb1_name: str, pdb2_name: str) -> None:
        files_list = glob.glob(str(pred_path) + "/*_unrelaxed*pdb")
        print(files_list)
        for fl in files_list:
            fl_name = fl.replace(".pdb", "")
            predicted_name = fl_name.split("/")[1]
            convert_cmd = (
                "awk '!/TER/' "
                + fl
                + " > "
                + fl_name.split("/")[0]
                + "/"
                + "rmTER_"
                + predicted_name
                + ".pdb"
            )
            print(convert_cmd)
            os.system(convert_cmd)

        convert_pdb2 = "awk '!/TER/' " + pdb2_name + ".pdb > " + pdb2_name + "_rmTER.pdb"
        print(convert_pdb2)
        os.system(convert_pdb2)

        # Extract a single chain from each multimer prediction
        for fl in files_list:
            fl_name = fl.replace(".pdb", "")
            predicted_name = fl_name.split("/")[1]

            line_cnt = 0
            for _ in range(0, 2):
                output_file_name = fl_name.split("/")[0] + "/" + "single_" + predicted_name + ".pdb"

                if line_cnt == 0:
                    with open(fl, "r") as infile, open(output_file_name, "w") as outfile:
                        for line in infile:
                            outfile.write(line)
                            line_cnt = line_cnt + 1
                            if "TER " in line:
                                line_cnt = line_cnt + 1
                                break


# Backwards-compatible alias used throughout the codebase
convert_m2s = ConvertM2S
