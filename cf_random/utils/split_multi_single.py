#!/bin/env python3
# -*- coding: utf-8 -*-
"""Split multimer PDBs into single-chain PDB files.

Utility to extract individual chains from multimer prediction files.
"""

import glob
import linecache


class SplitMultiToChains:
    """Splits multimer PDB files into single-chain PDB files."""

    def __init__(self, pred_path: str):

        chain_char = [
            "A",
            "B",
            "C",
            "D",
            "E",
            "F",
            "G",
            "H",
            "I",
            "J",
            "K",
            "L",
            "M",
            "N",
            "O",
            "P",
            "Q",
            "R",
            "S",
            "T",
            "U",
            "V",
            "W",
            "X",
            "Y",
            "Z",
        ]

        files_list = glob.glob(str(pred_path) + "/*_unrelaxed*pdb")

        for fl in files_list:
            ter_count = 0
            with open(fl, "r", encoding="utf-8") as file:
                for line in file:
                    ter = line.split()
                    ter_count += ter.count("TER")

            line_cnt = 0

            fl_name = fl.replace(".pdb", "")
            for i in range(0, ter_count):
                output_file_name = fl_name + "_chain_" + chain_char[i] + ".pdb"

                if line_cnt == 0:
                    with (
                        open(fl, "r", encoding="utf-8") as infile,
                        open(output_file_name, "w", encoding="utf-8") as outfile,
                    ):
                        for line in infile:
                            outfile.write(line)
                            line_cnt = line_cnt + 1
                            if "TER " in line:
                                line_cnt = line_cnt + 1
                                break

                else:
                    with (
                        open(fl, "r", encoding="utf-8") as infile,
                        open(output_file_name, "w", encoding="utf-8") as outfile,
                    ):
                        for line in infile:
                            linecache.getline(fl, line_cnt)
                            outfile.write(linecache.getline(fl, line_cnt))
                            line_cnt = line_cnt + 1
                            if linecache.getline(fl, line_cnt) == "TER ":
                                line_cnt = line_cnt + 1
                                break
