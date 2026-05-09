#!/bin/env python3
# -*- coding: utf-8 -*-
"""Split protein chains into single-chain PDB files.

Simple CLI utility to extract individual chains from a multi-chain PDB.
"""

import argparse
import linecache

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


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb1", type=str, help="PDB structure for the target crystal structure")
    args = parser.parse_args()

    pdb1 = args.pdb1
    pdb1_name = pdb1.replace(".pdb", "")

    TER_count = 0
    with open(pdb1, "r") as file:
        for line in file:
            TER = line.split()
            TER_count += TER.count("TER")

    line_cnt = 0
    for i in range(0, TER_count):
        output_file_name = pdb1_name + "_" + chain_char[i] + ".pdb"

        if line_cnt == 0:
            with open(pdb1, "r") as infile, open(output_file_name, "w") as outfile:
                for line in infile:
                    outfile.write(line)
                    line_cnt = line_cnt + 1
                    if "TER " in line:
                        line_cnt = line_cnt + 1
                        break

        else:
            with open(pdb1, "r") as infile, open(output_file_name, "w") as outfile:
                for line in infile:
                    linecache.getline(pdb1, line_cnt)
                    outfile.write(linecache.getline(pdb1, line_cnt))
                    line_cnt = line_cnt + 1
                    if linecache.getline(pdb1, line_cnt) == "TER ":
                        line_cnt = line_cnt + 1
                        break
