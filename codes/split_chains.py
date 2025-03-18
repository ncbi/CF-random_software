#!/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 17 09:29:00 2024

Splitting protein chain as each single file

@author: Myeongsang (Samuel) Lee
"""
import re
import Bio
import os
from os import listdir
from os.path import isfile, join
import sys
from pathlib import Path
import numpy as np
from numpy import genfromtxt
import glob
import linecache
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("--pdb1", type=str, help='PDB structure for the target crystal structure')
args = parser.parse_args() 


chain_char = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P',
        'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z']

pdb1 = args.pdb1; pdb1_name = pdb1.replace('.pdb','')


TER_count = 0
with open(pdb1, 'r') as file:
    for line in file:
        TER = line.split()
        TER_count += TER.count("TER")





line_cnt = 0
for i in range(0, TER_count):
    output_file_name = pdb1_name + '_' + chain_char[i] + '.pdb'

    if line_cnt == 0:
        with open(pdb1, 'r') as infile, open(output_file_name, 'w') as outfile:
            for line in infile:
                outfile.write(line)
                line_cnt = line_cnt + 1
                if "TER " in line:
                    line_cnt = line_cnt + 1
                    break

    else:
        with open(pdb1, 'r') as infile, open(output_file_name, 'w') as outfile:
            for line in infile:
                linecache.getline(pdb1, line_cnt)
                outfile.write(linecache.getline(pdb1, line_cnt))
                line_cnt = line_cnt + 1
                if linecache.getline(pdb1, line_cnt) == "TER ":
                    line_cnt = line_cnt + 1
                    break



        

