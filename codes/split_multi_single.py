#!/bin/env python3
# -*- coding: utf-8 -*-
"""

Converting the multimer PDB to a single PDB file

Created on Tue Dec 24 14:51:00 2025
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
import random
import linecache
import argparse



class split_multi_to_chains():
    def __init__(self, pred_path):


        chain_char = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P',
                'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z']


        current_dir = os.getcwd() + '/' 
        data_dir = Path(pred_path) # Path to the predicted models

        files_list = (glob.glob(str(pred_path) + "/*_unrelaxed*pdb"))

        for fl in files_list:
            TER_count = 0
            with open(fl, 'r') as file:
                for line in file: 
                    TER = line.split()
                    TER_count += TER.count("TER")

            line_cnt = 0

            fl_name = fl.replace('.pdb','')
            for i in range(0, TER_count):
                output_file_name = fl_name + '_chain_' + chain_char[i] + '.pdb' 

                if line_cnt == 0:
                    with open(fl, 'r') as infile, open(output_file_name, 'w') as outfile:
                        for line in infile:
                            outfile.write(line)
                            line_cnt = line_cnt + 1
                            if "TER " in line:
                                line_cnt = line_cnt + 1
                                break

                else:
                    with open(fl, 'r') as infile, open(output_file_name, 'w') as outfile:
                        for line in infile:
                            linecache.getline(fl, line_cnt)
                            outfile.write(linecache.getline(fl, line_cnt))
                            line_cnt = line_cnt + 1
                            if linecache.getline(fl, line_cnt) == "TER ":
                                line_cnt = line_cnt + 1
                                break






