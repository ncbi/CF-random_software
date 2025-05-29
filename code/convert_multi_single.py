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
import argparse





class convert_m2s():
    def __init__(self, pred_path, pdb1_name, pdb2_name):
        current_dir = os.getcwd() + '/' 
        data_dir = Path(pred_path) # Path to the predicted models

        files_list = (glob.glob(str(pred_path) + "/*_unrelaxed*pdb"))
        print(files_list)
        for fl in files_list:
            fl_name = fl.replace('.pdb','')
            predicted_name = fl_name.split('/')[1]
            #convert = "awk '!/TER/' " + fl + " > " + fl_name + "_converted.pdb"
            convert = "awk '!/TER/' " + fl + " > " + fl_name.split('/')[0] + '/' + "rmTER_" + predicted_name + ".pdb"
            print(convert)
            os.system(convert)


        convert_pdb2 = "awk '!/TER/' " + pdb2_name + ".pdb > " + pdb2_name + "_rmTER.pdb"
        print(convert_pdb2); os.system(convert_pdb2)


        ##### extract a single chain from multimer
        TER_count = 0

        for fl in files_list:
            fl_name = fl.replace('.pdb','')
            predicted_name = fl_name.split('/')[1]

            with open(fl, 'r') as file:
                for line in file:
                    TER = line.split()
                    TER_count += TER.count("TER")


            line_cnt = 0
            for i in range(0, 2):
                output_file_name = fl_name.split('/')[0] + '/' + "single_" + predicted_name + ".pdb"

                if line_cnt == 0:
                    with open(fl, 'r') as infile, open(output_file_name, 'w') as outfile:
                        for line in infile:
                            outfile.write(line)
                            line_cnt = line_cnt + 1
                            if "TER " in line:
                                line_cnt = line_cnt + 1
                                break

        #line_cnt = 0
        ##for i in range(0, TER_count):
        #for i in range(0, 2):
        #    output_file_name = pdb2_name[0:4] + '_multi.pdb'

        #    if line_cnt == 0:
        #        with open(pdb2, 'r') as infile, open(output_file_name, 'w') as outfile:
        #            for line in infile:
        #                outfile.write(line)
        #                line_cnt = line_cnt + 1
        #                if "TER " in line:
        #                    line_cnt = line_cnt + 1
        #                    break

        #pdb2_name_multi = output_file_name.replace('.pdb','')


