#!/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 21 14:51:00 2024
 
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
import matplotlib.pyplot as plt
import glob
import random
import argparse



class CF_MSA_max():
    def __init__(self, search_dir, output_dir, pdb_name, rseed):

        #command = 'colabfold_batch --amber --use-gpu-relax --model-type alphafold2 --num-seeds 5 --random-seed ' + str(rseed) + search_dir + output_dir
        command = 'colabfold_batch --amber --use-gpu-rela --num-seeds 5 --random-seed ' + str(rseed) + search_dir + output_dir
        print(command)
        os.system(command)
        



class CF_MSA_var():
    def __init__(self, pdb1_name, search_dir, output_dir, rseed):
        #### shallow MSA section
        #### Global viarlable
        max_msa = 1
        ext_msa = 2
        random_seed = np.array(rseed) ## needed to remove future

        self.pdb1_name = pdb1_name



        for ran_seed in random_seed:
            max_msa = 1
            ext_msa = 2

            TMscores_random = []

            for multi in (1, 2, 2, 2, 2, 2, 2):
                max_msa = max_msa * multi
                ext_msa = ext_msa * multi

                #### Colabfold part
                #command = 'colabfold_batch --amber --use-gpu-relax --model-type alphafold2 --num-seed 5 --max-seq ' + str(max_msa) + ' --max-extra-seq ' + str(ext_msa) + search_dir + output_dir + str(ran_seed) + '_max_' + str(max_msa) + '_ext_' + str(ext_msa)
                command = 'colabfold_batch --amber --use-gpu-relax --num-seed 5 --max-seq ' + str(max_msa) + ' --max-extra-seq ' + str(ext_msa) + search_dir + output_dir + str(ran_seed) + '_max_' + str(max_msa) + '_ext_' + str(ext_msa)
                print(command)
                os.system(command)



            fin_pred_dir = pdb1_name + '_predicted_models_rand_' + str(ran_seed) + '_max_*'
            gen_dir = 'blind_prediction/' + pdb1_name

            if not os.path.exists(gen_dir):
                os.makedirs(gen_dir)
                mv_command = 'mv ' + fin_pred_dir + ' blind_prediction/' + pdb1_name
                print(mv_command); os.system(mv_command)
            else:
                mv_command = 'mv ' + fin_pred_dir + ' blind_prediction/' + pdb1_name
                print(mv_command); os.system(mv_command)
                




class prediction_all_blind():
    def __init__(self, pdb1_name, search_dir):
        
        pre_random_seed = np.random.randint(0, 16, 1)
        random_seed = ''.join(map(str, pre_random_seed))
        print(random_seed)
        output_dir = ' ' + pdb1_name + '_predicted_models_full_rand_' + str(random_seed)
        print(output_dir)


        ##### Perform predction with full-length MSA
        MSA_full = CF_MSA_max(search_dir, output_dir, pdb1_name, random_seed)
        pwd = os.getcwd() + '/'


        # Directory section
        gen_dir = 'blind_prediction/' + pdb1_name

        if not os.path.exists(gen_dir):
            os.mkdir(gen_dir)


        pred_dir = pdb1_name + '_predicted_models_full_rand_' + str(random_seed) + '/'
        mv_folder_cmd = 'mv ' + pred_dir + ' blind_prediction/' + pdb1_name
        print(mv_folder_cmd); os.system(mv_folder_cmd)



        ##### check out varied-MSA with (msa-max: 1, 2, 4, 8, 16, 32, 64) (msa-extra: 2, 4, 8, 16, 32, 64, 128)
        output_dir = ' ' + pdb1_name + '_predicted_models_rand_'
        random_seed = random.sample(range(100), 1)
        MSA_var = CF_MSA_var(pdb1_name, search_dir, output_dir, random_seed)



