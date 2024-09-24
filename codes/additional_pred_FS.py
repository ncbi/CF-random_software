#!/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 22 13:40:00 2024
 
@author: Myeongsang (Samuel) Lee
"""
import os
import sys
from pathlib import Path
import numpy as np
from numpy import genfromtxt
import glob
import random

from pred_cal_tmscore_FS import TM_score
from cal_plddt_ACFS import plddt_cal

# call calculating TM-scores of fs region
from pred_cal_tmscore_FS import *

class additional_prediction():
    def __init__(self, search_dir, pdb1, pdb1_name, pdb2, pdb2_name, pwd, size):


        add_pred_dir = pwd + 'additional_sampling/' + pdb1_name
        selection = size
        if selection == 0:
            print("max = 1, max-extra = 2")
        elif selection == 1:
            print("max = 2, max-extra = 4")
        elif selection == 2:
            print("max = 4, max-extra = 8")
        elif selection == 3:
            print("max = 8, max-extra = 16")
        elif selection == 4:
            print("max = 16, max-extra = 32")
        elif selection == 5:
            print("max = 32, max-extra = 64")
        elif selection == 6:
            print("max = 64, max-extra = 128")


        TMscore_add_all = []
        TMscore_fs_add_all = []

        #if os.path.isdir( add_pred_dir ):
        #    print("Predictions were already done")


        #else:
        #add_dir = 'additional_sampling' 
        #os.mkdir(add_dir)
        add_pred_sub_dir = 'additional_sampling/' + pdb1_name
        sub_remove_cmd = 'rm -rf ' + add_pred_sub_dir
        os.system(sub_remove_cmd)

        #random_seed = random.sample(range(100), 16)

        if selection == 0: ## max = 1, max-extra = 2
            command = 'colabfold_batch --amber --use-gpu-relax --num-seeds 20 --max-seq 1 --max-extra-seq 2 ' + search_dir + ' ' + pwd  + 'additional_sampling/' + pdb1_name
            print(command)
            os.system(command)

            #pred_dir = 'additional_sampling/' + pdb1_name + '/'
            #print(pred_dir)
            #MSA_additional_TMscore = TM_score(pred_dir, pdb1, pdb1_name, pdb2, pdb2_name)
            #TMscore_add_all = np.append(TMscore_add_all, MSA_additional_TMscore.tmscores)
            #print(TMscore_add_all)

            #MSA_additional_TMscor_fs = TM_score_fs(pred_dir, pdb1, pdb1_name, pdb2, pdb2_name)
            #TMscore_fs_add_all = np.append(TMscore_fs_add_all, MSA_additional_TMscor_fs.tmscores_fs)
            #print(TMscore_fs_add_all)

        elif selection == 1: ## max = 2, max-extra = 4
            command = 'colabfold_batch --amber --use-gpu-relax --num-seeds 20 --max-seq 2 --max-extra-seq 4 ' + search_dir + ' ' + pwd  + 'additional_sampling/' + pdb1_name
            print(command)
            os.system(command)

        elif selection == 2: ## max = 4, max-extra = 8
            command = 'colabfold_batch --amber --use-gpu-relax --num-seeds 20 --max-seq 4 --max-extra-seq 8 ' + search_dir + ' ' + pwd  + 'additional_sampling/' + pdb1_name
            print(command)
            os.system(command)

        elif selection == 3: ## max = 8, max-extra = 16
            command = 'colabfold_batch --amber --use-gpu-relax --num-seeds 20 --max-seq 8 --max-extra-seq 16 ' + search_dir + ' ' + pwd  + 'additional_sampling/' + pdb1_name
            print(command)
            os.system(command)

        elif selection == 4: ## max = 16, max-extra = 32
            command = 'colabfold_batch --amber --use-gpu-relax --num-seeds 20 --max-seq 16 --max-extra-seq 32 ' + search_dir + ' ' + pwd  + 'additional_sampling/' + pdb1_name
            print(command)
            os.system(command)

        elif selection == 5: ## max = 32, max-extra = 64
            command = 'colabfold_batch --amber --use-gpu-relax --num-seeds 20 --max-seq 32 --max-extra-seq 64 ' + search_dir + ' ' + pwd  + 'additional_sampling/' + pdb1_name
            print(command)
            os.system(command)

        elif selection == 6: ## max = 64, max-extra = 128
            command = 'colabfold_batch --amber --use-gpu-relax --num-seeds 20 --max-seq 64 --max-extra-seq 128 ' + search_dir + ' ' + pwd  + 'additional_sampling/' + pdb1_name
            print(command)
            os.system(command)


        pred_dir = 'additional_sampling/' + pdb1_name + '/'
        print(pred_dir)

        MSA_additional_TMscore = TM_score(pred_dir, pdb1, pdb1_name, pdb2, pdb2_name)
        TMscore_add_all = np.append(TMscore_add_all, MSA_additional_TMscore.tmscores)
        print(TMscore_add_all)

        MSA_additional_TMscor_fs = TM_score_fs(pred_dir, pdb1, pdb1_name, pdb2, pdb2_name)
        TMscore_fs_add_all = np.append(TMscore_fs_add_all, MSA_additional_TMscor_fs.tmscores_fs)
        print(TMscore_fs_add_all)

        #TMscore_add_all_reshape = TMscore_add_all.reshape(150, 5) ## for testing
        #np.savetxt('TMScore_additional-MSA_' + pdb1_name  + '.csv', TMscore_add_all_reshape, fmt='%2.3f')

        #TMscore_fs_add_all_reshape = TMscore_fs_add_all.reshape(150, 5)
        #np.savetxt('TMScore_fs_additional-MSA_' + pdb1_name  + '.csv', TMscore_fs_add_all_reshape, fmt='%2.3f')
        
        TMscore_add_all_reshape = TMscore_add_all.reshape(40, 5)
        np.savetxt('TMScore_additional-MSA_' + pdb1_name  + '.csv', TMscore_add_all_reshape, fmt='%2.3f')

        TMscore_fs_add_all_reshape = TMscore_fs_add_all.reshape(40, 5)
        np.savetxt('TMScore_fs_additional-MSA_' + pdb1_name  + '.csv', TMscore_fs_add_all_reshape, fmt='%2.3f')


