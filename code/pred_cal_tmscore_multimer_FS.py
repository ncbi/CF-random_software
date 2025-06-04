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

# call related modules of tmtools after installation
from tmtools import tm_align
from tmtools.io import get_structure, get_residue_data
from tmtools.testing import get_pdb_path

# call calculating TM-scores of fs region
from cal_tmscore_fs_only import * 
from cal_tmscore_fs_multimer import *

# call converting the multimer as a single chain structure
from convert_multi_single import *

# call converting the multimer as a separated chains
from split_multi_single import *


class TM_score_monomer():
    def __init__(self, pred_dir, pdb1_name, pdb2_name):
        
        ## loading reference pdb for TM-score
        pwd = os.getcwd() + '/'
        tmscores_monomer = []

        files_list = (glob.glob(str(pred_dir) + "/*_unrelaxed*pdb"))
        print(files_list)

        ##### pdb1_name part
        pdb1_dir = pwd + pdb1_name
        r2 = get_structure(get_pdb_path(str(pdb1_dir)))
        coords2, seq2 = get_residue_data(r2)

        if len(files_list) == 0:
            tmscores_monomer = [0.0, 0.0, 0.0, 0.0, 0.0]
            return tmscores_monomer

        for model in files_list:
            #modelpath = Path(model)
            #model  = str(modelpath.parent) + "/" + modelpath.stem
            model = model.replace('.pdb','')
            model = pwd + model
            s = get_structure(get_pdb_path(model))
            coords1, seq1 = get_residue_data(s)
            res = tm_align(coords1, coords2, seq1, seq2)
            tmscore = round(res.tm_norm_chain1,5) # wrt to model
            tmscores_monomer.append(tmscore)


        ##### pdb2_name part
        pdb2_dir = pwd + pdb2_name
        r3 = get_structure(get_pdb_path(str(pdb2_dir)))
        coords2, seq2 = get_residue_data(r3)

        if len(files_list) == 0:
            tmscores_monomer = [0.0, 0.0, 0.0, 0.0, 0.0]
            return tmscores_monomer

        for model in files_list:
            #modelpath = Path(model)
            #model  = str(modelpath.parent) + "/" + modelpath.stem
            #model = model.replace('_converted.pdb','_converted')
            model = model.replace('.pdb','')
            model = pwd + model
            s = get_structure(get_pdb_path(model))
            coords1, seq1 = get_residue_data(s)
            res = tm_align(coords1, coords2, seq1, seq2)
            tmscore = round(res.tm_norm_chain1,5) # wrt to model
            tmscores_monomer.append(tmscore)

        print(tmscores_monomer)
        self.tmscores_monomer = tmscores_monomer


class TM_score_multimer():
    def __init__(self, pred_dir, pdb1_name, pdb2_name):

        ## loading reference pdb for TM-score
        pwd = os.getcwd() + '/'
        tmscores_multimer = []

        ##### convert the multimer file as a single structure
        check_files_list = (glob.glob(str(pred_dir) + "/rmTER*_unrelaxed*.pdb"))
        #check_files_list = (glob.glob(str(pred_dir) + "/*_unrelaxed*pdb"))
        print(check_files_list)
        if not check_files_list:
            convert_m2s(pred_dir, pdb1_name, pdb2_name)
            files_list = (glob.glob(str(pred_dir) + "/rmTER*_unrelaxed*.pdb"))
            #files_list = (glob.glob(str(pred_dir) + "/*_unrelaxed*pdb"))
            print(files_list)
        else:
            files_list = (glob.glob(str(pred_dir) + "/rmTER*_unrelaxed*.pdb"))
            #files_list = (glob.glob(str(pred_dir) + "/*_unrelaxed*pdb"))
            print(files_list)


        ##### pdb2_name part
        pdb2_dir = pwd + pdb2_name + '_rmTER'
        r3 = get_structure(get_pdb_path(str(pdb2_dir)))
        coords2, seq2 = get_residue_data(r3)

        if len(files_list) == 0:
            tmscores_multimer = [0.0, 0.0, 0.0, 0.0, 0.0]
            return tmscores_multimer

        for model in files_list:
            #modelpath = Path(model)
            #model  = str(modelpath.parent) + "/" + modelpath.stem
            #model = model.replace('_converted.pdb','_converted')
            model = model.replace('.pdb','')
            model = pwd + model
            s = get_structure(get_pdb_path(model))
            coords1, seq1 = get_residue_data(s)
            res = tm_align(coords1, coords2, seq1, seq2)
            tmscore = round(res.tm_norm_chain1,5) # wrt to model
            tmscores_multimer.append(tmscore)

        print(tmscores_multimer)


        ##### pdb1_name part
        pdb1_dir = pwd + pdb1_name
        r2 = get_structure(get_pdb_path(str(pdb1_dir)))
        coords2, seq2 = get_residue_data(r2)

        if len(files_list) == 0:
            tmscores_multimer = [0.0, 0.0, 0.0, 0.0, 0.0]
            return tmscores_multimer

        for model in files_list:
            #modelpath = Path(model)
            #model  = str(modelpath.parent) + "/" + modelpath.stem
            model = model.replace('.pdb','')
            #model = model.replace('.pdb','')
            model = pwd + model
            s = get_structure(get_pdb_path(model))
            coords1, seq1 = get_residue_data(s)
            res = tm_align(coords1, coords2, seq1, seq2)
            tmscore = round(res.tm_norm_chain1,5) # wrt to model
            tmscores_multimer.append(tmscore)

        self.tmscores_multimer = tmscores_multimer






class CF_MSA_max():
    def __init__(self, search_dir, output_dir, pdb_name, rseed, num_seeds, model_type):

        command = 'colabfold_batch --num-seeds ' + str(num_seeds) + ' --model-type alphafold2_ptm --random-seed ' + str(rseed) + search_dir + output_dir
        print(command)
        os.system(command)
        



class CF_MSA_var():
    def __init__(self, pdb1_name, pdb2_name, search_dir, output_dir, rseed, num_seeds, model_type):
        #### shallow MSA section
        #### Global viarlable
        max_msa = 1; ext_msa = 2
        random_seed = rseed
        self.pdb1_name = pdb1_name; self.pdb2_name = pdb2_name

        for multi in (1, 2, 2, 2, 2, 2, 2):
            max_msa = int(max_msa * multi)
            ext_msa = int(ext_msa * multi)
        
            #### Colabfold part
            command = 'colabfold_batch --num-seeds ' + str(num_seeds) + ' --model-type ' + str(model_type) + ' --max-seq ' + str(max_msa) + ' --max-extra-seq ' + str(ext_msa) + search_dir + output_dir + str(random_seed) + '_max_' + str(max_msa) + '_ext_' + str(ext_msa)
            print(command); os.system(command)
        
       

    def cal_TM_score_multi(self, pdb1_name, pdb2_name, num_seeds, search_dir, output_dir, rseed, pdb1, pdb2):

        max_msa = 1; ext_msa = 2
        multi_size = 0; random_seed = rseed
        TMscore_multi = []; TMscore_multi_average = np.zeros((7, 1)) 
        TMscore_multi_fs = []; TMscore_multi_fs_average = np.zeros((7, 1))

        for multi in (1, 2, 2, 2, 2, 2, 2):
            max_msa = int(max_msa * multi)
            ext_msa = int(ext_msa * multi)

            fin_pred_dir = pdb1_name + '_predicted_models_rand_' + str(rseed) + '_max_' + str(max_msa) + '_ext_' + str(ext_msa) + '/'
            fin_pred_dir_all = pdb1_name + '_predicted_models_rand_' + str(rseed) + '_max_*'
            pred_files_list = (glob.glob(str(fin_pred_dir) + "/*_unrelaxed*pdb"))
     
            if len(pred_files_list) == 0:
                print("The TMscore list is empty")
                tmp = np.zeros((1, 25))
                TMscore_multi = np.append(TMscore_multi, tmp)
                TMscore_multi_fs = np.append(TMscore_multi_fs, tmp)
            else:
                run_TMscore_multi = TM_score_multimer(fin_pred_dir, pdb1_name, pdb2_name)
                TMscore_multi = np.append(TMscore_multi, run_TMscore_multi.tmscores_multimer); print(TMscore_multi)
                ##### for measuring the fold-switching region in multimer, just measure the TM-score of 
                ##### the first chain in between predicted and reference file 
                pdb2 = pdb2_name + '_rmTER.pdb'
                fin_fs_pred_dir = pdb1_name + '_predicted_models_rand_' + str(rseed) + '_max_' + str(max_msa) + '_ext_' + str(ext_msa) + '/'
                print(fin_fs_pred_dir)
                run_TMscore_multi_fs = TM_score_fs_multi(fin_fs_pred_dir, pdb1, pdb1_name, pdb2, pdb2_name) 
                TMscore_multi_fs = np.append(TMscore_multi_fs, run_TMscore_multi_fs.tmscores_fs ); print(TMscore_multi_fs)

        TMscore_multi = TMscore_multi.reshape(7 * 2, num_seeds * 5)
        np.savetxt('TMScore_random-MSA_' + pdb1_name  + '.csv', TMscore_multi, fmt='%2.3f')
        TMscore_multi = TMscore_multi[::2]
        TMscore_multi_fs = TMscore_multi_fs.reshape(7 * 2, num_seeds * 5)
        np.savetxt('TMScore_fs_random-MSA_' + pdb1_name  + '.csv', TMscore_multi_fs, fmt='%2.3f')
        TMscore_multi_fs = TMscore_multi_fs[1::2]

        print("TMscore multimer:"); print(TMscore_multi)
        print("TMscore fold-switching region in multimer:"); print(TMscore_multi_fs)

        if np.any(TMscore_multi > 0.4) and np.any(TMscore_multi_fs > 0.4):
            #for ii in range(0, int(TMscore_multi.shape[0] - 1)):
            tmp_cnt = 0
            for i in range(0, int(TMscore_multi.shape[0] - 1)):
            #for i in range(0, 13, 2):
                #TMscore_multi_average[ii] = np.average(TMscore_multi[ii])
                #TMscore_multi_fs_average[ii] = np.average(TMscore_multi_fs[ii])
                TMscore_multi_average[tmp_cnt] = np.average(TMscore_multi[i])
                TMscore_multi_fs_average[tmp_cnt] = np.average(TMscore_multi_fs[i])
                tmp_cnt = tmp_cnt + 1

            location = np.argmax(np.max(TMscore_multi_average, axis=1))
            print("The selected size of shallow random MSA is: ", np.argmax(np.max(TMscore_multi_fs_average, axis=1)))
            self.size_selection = int(location)

            mv_command = 'mv ' + fin_pred_dir_all + ' multimer_prediction/' + pdb1_name
            print(mv_command); os.system(mv_command)
     

        else:
            print("All calculated TMscores are not satisfying the creteria")
            print("All process is done.")
            mv_command = 'mv ' + fin_pred_dir_all + ' failed_prediction/'; os.system(mv_command)
            sys.exit()
                




class prediction_all_multimer_FS():
    def __init__(self, pdb1_name, pdb2_name, search_dir, nMSA, model_type, search_multi_dir, pdb1, pdb2):
        ### note: pdb1_name should be nomomer and pdb2_name should be multimer
        num_seeds = 5 + nMSA
        TER_count = 0
        pwd = os.getcwd() + '/'
        rm_converted_pdb = 'rm ' + pdb2_name + '_rmTER.pdb'; os.system(rm_converted_pdb)


        ##############################################################
        ##### Predicting all CF-random runs before calculate TM-scores
        ##### Predicting the monomer with deep MSA         
        #pre_random_seed = np.arange(0, 10, 1)
        pre_random_seed = random.sample(range(10), 1)
        random_seed_full_MSA = ''.join(map(str, pre_random_seed))
        output_dir = ' ' + pdb1_name + '_predicted_models_full_rand_' + str(random_seed_full_MSA)

        ##### Perform predction with full-length MSA
        MSA_full = CF_MSA_max(search_dir, output_dir, pdb1_name, random_seed_full_MSA, num_seeds, model_type)

        ##### Predicting the multimer with shallow random MSAs
        ##### check out varied-MSA with (msa-max: 1, 2, 4, 8, 16, 32, 64) (msa-extra: 2, 4, 8, 16, 32, 64, 128)
        output_dir = ' ' + pdb1_name + '_predicted_models_rand_'
        random_seed = random.sample(range(100), 1)
        random_seed = ''.join(map(str, random_seed))
        search_dir_update = ' ' + search_multi_dir.replace(' ','') + ' '

        MSA_var = CF_MSA_var(pdb1_name, pdb2_name, search_dir_update, output_dir, random_seed, num_seeds, model_type)
        #self.size_selection= MSA_var.size_selection


        ################################################################
        ##### Calculating all TM-scores including monomer and multimer
        ##### TM-score calculation for monoemr
        TMscore_monomer = []

        # Directory section
        gen_dir = 'multimer_prediction/' + pdb1_name

        if not os.path.exists(gen_dir):
            os.mkdir(gen_dir)

        pred_dir  = pdb1_name + '*predicted_models_full*'



        ##### Calculating the TM-score of fold-switching region
        ##### Extracting a signle chain from a multimer
        TER_count = 0
        with open(pdb2, 'r') as file:
            for line in file:
                TER = line.split()
                TER_count += TER.count("TER")
        
        line_cnt = 0
        #for i in range(0, TER_count):
        for i in range(0, 2):
            #output_file_name = pdb2_name[0:4] + '_multi.pdb'
            output_file_name = pdb2_name.split('_')[0] + '_multi.pdb'
        
            if line_cnt == 0:
                with open(pdb2, 'r') as infile, open(output_file_name, 'w') as outfile:
                    for line in infile:
                        outfile.write(line)
                        line_cnt = line_cnt + 1
                        if "TER " in line:
                            line_cnt = line_cnt + 1
                            break

        pdb2_name_multi = output_file_name.replace('.pdb','')


        ##### Calculate TM-score of monomer
        run_TMscore = TM_score_monomer(pred_dir, pdb1_name, pdb2_name)
        TMscore_monomer = np.array(run_TMscore.tmscores_monomer) 
        TMscore_monomer = TMscore_monomer.reshape(2, num_seeds * 5); print(TMscore_monomer)
        ##### Calculate TM-score of fold-switching region
        run_fs_TMscore = TM_score_fs(pred_dir, pdb1, pdb1_name, output_file_name, pdb2_name_multi)
        TMscore_monomer_fs = np.array(run_fs_TMscore.tmscores_fs)
        TMscore_monomer_fs = TMscore_monomer_fs.reshape(2, num_seeds * 5); print(TMscore_monomer_fs)

        if np.any(TMscore_monomer[0, :] >= 0.5) and np.any(TMscore_monomer_fs[0, :] >= 0.4):
            pred_dir = pdb1_name + '_predicted_models_full_rand_' + str(random_seed_full_MSA) + '/'
            mv_folder_cmd = 'mv ' + pred_dir + ' multimer_prediction/' + pdb1_name
            print(mv_folder_cmd); os.system(mv_folder_cmd)
            np.savetxt('TMScore_full-MSA_' + pdb1_name + '.csv', TMscore_monomer, fmt='%2.3f')

            MSA_var.cal_TM_score_multi(pdb1_name, pdb2_name_multi, num_seeds, search_dir_update, output_dir, random_seed, pdb1, output_file_name)
            #TMscore_multi_selection = MSA_var.cal_TM_score_multi(pdb1_name, pdb2_name_multi, num_seeds, search_dir_update, output_dir, random_seed, pdb1, output_file_name)
            print(MSA_var.size_selection); self.size_selection = MSA_var.size_selection



        else:
            pred_dir = pdb1_name + '_predicted_models*_rand_*/'
            mv_command = 'mv ' + pred_dir + ' failed_prediction/'; 
            print(mv_command); os.system(mv_command)
            print("Deep MSA cannot find the monomer")
            sys.exit()




        ##### TM-score calculation for multimer

