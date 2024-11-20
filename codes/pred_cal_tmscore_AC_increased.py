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




class TM_score_inc():
    def __init__(self, pred_dir, pdb1, pdb1_name, pdb2, pdb2_name):

        ## loading reference pdb for TM-score
        pwd = os.getcwd() + '/'
        tmscores = []
        
        #files_list = sorted(glob.glob(str(pred_dir) + "/*_relaxed*pdb"))
        files_list = (glob.glob(str(pred_dir) + "/*_unrelaxed*pdb"))
        print(files_list)


        ##### pdb1_name part
        pdb1_dir = pwd + pdb1_name
        r2 = get_structure(get_pdb_path(str(pdb1_dir)))
        coords2, seq2 = get_residue_data(r2)

        if len(files_list) == 0:
            tmscores = [0.0, 0.0, 0.0, 0.0, 0.0]
            return tmscores
        
        for model in files_list:
            #modelpath = Path(model)
            #model  = str(modelpath.parent) + "/" + modelpath.stem
            model = model.replace('.pdb','')
            model = pwd + model
            s = get_structure(get_pdb_path(model))
            coords1, seq1 = get_residue_data(s)
            res = tm_align(coords1, coords2, seq1, seq2)
            tmscore = round(res.tm_norm_chain1,5) # wrt to model
            tmscores.append(tmscore)

        print(tmscores)
        #print(tmscores[0:5])
        ##### pdb2_name part
        pdb2_dir = pwd + pdb2_name
        r3 = get_structure(get_pdb_path(str(pdb2_dir)))
        coords2, seq2 = get_residue_data(r3)


        for model in files_list:
            #modelpath = Path(model)
            #model  = str(modelpath.parent) + "/" + modelpath.stem
            model = model.replace('.pdb','')
            model = pwd + model
            s = get_structure(get_pdb_path(model))
            coords1, seq1 = get_residue_data(s)
            res = tm_align(coords1, coords2, seq1, seq2)
            tmscore = round(res.tm_norm_chain1,5) # wrt to model
            tmscores.append(tmscore)

        print(tmscores)
        self.tmscores = tmscores
        
        
        #if np.sum(tmscores[0:5]) > np.sum(tmscores[5:10]):
        #    self.ref_pdbname = pdb1_name
        #    self.alt_pdbname = pdb2_name
        #else:
        #    print("Full-length MSA is not matching with crystal structure (pdb1)")
        #    sys.exit()
        #    #self.ref_pdbname = pdb2_name
        #    #self.alt_pdbname = pdb1_name

            


class CF_MSA_max():
    def __init__(self, search_dir, output_dir, pdb_name, rseed):

        command = 'colabfold_batch --num-seeds 25 --random-seed ' + str(rseed) + search_dir + output_dir
        print(command)
        os.system(command)
        
        
    def additional_max(self, search_dir, output_dir, rseed, pdb1, pdb1_name, pdb2, pdb2_name, pwd):
        ### run additional full-MSA prediction with different random number
        gen_dir = 'failed_prediction/' + pdb1_name

        if not os.path.exists(gen_dir):
            os.mkdir(gen_dir)
            
            
        mv_folder_cmd = 'mv ' + pdb1_name + '_predicted_models_full_rand_' + rseed + ' failed_prediction/' + pdb1_name
        print(mv_folder_cmd); os.system(mv_folder_cmd)
        
        
        print("Full-MSA prediction is not tightly aligned to crystal structure")
        print("putting additional 15 random seeds")
        print("               ")
        add_trial = random.sample(range(17, 100), 5)
        print(add_trial)
        
        
        for add_rseed in add_trial:
            
            output_dir = ' ' + pdb1_name + '_predicted_models_full_rand_' + str(add_rseed)
            MSA_full = CF_MSA_max(search_dir, output_dir, pdb1_name, add_rseed)
            
            pred_dir = pdb1_name + '_predicted_models_full_rand_' + str(add_rseed) + '/'
            MSA_full_TMscore = TM_score_inc(pred_dir, pdb1, pdb1_name, pdb2, pdb2_name)
            
            
            ##### check-out TM-scores of prediction with full-length-MSA (whole)
            pred_dir = pdb1_name + '_predicted_models_full_rand_' + str(add_rseed) + '/'
            MSA_full_TMscore = TM_score_inc(pred_dir, pdb1, pdb1_name, pdb2, pdb2_name)
            full_TMscore = np.array(MSA_full_TMscore.tmscores)
            full_TMscore = full_TMscore.reshape(2, 125)
            


            if np.any(full_TMscore[0, :] > 0.5):
                ref_name = pdb1_name; alt_name = pdb2_name
                self.ref_name = ref_name; self.alt_name = alt_name
                break
            elif np.any(full_TMscore[1, :] > 0.5):
                ref_name = pdb2_name; alt_name = pdb1_name
                self.ref_name = ref_name; self.alt_name = alt_name
                break
            else:
                print("Full-MSA prediction is not tightly aligned to crystal structure with additional 15 seeds")
                mv_folder_cmd = 'mv ' + pdb1_name + '_predicted_models_full_rand_' + str(add_rseed) + ' failed_prediction/' + pdb1_name
                print(mv_folder_cmd); os.system(mv_folder_cmd)
                print("                  ")


        print("                  ")
        print("Predictions using full-length MSA with additional randon-seed are all failed")
        print("Prediction is done")
        sys.exit()



class CF_MSA_var():
    def __init__(self, pdb1, pdb1_name, pdb2, pdb2_name, search_dir, output_dir, rseed, ref_name, alt_name):
        #### shallow MSA section
        #### Global viarlable
        max_msa = 1
        ext_msa = 2
        random_seed = np.array(rseed) ## needed to remove future

        self.pdb1_name = pdb1_name
        self.alt_name = alt_name
        self.ref_name = ref_name

        TMscores_random = [] ## whole structure


        for ran_seed in random_seed:
            max_msa = 1
            ext_msa = 2

            TMscores_random = []

            for multi in (1, 2, 2, 2, 2, 2, 2):
                max_msa = max_msa * multi
                ext_msa = ext_msa * multi

                #### Colabfold part
                command = 'colabfold_batch --num-seed 25 --max-seq ' + str(max_msa) + ' --max-extra-seq ' + str(ext_msa) + search_dir + output_dir + str(ran_seed) + '_max_' + str(max_msa) + '_ext_' + str(ext_msa)
                print(command)
                os.system(command)


                #### TMscore whole part
                pred_dir = pdb1_name + '_predicted_models_rand_' + str(ran_seed) + '_max_' + str(max_msa) + '_ext_' + str(ext_msa) + '/'
                MSA_shallow_TMscore = TM_score_inc(pred_dir, pdb1, pdb1_name, pdb2, pdb2_name)

                TMscores_random = np.append(TMscores_random, MSA_shallow_TMscore.tmscores)
                print(TMscores_random)




            self.TMscores_random = TMscores_random
            fin_pred_dir = pdb1_name + '_predicted_models_rand_' + str(ran_seed) + '_max_*'

            TMscores_random_reshape = self.TMscores_random.reshape(14, 125)

            print(TMscores_random_reshape)

            TMscores_random_alter = np.zeros((7, 125))

            #### finding alternative pdb_name


            if alt_name == pdb2_name:            
                #for i in 1, 3, 5, 7, 9, 11, 13 in TM_scores:
                tmp_cnt = 0
                for i in range(1, 14, 2):
                    print(TMscores_random_reshape[i, :])
                    TMscores_random_alter[tmp_cnt, :] = TMscores_random_reshape[i, :]
                    tmp_cnt = tmp_cnt + 1
            else:
                #for i in 0, 2, 4, 6, 8, 10, 12 in TM_scores:
                tmp_cnt = 0
                for i in range(0, 13, 2):
                    print(TMscores_random_reshape[i, :])
                    TMscores_random_alter[tmp_cnt, :] = TMscores_random_reshape[i, :]
                    tmp_cnt = tmp_cnt + 1




            print("                   ")
            print("Confirming the TM-score with alternative conformation is good or not")
            print(TMscores_random_alter)
            print("                   ")
            

            if np.any(TMscores_random_reshape > 0.5):
                # save all TM-scores from random MSA (1-2, 2-4, 4-8.... in order)
                #TMscores_random_reshape = TMscores_random.reshape(14, 5)
                np.savetxt('TMScore_random-MSA_' + pdb1_name  + '.csv', TMscores_random_alter, fmt='%2.3f')
                
                
                gen_dir = 'successed_prediction/' + pdb1_name
                if not os.path.exists(gen_dir):
                    os.makedirs(gen_dir)
                    mv_command = 'mv ' + fin_pred_dir + ' successed_prediction/' + pdb1_name
                    print(mv_command)
                    os.system(mv_command)
                    #self.TMscores_random_alter = TMscores_random_alter
                    self.TMscores_random_alter = TMscores_random_reshape
                    break
                else:
                    mv_command = 'mv ' + fin_pred_dir + ' successed_prediction/' + pdb1_name
                    print(mv_command)
                    os.system(mv_command)
                    #self.TMscores_random_alter = TMscores_random_alter
                    self.TMscores_random_alter = TMscores_random_reshape
                    break

            else:
                gen_dir = 'failed_prediction/' + pdb1_name
                if not os.path.exists(gen_dir):
                    os.makedirs(gen_dir)
                    mv_command = 'mv ' + fin_pred_dir + ' failed_prediction/' + pdb1_name
                    print(mv_command)
                    os.system(mv_command)
                    
                    print("Full-MSA prediction is not tightly aligned to crystal structure with additional seeds")
                    print("Predcition is done")

                else:
                    mv_command = 'mv ' + fin_pred_dir + ' failed_prediction/' + pdb1_name
                    print(mv_command)
                    os.system(mv_command)

                    print("Full-MSA prediction is not tightly aligned to crystal structure with additional seeds")
                    print("Predcition is done")



            self.TMscores_random_alter = TMscores_random_reshape


            #self.TMscores_random_alter = TMscores_random_alter
            if ran_seed == random_seed[-1]:
                print("Can't find the size of random shallow MSA with additional random seeds")
                print("Prediction is done")
                sys.exit()


    def select_size(self, TMscores_random_alter, pdb1_name, pdb2_name, alt_name):
        
        TMscores_random_reshape = TMscores_random_alter.reshape(14, 125)
        TMscores_random_locat = np.zeros((7, 125))
        
        #### finding locatnative pdb_name
        
        if alt_name == pdb2_name:
            #for i in 1, 3, 5, 7, 9, 11, 13 in TM_scores:
            tmp_cnt = 0
            for i in range(1, 14, 2):
                print(TMscores_random_reshape[i, :])
                TMscores_random_locat[tmp_cnt, :] = TMscores_random_reshape[i, :]
                tmp_cnt = tmp_cnt + 1
        else:
            #for i in 0, 2, 4, 6, 8, 10, 12 in TM_scores:
            tmp_cnt = 0
            for i in range(0, 13, 2):
                print(TMscores_random_reshape[i, :])
                TMscores_random_locat[tmp_cnt, :] = TMscores_random_reshape[i, :]
                tmp_cnt = tmp_cnt + 1


        TMscore_data = TMscores_random_locat
        TMscore_data = TMscores_random_locat.reshape(7, 125)
        TMscore_data_sum = np.zeros((7, 1))


        
        print((TMscore_data.shape[0]))


        for ii in range(0, int(TMscore_data.shape[0])):
            TMscore_data_sum[ii] = np.sum(TMscore_data[ii])

        
        location = np.argmax(np.max(TMscore_data_sum, axis=1))
        print(np.argmax(np.max(TMscore_data_sum, axis=1)))





        print("Selecting...")

        TMscore_data = TMscores_random_alter
        TMscore_data = TMscores_random_alter.reshape(14, 125)


        location_org = location


        if alt_name == pdb2_name:
            location = (location * 2) + 1
        else:
            location = (location * 2) 




        
        if alt_name == pdb2_name and np.any(TMscore_data[location, :] >= 0.5):
            print(TMscore_data[location, :])
            selection = int((location - 1) / 2)
            self.selection = selection

        elif alt_name == pdb1_name and np.any(TMscore_data[location, :] >= 0.5):
            print(TMscore_data[location, :])
            selection = int(location / 2); 
            self.selection = selection

        else:
            print("Predictions are bad")
            print("Predictions of whole structure are bad")
            rm_folder_cmd = 'rm -rf successed_prediction/' + self.pdb1_name + '/'
            print(rm_folder_cmd)
            os.system(rm_folder_cmd)
            sys.exit()




class prediction_all_AC_inc():
    def __init__(self, pdb1, pdb1_name, pdb2, pdb2_name, search_dir):
        
        pre_random_seed = np.random.randint(0, 16, 1)
        random_seed = ''.join(map(str, pre_random_seed))
        print(random_seed)
        output_dir = ' ' + pdb1_name + '_predicted_models_full_rand_' + str(random_seed)
        print(output_dir)


        ##### Perform predction with full-length MSA
        MSA_full = CF_MSA_max(search_dir, output_dir, pdb1_name, random_seed)

        pwd = os.getcwd() + '/'

        ##### check-out TM-scores of prediction with full-length-MSA (whole)
        pred_dir = pdb1_name + '_predicted_models_full_rand_' + str(random_seed) + '/'
        MSA_full_TMscore = TM_score_inc(pred_dir, pdb1, pdb1_name, pdb2, pdb2_name)
        full_TMscore = np.array(MSA_full_TMscore.tmscores)
        full_TMscore = full_TMscore.reshape(2, 125)
        #print(full_TMscore[0, :]) ;        #print(full_TMscore[1, :])


        print("TM-score of fold1")
        print(full_TMscore[0, :])
        print("                ")
        print("TM-score of fold2")
        print(full_TMscore[1, :])
        print("                ")

        
        ##### check-out the 1st prediction results are good or not
        if np.average(full_TMscore[0, :]) > np.average(full_TMscore[1, :]):
            if np.any(full_TMscore[0, :] > 0.5):
                ref_name = pdb1_name; alt_name = pdb2_name
            else:
                print("Performing predictions with full-MSA using additional 5 seeds")
                MSA_full_add = MSA_full.additional_max(search_dir, output_dir,random_seed, pdb1, pdb1_name, pdb2, pdb2_name, pwd)
                ref_name = MSA_full_add.ref_name; alt_name = MSA_full_add.alt_name

        else:
            if np.any(full_TMscore[1, :] > 0.5): 
                ref_name = pdb2_name; alt_name = pdb1_name
            elif np.any(full_TMscore[0, :] > 0.5):
                ref_name = pdb1_name; alt_name = pdb2_name
            else:
                print("Performing predictions with full-MSA using additional 5 seeds")
                MSA_full_add = MSA_full.additional_max(search_dir, output_dir,random_seed, pdb1, pdb1_name, pdb2, pdb2_name, pwd)
                ref_name = MSA_full_add.ref_name; alt_name = MSA_full_add.alt_name





        print("Reference structure: ", ref_name)
        print("Alternative structure: ", alt_name)

        # save TM-score from full-length MSA
        np.savetxt('TMScore_full-MSA_' + pdb1_name + '.csv', full_TMscore, fmt='%2.3f')
        # save TM-score of fold-switching region from full-length MSA

        # Directory section
        gen_dir = 'successed_prediction/' + pdb1_name

        if not os.path.exists(gen_dir):
            os.mkdir(gen_dir)

        mv_folder_cmd = 'mv ' + pred_dir + ' successed_prediction/' + pdb1_name
        print(mv_folder_cmd); os.system(mv_folder_cmd)

        print("Full-MSA prediction is tightly aligned to crystal structure"); print("               ")

        ##### check out varied-MSA with (msa-max: 1, 2, 4, 8, 16, 32, 64) (msa-extra: 2, 4, 8, 16, 32, 64, 128)
        output_dir = ' ' + pdb1_name + '_predicted_models_rand_'
        random_seed = random.sample(range(100), 5)
        MSA_var = CF_MSA_var(pdb1, pdb1_name, pdb2, pdb2_name, search_dir, output_dir, random_seed, ref_name, alt_name)


        print("                                     ")
        print("Finding optimal size of ramdon MSA...")
        MSA_var.select_size(MSA_var.TMscores_random_alter, pdb1_name, pdb2_name, alt_name)
        #size_selection = MSA_var.selection
        #print(size_selection)
        #self.size_selection = size_selection

        
        if np.all(MSA_var.TMscores_random_alter) < 0.5:
            print("All predictions are failed")
            sys.exit()
    
        else:
            size_selection = MSA_var.selection
            print(size_selection)
            self.size_selection = size_selection

