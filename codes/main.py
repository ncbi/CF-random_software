#!/upyMolsr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  9 14:51:00 2024
 
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
import argparse


from pred_cal_tmscore_FS import *
from pred_cal_tmscore_blind import *
from pred_cal_tmscore_AC import *
from additional_pred_FS import * 
from additional_pred_AC import *
from foldseek_run import *
from CF_random_blind import *
from cal_plddt_ACFS import *
from PLOT_AC import *
from PLOT_FS import *

if __name__ == "__main__":

    import warnings
    warnings.filterwarnings('ignore')

    ######################################################################################################
    ###### initiallization pdb format (removing HETATM)
    #os.system("for i in *pdb;do echo $i;sed -i '/HETATM/d' $i;done")



    ######################################################################################################
    ###### initiallization and input
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb1", type=str, help='PDB structure for the target crystal structure (target to be predicted)')
    parser.add_argument("--pdb2", type=str, help='PDB structure for the alternative crystal structure')
    parser.add_argument("--fname", type=str, help='put folder name after colabsearch' )
    parser.add_argument("--pname", type=str, help='temporary name for predicting blind mode' )
    parser.add_argument("--option", type=str, help='select prediction mode AC and FS e.g. AC = alterantive conformation or FS = fold-switching')
    args = parser.parse_args()



    if args.option == "blind":
        if args.pdb1 is None:
            pdb1_name = args.pname
            print("work name:", pdb1_name)
        else:
            pdb1_name = args.fname
            pdb1_name = pdb1_name.replace('/','')
            print("work name:", pdb1_name)
    elif args.pdb1 is None:
        pdb1_name = args.fname
    else:
        pdb1 = args.pdb1; pdb2 = args.pdb2
        pdb1_name = pdb1.replace('.pdb','');  pdb2_name = pdb2.replace('.pdb','')



    pwd = os.getcwd() + '/'
    search_dir = ' ' + pwd + args.fname

    blind = 'blind_prediction'
    success = 'successed_prediction'
    fail = 'failed_prediction'



    if not os.path.exists(success):
        os.mkdir(success)
    else:
        succ_dir_count = 0
        for root_dir, cur_dir, files in os.walk(pwd + success + '/' + pdb1_name + '/'):
            succ_dir_count += len(cur_dir)
        
    if os.path.exists(success + '/' + pdb1_name):
        if succ_dir_count >= 8:
            print("Prediction was already done")
        else:
            print("Folder is already created and cleaning existed subfolders")
            rm_pre_folders = 'rm -rf ' + success + '/' + pdb1_name + '/'
            os.system(rm_pre_folders)
    else:
        pass



    if not os.path.exists(fail):
        os.mkdir(fail)
    else:
        fail_dir_count = 0
        for root_dir, cur_dir, files in os.walk(pwd + fail + '/' + pdb1_name + '/'):
            fail_dir_count += len(cur_dir)

    if os.path.exists(fail + '/' + pdb1_name):
        if fail_dir_count >= 8:
            print("Prediction was already done")
        else:
            print("Folder is already created and cleaning existed subfolders")
            rm_pre_folders = 'rm -rf ' + fail + '/' + pdb1_name + '/'
            os.system(rm_pre_folders)
    else:
        pass








    if args.option == "AC":
        print("Predicting alternative conformations")
        ######################################################################################################
        ###### running prediction using full- and shallow random-MSA
        if os.path.exists(success + '/' + pdb1_name) and succ_dir_count >= 8:
            print("Predictions including full- and random-MSA were already done")
        else:
            pred_1st_all = prediction_all_AC(pdb1, pdb1_name, pdb2, pdb2_name, search_dir)
            shallow_MSA_size = []
            shallow_MSA_size = np.append(shallow_MSA_size, pred_1st_all.size_selection)
            print("               ")
            print("Specific size of shallow random MSA is similar to full-MSA")
            print(shallow_MSA_size)
            np.savetxt('selected_MSA-size_' + pdb1_name + '.csv', shallow_MSA_size)

        ######################################################################################################
        ###### additional prediction based on selected size 
        add_dir = 'additional_sampling/'
        shallow_MSA_size = genfromtxt("selected_MSA-size_" + pdb1_name + ".csv")

        add_pred_dir = add_dir + pdb1_name

        if os.path.isdir( add_pred_dir ):
            print("Checking predictions of ensembles after selection")
            num_files = len([f for f in os.listdir(add_pred_dir)if os.path.isfile(os.path.join(add_pred_dir, f))])
            print(num_files)
            
            if num_files > 200:
                print("Predictions were already done")
            else:
                additional_prediction_AC(search_dir, pdb1, pdb1_name, pdb2, pdb2_name, pwd, shallow_MSA_size)

        else:
            additional_prediction_AC(search_dir, pdb1, pdb1_name, pdb2, pdb2_name, pwd, shallow_MSA_size)

        ######################################################################################################
        ##### calculate plddt including initial and additional predictions
        if os.path.exists( pwd + add_dir ):
            list_org_samplings = glob.glob( str(pwd) + str(success) + '/' + str(pdb1_name) + '/*full_rand*/')
            list_ran_samplings = glob.glob( str(pwd) + str(success) + '/' + str(pdb1_name) + '/*max*/')
            list_add_samplings = glob.glob( str(pwd) + str(add_dir) + str(pdb1_name))

            full = 'full-MSA'
            random = 'random-MSA'
            addition = 'additional-MSA'
            plddt_cal(list_org_samplings, full, pdb1_name)
            plddt_cal(list_ran_samplings, random, pdb1_name)
            plddt_cal(list_add_samplings, addition, pdb1_name)

        ######################################################################################################
        ##### plot the 2D-scatter plot of TM-scores with pLDDT
        plot_2D_scatter_AC(full, random, addition, pdb1, pdb1_name, pdb2, pdb2_name)









    elif args.option == "FS":
        print("Predicting fold-swithcing models")
        ######################################################################################################
        ###### running prediction using full- and shallow random-MSA
        if os.path.exists(success + '/' + pdb1_name) and succ_dir_count >= 8:
            print("Predictions including full- and random-MSA were already done")
        else:
            pred_1st_all = prediction_all(pdb1, pdb1_name, pdb2, pdb2_name, search_dir)
            shallow_MSA_size = []
            shallow_MSA_size = np.append(shallow_MSA_size, pred_1st_all.size_selection)
            print("               ")
            print("Specific size of shallow random MSA is similar to full-MSA")
            print(shallow_MSA_size)
            np.savetxt('selected_MSA-size_' + pdb1_name + '.csv', shallow_MSA_size)
    
    
        ######################################################################################################
        ###### additional prediction based on selected size 
        add_dir = 'additional_sampling/'
        shallow_MSA_size = genfromtxt("selected_MSA-size_" + pdb1_name + ".csv")
    
        add_pred_dir = add_dir + pdb1_name
    
    
    
        if os.path.isdir( add_pred_dir ):
            print("Checking predictions of ensembles after selection")
            num_files = len([f for f in os.listdir(add_pred_dir)if os.path.isfile(os.path.join(add_pred_dir, f))])
            print(num_files)
    
            if num_files > 200:
                print("Predictions were already done")
            else:
                additional_prediction(search_dir, pdb1, pdb1_name, pdb2, pdb2_name, pwd, shallow_MSA_size)
        else:
            additional_prediction(search_dir, pdb1, pdb1_name, pdb2, pdb2_name, pwd, shallow_MSA_size)
    
        ######################################################################################################
        ##### calculate plddt including initial and additional predictions
        if os.path.exists( pwd + add_dir ):
            list_org_samplings = glob.glob( str(pwd) + str(success) + '/' + str(pdb1_name) + '/*full_rand*/')
            list_ran_samplings = glob.glob( str(pwd) + str(success) + '/' + str(pdb1_name) + '/*max*/') 
            list_add_samplings = glob.glob( str(pwd) + str(add_dir) + str(pdb1_name))
            #print(list_org_samplings)
            #print(list_add_samplings)
            full = 'full-MSA'
            random = 'random-MSA'
            addition = 'additional-MSA'
            plddt_cal(list_org_samplings, full, pdb1_name)
            plddt_cal(list_ran_samplings, random, pdb1_name)
            plddt_cal(list_add_samplings, addition, pdb1_name)
    
    
        ######################################################################################################
        ##### plot the 2D-scatter plot of TM-scores with pLDDT
        plot_2D_scatter(full, random, addition, pdb1, pdb1_name, pdb2, pdb2_name)









    elif args.option == "blind":
        print("Predicting fold-swithcing proteins without crystal structures of pdbs")
        ######################################################################################################
        ###### check previous predictions were performed or not
        if not os.path.exists(blind):
            os.mkdir(blind)
        else:
            blind_dir_count = 0
            for root_dir, cur_dir, files in os.walk(pwd + blind + '/' + pdb1_name + '/'):
                blind_dir_count += len(cur_dir)
            
        if os.path.exists(blind + '/' + pdb1_name):
            if blind_dir_count >= 8:
                print("Prediction was already done")
            else:
                print("Folder is already created and cleaning existed subfolders")
                rm_pre_folders = 'rm -rf ' + blind + '/' + pdb1_name + '/'
                os.system(rm_pre_folders)
        else:
            pass



        ###### running prediction using full- and shallow random-MSA
        blind_pred_path = 'blind_prediction/' + pdb1_name
        print(blind_pred_path)

        if os.path.exists(blind + '/' + pdb1_name) and blind_dir_count >= 8:
            print("Predictions including full- and random-MSA were already done")


            fseek_file_count = 0
            for root_dir, cur_dir, files in os.walk(pwd + blind + '/' + pdb1_name + '/'):
                fseek_file_count += len(files)

            print(fseek_file_count)
            #if fseek_file_count == 856: ##(107 * 8) 107 includes foldseek file and 8 means the numbers of prediction folders
            if fseek_file_count >= 640: ##672
                print("    "); print("Foldseek search was done")
                #### performing the PCA calculation with RMSD
                #PCA_rmsd(pdb1_name, blind_pred_path)
                blind_screening(pdb1_name, blind_pred_path)
            else:
                running_foldseek_all(pdb1_name)

                #### performing the PCA calculation with RMSD
                #PCA_rmsd(pdb1_name, blind_pred_path)
                blind_screening(pdb1_name, blind_pred_path)



        else:
            prediction_all_blind(pdb1_name, search_dir)
            print("               ")
            print("Finished running for prediction using full- and shallow random-MSAs")
            
            print("               ")
            print("Running Foldseek to find the relatedcrystal structures")
            running_foldseek_all(pdb1_name)

            #### performing the PCA calculation with RMSD
            #PCA_rmsd(pdb1_name, blind_pred_path)
            blind_screening(pdb1_name, blind_pred_path)





    else:
        print("Please type correct option")
