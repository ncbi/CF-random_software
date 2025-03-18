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
from cal_plddt_ACFS import *
from PLOT_AC import *
from PLOT_FS import *
from search_w_foldseek_cluster import *

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
    parser.add_argument("--fname", type=str, help='put MSA folder name after colabsearch' )
    parser.add_argument("--fmname", type=str, help='put multimer MSA folder name after colabsearch' )
    parser.add_argument("--pname", type=str, help='job name for predicting blind mode' )
    parser.add_argument("--nMSA", type=str, help='number of samples for predicting the structure with MSA')
    parser.add_argument("--nENS", type=str, help='number of samples for predicting the structure for ensemble generation')
    parser.add_argument("--option", type=str, help='select prediction mode inAC, AC and FS e.g. AC = alterantive conformation or FS = fold-switching or inAC = increased sampling for predicting alternative conformation')
    parser.add_argument("--type", type=str, help='select model-type of Colabfold e.g. ptm, monomer, and , multimer')
    args = parser.parse_args()



    blind = 'blind_prediction'
    success = 'successed_prediction'
    fail = 'failed_prediction'
    multi = 'multimer_prediction'
    pwd = os.getcwd() + '/'


    if args.option == "blind":
        if args.pdb1 is None:
            pdb1_name = args.pname
            print("work name:", pdb1_name)
        elif args.pdb1 is None and args.pname is None:
            pdb1_name = args.fname
            pdb1_name = pdb1_name.replace('/','')
            print("work name:", pdb1_name)
        else:
            pdb1_name = args.fname
            pdb1_name = pdb1_name.replace('/','')
            print("work name:", pdb1_name)
    elif args.pdb1 is None:
        pdb1_name = args.fname
    elif args.pdb1 is not None and args.pdb2 is not None:
        pdb1 = args.pdb1; pdb2 = args.pdb2
        pdb1_name = pdb1.replace('.pdb','');  pdb2_name = pdb2.replace('.pdb','')
        print(pdb1_name, pdb2_name)


    #if int(args.nMSA) == 0 and int(args.nENS) == 0:
    if args.nMSA is None and args.nENS is None: 
        nMSA = 0; nENS = 0;
    elif args.nMSA is not None and args.nENS is not None: 
        nMSA = int(args.nMSA); nENS = int(args.nENS)
    elif args.nMSA is None and args.nENS is not None: 
        nMSA = 0; nENS = int(args.nENS)
    elif args.nMSA is not None and args.nENS is None:
        nMSA = int(args.nMSA); nENS = 0
    else:
        print("Please put correct option of nMSA or nENS")
        exit()


    if args.fname is None and args.fmname is None:
        print("Please put MSA folder and file for monomer prediction")
        sys.exit()
    elif args.fname is None and args.fmname is not None:
        print("Please put MSA folder and file for monomer prediction")
        sys.exit()
    elif args.fname is not None and args.fmname is None:
        search_dir = ' ' + pwd + args.fname; search_multi_dir = 0
    elif args.fname is not None and args.fmname is not None:
        search_dir = ' ' + pwd + args.fname; search_multi_dir = ' ' + pwd + args.fmname;




    ### model-type identification
    model_type = []
    if args.type is None or args.type == "ptm":
        model_type = "alphafold2_ptm"
    elif args.type == "monomer":
        model_type = "alphafold2"
    elif args.type == "multimer":
        ### check how many chains in a multimer
        TER_count = 0 
        with open(pdb1, 'r') as file:
            for line in file:
                TER = line.split()
                TER_count += TER.count("TER")

            print(TER_count, " of chains in this multimer file.")
            model_type = "alphafold2_multimer_v3"

            if not os.path.exists(multi):
                os.mkdir(multi)

    else:
        print("Please put correct model-type option")
        exit()



    pwd = os.getcwd() + '/'
    search_dir = ' ' + pwd + args.fname





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



    if not os.path.exists(multi):
        os.mkdir(multi)
    else:
        succ_dir_count = 0
        for root_dir, cur_dir, files in os.walk(pwd + multi + '/' + pdb1_name + '/'):
            succ_dir_count += len(cur_dir)

    if os.path.exists(multi + '/' + pdb1_name):
        if succ_dir_count >= 8:
            print("Prediction was already done")
        else:
            print("Folder is already created and cleaning existed subfolders")
            rm_pre_folders = 'rm -rf ' + multi + '/' + pdb1_name + '/'
            os.system(rm_pre_folders)
    else:
        pass






    if args.option == "AC":
        print("Predicting alternative conformations")
        ######################################################################################################
        ###### running prediction using full- and shallow random-MSA
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





        if os.path.exists(success + '/' + pdb1_name) and succ_dir_count >= 8:
            print("Predictions including full- and random-MSA were already done")
        elif os.path.exists(multi + '/' + pdb1_name) and succ_dir_count >= 8:
            print("Predictions including full- and random-MSA were already done")
        else:

            rm_pre_folders = 'rm -rf ' + success + '/' + pdb1_name + '/'; os.system(rm_pre_folders)
            rm_pre_folders = 'rm -rf ' + multi + '/' + pdb1_name + '/'  ; os.system(rm_pre_folders)


            pred_1st_all = prediction_all_AC(pdb1, pdb1_name, pdb2, pdb2_name, search_dir, nMSA, model_type, search_multi_dir)
            shallow_MSA_size = []
            shallow_MSA_size = np.append(shallow_MSA_size, pred_1st_all.size_selection)
            print("               ")
            print("Specific size of shallow random MSA is similar to full-MSA")
            print(shallow_MSA_size)
            np.savetxt('selected_MSA-size_' + pdb1_name + '.csv', shallow_MSA_size)


        ######################################################################################################
        ##### calculate plddt of initial predictions
        if model_type == "alphafold2_multimer_v3":
            list_org_samplings = glob.glob( str(pwd) + str(multi) + '/' + str(pdb1_name) + '/*full_rand*/')
            list_ran_samplings = glob.glob( str(pwd) + str(multi) + '/' + str(pdb1_name) + '/*max*/')

            full = 'full-MSA'; random = 'random-MSA' ; 
            plddt_cal(list_org_samplings, full, pdb1_name, nMSA, nENS, model_type)
            plddt_cal(list_ran_samplings, random, pdb1_name, nMSA, nENS, model_type)

        else:
            list_org_samplings = glob.glob( str(pwd) + str(success) + '/' + str(pdb1_name) + '/*full_rand*/')
            list_ran_samplings = glob.glob( str(pwd) + str(success) + '/' + str(pdb1_name) + '/*max*/')

            full = 'full-MSA'; random = 'random-MSA' ; 
            plddt_cal(list_org_samplings, full, pdb1_name, nMSA, nENS, model_type)
            plddt_cal(list_ran_samplings, random, pdb1_name, nMSA, nENS, model_type)

        ######################################################################################################
        ##### plot the 2D-scatter plot of TM-scores with pLDDT
        plot_2D_scatter_AC(full, random, pdb1, pdb1_name, pdb2, pdb2_name, nMSA, nENS, model_type)






    elif args.option == "FS":
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


        print("Predicting fold-swithcing models")
        ######################################################################################################
        ###### running prediction using full- and shallow random-MSA
        if os.path.exists(success + '/' + pdb1_name) and succ_dir_count >= 8:
            print("Predictions including full- and random-MSA were already done")
        elif os.path.exists(multi + '/' + pdb1_name) and succ_dir_count >= 8:
            print("Predictions including full- and random-MSA were already done")
        else:
            pred_1st_all = prediction_all(pdb1, pdb1_name, pdb2, pdb2_name, search_dir, nMSA, model_type, search_multi_dir)
            shallow_MSA_size = []
            if args.type != "multimer":
                shallow_MSA_size = np.append(shallow_MSA_size, pred_1st_all.size_selection)
            else:
                shallow_MSA_size = np.append(shallow_MSA_size, pred_1st_all.size_selection)
            print("               ")
            print("Specific size of shallow random MSA is similar to full-MSA")
            print(shallow_MSA_size)
            np.savetxt('selected_MSA-size_' + pdb1_name + '.csv', shallow_MSA_size)
    
    
    
        ######################################################################################################
        ##### calculate plddt of initial predictions
        if model_type == "alphafold2_multimer_v3":
            list_org_samplings = glob.glob( str(pwd) + str(multi) + '/' + str(pdb1_name) + '/*full_rand*/')
            list_ran_samplings = glob.glob( str(pwd) + str(multi) + '/' + str(pdb1_name) + '/*max*/')

            full = 'full-MSA'; random = 'random-MSA' ; 
            plddt_cal(list_org_samplings, full, pdb1_name, nMSA, nENS, model_type)
            plddt_cal(list_ran_samplings, random, pdb1_name, nMSA, nENS, model_type)

        else:
            list_org_samplings = glob.glob( str(pwd) + str(success) + '/' + str(pdb1_name) + '/*full_rand*/')
            list_ran_samplings = glob.glob( str(pwd) + str(success) + '/' + str(pdb1_name) + '/*max*/')

            full = 'full-MSA'; random = 'random-MSA' ; 
            plddt_cal(list_org_samplings, full, pdb1_name, nMSA, nENS, model_type)
            plddt_cal(list_ran_samplings, random, pdb1_name, nMSA, nENS, model_type)
        
        
        
        
        ######################################################################################################
        ##### plot the 2D-scatter plot of TM-scores with pLDDT
        if model_type == "alphafold2_multimer_v3":
            plot_2D_scatter_AC(full, random, pdb1, pdb1_name, pdb2, pdb2_name, nMSA, nENS, model_type)
        else:
            plot_2D_scatter(full, random, pdb1, pdb1_name, pdb2, pdb2_name, nMSA, nENS)






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
                blind_screening(pdb1_name, blind_pred_path)
            else:
                #running_foldseek_all(pdb1_name)

                #### performing the PCA calculation with RMSD
                blind_screening(pdb1_name, blind_pred_path)



        else:
            prediction_all_blind(pdb1_name, search_dir, nMSA, model_type)
            print("               ")
            print("Finished running for prediction using full- and shallow random-MSAs")
            
            print("               ")
            print("Running Foldseek to find the relatedcrystal structures")
            #running_foldseek_all(pdb1_name)

            #### performing the PCA calculation with RMSD
            blind_screening(pdb1_name, blind_pred_path)





    else:
        print("Please type correct option")
