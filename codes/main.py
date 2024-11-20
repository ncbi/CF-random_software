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
    parser.add_argument("--pname", type=str, help='job name for predicting blind mode' )
    parser.add_argument("--nMSA", type=str, help='number of samples for predicting the structure with MSA')
    parser.add_argument("--nENS", type=str, help='number of samepls for ensemble generation')
    parser.add_argument("--option", type=str, help='select prediction mode inAC, AC and FS e.g. AC = alterantive conformation or FS = fold-switching or inAC = increased sampling for predicting alternative conformation')
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
    elif args.pdb1 is not None and args.pdb2 is not None:
        pdb1 = args.pdb1; pdb2 = args.pdb2
        pdb1_name = pdb1.replace('.pdb','');  pdb2_name = pdb2.replace('.pdb','')


    #if int(args.nMSA) == 0 and int(args.nENS) == 0:
    if args.nMSA is None and args.nENS is None: 
        nMSA = 0; nENS = 0;
    #elif int(args.nMSA) > 0 and int(args.nENS) > 0:
    elif args.nMSA is not None and args.nENS is not None: 
        nMSA = int(args.nMSA); nENS = int(args.nENS)
    elif args.nMSA is None and args.nENS is not None: 
        nMSA = 0; nENS = int(args.nENS)
    elif args.nMSA is not None and args.nENS is None:
        nMSA = int(args.nMSA); nENS = 0
    else:
        print("Please put correct option of nMSA or nENS")
        exit()




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
            pred_1st_all = prediction_all_AC(pdb1, pdb1_name, pdb2, pdb2_name, search_dir, nMSA)
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
                additional_prediction_AC(search_dir, pdb1, pdb1_name, pdb2, pdb2_name, pwd, shallow_MSA_size, nENS)

        else:
            additional_prediction_AC(search_dir, pdb1, pdb1_name, pdb2, pdb2_name, pwd, shallow_MSA_size, nENS)

        ######################################################################################################
        ##### calculate plddt including initial and additional predictions
        if os.path.exists( pwd + add_dir ):
            list_org_samplings = glob.glob( str(pwd) + str(success) + '/' + str(pdb1_name) + '/*full_rand*/')
            list_ran_samplings = glob.glob( str(pwd) + str(success) + '/' + str(pdb1_name) + '/*max*/')
            list_add_samplings = glob.glob( str(pwd) + str(add_dir) + str(pdb1_name))

            full = 'full-MSA'
            random = 'random-MSA'
            addition = 'additional-MSA'
            plddt_cal(list_org_samplings, full, pdb1_name, nMSA, nENS)
            plddt_cal(list_ran_samplings, random, pdb1_name, nMSA, nENS)
            plddt_cal(list_add_samplings, addition, pdb1_name, nMSA, nENS)

        ######################################################################################################
        ##### plot the 2D-scatter plot of TM-scores with pLDDT
        plot_2D_scatter_AC(full, random, addition, pdb1, pdb1_name, pdb2, pdb2_name, nMSA, nENS)


        ######################################################################################################
        ##### merging the plDDT and TM-scores.
        TMs_whole_full = 'TMScore_full-MSA_' + pdb1_name + '.csv'
        TMs_whole_addi = 'TMScore_additional-MSA_' + pdb1_name + '.csv'
        TMs_whole_rand = 'TMScore_random-MSA_' + pdb1_name + '.csv'
        plDDT_full = 'plddt_full-MSA_' + pdb1_name + '.csv'
        plDDT_rand = 'plddt_random-MSA_' + pdb1_name + '.csv'
        plDDT_addi = 'plddt_additional-MSA_' + pdb1_name + '.csv'

        ## deep MSA
        f = open(TMs_whole_full); text = f.read(); f.close()
        f = open(TMs_whole_full, 'w'); f.write("## TM-score of whole structure from deep MSA\n")
        f.write(text); f.close()

        f = open(plDDT_full); text = f.read(); f.close()
        f = open(plDDT_full, 'w'); f.write("## plDDT of whole structure from deep MSA\n")
        f.write(text); f.close()

        # merge seprated files as one file
        deep_output = 'cat ' + TMs_whole_full + ' ' + plDDT_full + ' > TMs_plDDT_full_all_' + pdb1_name + '.csv'
        os.system(deep_output)
        # delete the separated files
        rm_TMs_whole_full = 'mv ' + TMs_whole_full + ' ' + TMs_whole_full + '_backup'; rm_plDDT_full = 'mv ' + plDDT_full + ' ' + plDDT_full + '_backup'
        os.system(rm_TMs_whole_full); os.system(rm_plDDT_full)


        ## random MSAs
        f = open(TMs_whole_rand); text = f.read(); f.close()
        f = open(TMs_whole_rand, 'w'); f.write("## TM-score of whole structure from random MSAs\n")
        f.write(text); f.close()

        f = open(plDDT_rand); text = f.read(); f.close()
        f = open(plDDT_rand, 'w'); f.write("## plDDT of whole structure from random MSAs\n")
        f.write(text); f.close()

        # merge seprated files as one file
        deep_output = 'cat ' + TMs_whole_rand + ' ' + plDDT_rand + ' > TMs_plDDT_rand_all_' + pdb1_name + '.csv'
        os.system(deep_output)
        # delete the separated files
        rm_TMs_whole_rand = 'mv ' + TMs_whole_rand +  ' ' + TMs_whole_rand + '_backup'; rm_plDDT_rand = 'mv ' + plDDT_rand + ' ' + plDDT_rand + '_backup'
        os.system(rm_TMs_whole_rand); os.system(rm_plDDT_rand)


        ## ensemble generation
        f = open(TMs_whole_addi); text = f.read(); f.close()
        f = open(TMs_whole_addi, 'w'); f.write("## TM-score of whole structure from ensemble generation\n")
        f.write(text); f.close()

        f = open(plDDT_addi); text = f.read(); f.close()
        f = open(plDDT_addi, 'w'); f.write("## plDDT of whole structure from deep generation\n")
        f.write(text); f.close()

        # merge seprated files as one file
        deep_output = 'cat ' + TMs_whole_addi + ' ' + plDDT_addi + ' > TMs_plDDT_addi_all_' + pdb1_name + '.csv'
        os.system(deep_output)
        # delete the separated files
        rm_TMs_whole_addi = 'mv ' + TMs_whole_addi + ' ' + TMs_whole_addi + '_backup'; rm_plDDT_addi = 'mv ' + plDDT_addi + ' ' + plDDT_addi + '_backup'
        os.system(rm_TMs_whole_addi); os.system(rm_plDDT_addi)








    elif args.option == "FS":
        print("Predicting fold-swithcing models")
        ######################################################################################################
        ###### running prediction using full- and shallow random-MSA
        if os.path.exists(success + '/' + pdb1_name) and succ_dir_count >= 8:
            print("Predictions including full- and random-MSA were already done")
        else:
            pred_1st_all = prediction_all(pdb1, pdb1_name, pdb2, pdb2_name, search_dir, nMSA)
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
                additional_prediction(search_dir, pdb1, pdb1_name, pdb2, pdb2_name, pwd, shallow_MSA_size, nENS)
        else:
            additional_prediction(search_dir, pdb1, pdb1_name, pdb2, pdb2_name, pwd, shallow_MSA_size, nENS)
    
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
            plddt_cal(list_org_samplings, full, pdb1_name, nMSA, nENS)
            plddt_cal(list_ran_samplings, random, pdb1_name, nMSA, nENS)
            plddt_cal(list_add_samplings, addition, pdb1_name, nMSA, nENS)
    
    
        ######################################################################################################
        ##### plot the 2D-scatter plot of TM-scores with pLDDT
        plot_2D_scatter(full, random, addition, pdb1, pdb1_name, pdb2, pdb2_name, nMSA, nENS)


        ######################################################################################################
        ##### merging the plDDT and TM-scores.
        TMs_whole_full = 'TMScore_full-MSA_' + pdb1_name + '.csv'
        TMs_whole_addi = 'TMScore_additional-MSA_' + pdb1_name + '.csv'
        TMs_whole_rand = 'TMScore_random-MSA_' + pdb1_name + '.csv'
        TMs_fs_full = 'TMScore_fs_full-MSA_' + pdb1_name + '.csv'
        TMs_fs_addi = 'TMScore_fs_additional-MSA_' + pdb1_name + '.csv'
        TMs_fs_rand = 'TMScore_fs_random-MSA_' + pdb1_name + '.csv'
        plDDT_full = 'plddt_full-MSA_' + pdb1_name + '.csv'
        plDDT_rand = 'plddt_random-MSA_' + pdb1_name + '.csv'
        plDDT_addi = 'plddt_additional-MSA_' + pdb1_name + '.csv'

        ## deep MSA
        f = open(TMs_whole_full); text = f.read(); f.close()
        f = open(TMs_whole_full, 'w'); f.write("## TM-score of whole structure from deep MSA\n")
        f.write(text); f.close()

        f = open(TMs_fs_full); text = f.read(); f.close()
        f = open(TMs_fs_full, 'w'); f.write("## TM-score of fold-switching region from deep MSA\n")
        f.write(text); f.close()

        f = open(plDDT_full); text = f.read(); f.close()
        f = open(plDDT_full, 'w'); f.write("## plDDT of whole structure from deep MSA\n")
        f.write(text); f.close()

        # merge seprated files as one file
        deep_output = 'cat ' + TMs_whole_full + ' ' + TMs_fs_full + ' ' + plDDT_full + ' > TMs_plDDT_full_all_' + pdb1_name + '.csv' 
        os.system(deep_output)
        # delete the separated files
        rm_TMs_whole_full = 'mv ' + TMs_whole_full + ' ' + TMs_whole_full + '_backup'; 
        rm_TMs_fs_full = 'mv ' + TMs_fs_full + ' ' + TMs_fs_full + '_backup'; 
        rm_plDDT_full = 'mv ' + plDDT_full + ' ' + plDDT_full + '_backup'
        os.system(rm_TMs_whole_full); os.system(rm_TMs_fs_full); os.system(rm_plDDT_full)


        ## random MSAs 
        f = open(TMs_whole_rand); text = f.read(); f.close()
        f = open(TMs_whole_rand, 'w'); f.write("## TM-score of whole structure from random MSAs\n")
        f.write(text); f.close()

        f = open(TMs_fs_rand); text = f.read(); f.close()
        f = open(TMs_fs_rand, 'w'); f.write("## TM-score of fold-switching region from random MSAs\n")
        f.write(text); f.close()

        f = open(plDDT_rand); text = f.read(); f.close()
        f = open(plDDT_rand, 'w'); f.write("## plDDT of whole structure from random MSAs\n")
        f.write(text); f.close()
        
        # merge seprated files as one file
        deep_output = 'cat ' + TMs_whole_rand + ' ' + TMs_fs_rand + ' ' + plDDT_rand + ' > TMs_plDDT_rand_all_' + pdb1_name + '.csv'
        os.system(deep_output)
        # delete the separated files
        rm_TMs_whole_rand = 'mv ' + TMs_whole_rand + ' ' + TMs_whole_rand + '_backup'; 
        rm_TMs_fs_rand = 'mv ' + TMs_fs_rand + ' ' + TMs_fs_rand + '_backup'; 
        rm_plDDT_rand = 'mv ' + plDDT_rand + ' ' + plDDT_rand + '_backup'
        os.system(rm_TMs_whole_rand); os.system(rm_TMs_fs_rand); os.system(rm_plDDT_rand)

        
        ## ensemble generation
        f = open(TMs_whole_addi); text = f.read(); f.close()
        f = open(TMs_whole_addi, 'w'); f.write("## TM-score of whole structure from ensemble generation\n")
        f.write(text); f.close()

        f = open(TMs_fs_addi); text = f.read(); f.close()
        f = open(TMs_fs_addi, 'w'); f.write("## TM-score of fold-switching region from enselble generation\n")
        f.write(text); f.close()

        f = open(plDDT_addi); text = f.read(); f.close()
        f = open(plDDT_addi, 'w'); f.write("## plDDT of whole structure from deep generation\n")
        f.write(text); f.close()

        # merge seprated files as one file
        deep_output = 'cat ' + TMs_whole_addi + ' ' + TMs_fs_addi + ' ' + plDDT_addi + ' > TMs_plDDT_addi_all_' + pdb1_name + '.csv'
        os.system(deep_output)
        # delete the separated files
        rm_TMs_whole_addi = 'mv ' + TMs_whole_addi + ' ' + TMs_whole_addi + '_backup'; 
        rm_TMs_fs_addi = 'mv ' + TMs_fs_addi + ' ' + TMs_fs_addi + '_backup'; 
        rm_plDDT_addi = 'mv ' + plDDT_addi + ' ' + plDDT_addi + '_backup'
        os.system(rm_TMs_whole_addi); os.system(rm_TMs_fs_addi); os.system(rm_plDDT_addi)






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
            prediction_all_blind(pdb1_name, search_dir, nMSA)
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
