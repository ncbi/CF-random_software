#!/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 22 13:40:00 2024
 
@author: Myeongsang (Samuel) Lee
"""
import os
import sys
import textalloc as ta
from pathlib import Path
import numpy as np
from numpy import genfromtxt
from matplotlib import pyplot as plt
from adjustText import adjust_text
import glob

from cal_tmscore_fs_flmsa import *
from fs_seq_compare import *

class plot_2D_scatter():
    def __init__(self, full_cate, random_cate, addition_cate, pdb1, pdb1_name, pdb2, pdb2_name, nMSA, nENS):
        ##### load TM-scores both full- and ramdon-MSA
        TMs_full =  genfromtxt("TMScore_" + full_cate + "_" + pdb1_name + ".csv", delimiter = ' ' )
        TMs_random =  genfromtxt("TMScore_" + random_cate + "_" + pdb1_name + ".csv", delimiter = ' ' )
        TMs_addition = genfromtxt("TMScore_" + addition_cate + "_" + pdb1_name + ".csv", delimiter = ' ' )

        ############ load pLDDT scores both full- and ramdon-MSA
        plddt_full = genfromtxt("plddt_" + full_cate + "_" + pdb1_name + ".csv", delimiter = ' ' )
        plddt_random = genfromtxt("plddt_" + random_cate + "_" + pdb1_name + ".csv", delimiter = ' ' )
        plddt_addition = genfromtxt("plddt_" + addition_cate + "_" + pdb1_name + ".csv", delimiter = ' ' )


        #################################################################
        ########### getting the TM-score values of fold-switching region


        pwd = os.getcwd() + '/'

        pred_path = 'additional_sampling/rand_*/'
        fs_addition_TMs =  genfromtxt("TMScore_fs_" + addition_cate + "_" + pdb1_name + ".csv", delimiter = ' ' )
        fs_full_TMs = genfromtxt("TMScore_fs_" + full_cate + "_" + pdb1_name + ".csv", delimiter = ' ')
        TMs_fs_full = fs_full_TMs
        TMs_fs_addition = fs_addition_TMs
        TMs_fs_size = int((TMs_fs_addition.shape[0] / 2))



        ######### plotting the TM-score values as 2D scatter plot
        print("                                        ")
        print("Size of column: ", TMs_random.shape[-1])
        print("Size of row: ", TMs_random.shape[0])
        print("Dimension: ", TMs_random.ndim)

        print("                                        ")
        print(TMs_random)
        print("                                        ")
        print(TMs_full)
        print("                                        ")
        print(TMs_addition)
        #plt.figure(figsize=(11.5, 12.5))


        print("checking plddt")
        print(plddt_addition)
        print(plddt_full)
        



        TMs_fs_full_resh = np.reshape(TMs_fs_full, ((((nMSA + 5) * 2), 5)))

        f1 = np.concatenate((TMs_fs_addition[0:(nENS + 20), :], TMs_fs_full_resh[0:(nMSA + 5) :]), axis=0)
        print(f1)
        f2 = np.concatenate((TMs_fs_addition[(nENS + 20):(nENS + 20) * 2, :], TMs_fs_full_resh[(nMSA + 5):(nMSA + 5) * 2, :]), axis=0)
        print(f2)


        if np.all(f1 > f2) or np.all(f1 < f2):
            print("Prediction is biased")   
            sys.exit()
        else:
            print("Prediction is not biased")
            print("Double-checking secondary structure ratio between crystal structure and predictions")
            if np.all(f1 < 0.5) and np.all(f2 < 0.5):
                pred_dir = 'additional_sampling/' + pdb1_name + '/'
                fs_range(pdb1, pdb2, pdb1_name, pdb2_name, pred_dir)
                pred_dir = 'successed_prediction/' + pdb1_name + '/*full*'
                fs_range(pdb1, pdb2, pdb1_name, pdb2_name, pred_dir)


                #fs_compare_output.log
                comp_read = genfromtxt("fs_compare_output_" + pdb1_name + ".log")

                if comp_read == "success":
                    os.system('rm -rf output*log')
                    os.system('rm -f fs_compare_output.log')
                else:
                    sys.exit()


        plt.figure(0)

        for ii in range(0, int((TMs_addition.shape[0] / 2) - 1)):
            plt.scatter(TMs_addition[ii, :], TMs_addition[(ii + (nENS + 20)), :], c = plddt_addition[ii, :], cmap='plasma', vmin=50, vmax=100, s=35, marker="o")
            

        clb=plt.colorbar()
        clb.ax.tick_params(labelsize=15)


        plt.scatter(TMs_full[0, :], TMs_full[1, :], c = plddt_full, cmap='plasma', vmin=50, vmax=100, s=35, marker="o")


        x = [ 0 , 1 ]
        y = [ 0 , 1 ]

        plt.ylim(0, 1)
        plt.xlim(0, 1)
    

        plt.plot(x, y, linestyle='dashed', color = 'black')
        
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        
        plt.xlabel('TM-Score similar to fold1(' + pdb1_name + ')', fontsize=15); plt.ylabel('TM-score similar to fold2(' + pdb2_name + ')', fontsize=15)
        plt.savefig('TMscore_' + full_cate  + '_' + pdb1_name +  '.png', transparent = True)











        plt.figure(1)
        for ii in range(0, int((TMs_fs_addition.shape[0] / 2) - 1)):
            plt.scatter(TMs_fs_addition[ii, :], TMs_fs_addition[(ii + (nENS + 20)), :], c = plddt_addition[ii, :], cmap='plasma', vmin=50, vmax=100, s=35, marker="o")


        x = [ 0.0 , 1 ] ;  y = [ 0.0 , 1 ]
        plt.ylim(0.0, 1)
        plt.xlim(0.0, 1)

        dlb=plt.colorbar()
        dlb.ax.tick_params(labelsize=15)

        plt.scatter(TMs_fs_full[0, :], TMs_fs_full[1, :], c = plddt_full, cmap='plasma', vmin=50, vmax=100, s=35, marker="o")


        plt.plot(x, y, linestyle='dashed', color = 'black')

        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)

        plt.xlabel('TM-Score similar to fold1(' + pdb1_name + ')', fontsize=15); plt.ylabel('TM-score similar to fold2(' + pdb2_name + ')', fontsize=15)
        plt.savefig('TMscore_fs-region_' + full_cate  + '_' + pdb1_name +  '.png', transparent = True)





        ###############################################################
        ####### additional plot for labeling - checking purpose #######
        ###############################################################
        plt.figure(3)
        plt.figure(figsize=(30, 30), dpi=300)

        lbl_add = [] ; lbl_full = [];
        for lll in range(0, int(((TMs_fs_addition.shape[0]) * 5) / 2)):
            lbl_add.append(str(lll))
        print(lbl_add)



        TMs_fs_addition_re = []
        plddt_addition_re = []
        print(TMs_fs_addition)
        for ii in range(0, int(((TMs_fs_addition.shape[0])))):
            for iii in range(0, 5):
                TMs_fs_addition_re.append(TMs_fs_addition[ii, iii])

        


        for ii in range(0, int(((TMs_fs_addition.shape[0]) / 2))):
            for iii in range(0, 5):
                plddt_addition_re.append(plddt_addition[ii, iii])                



        total_n_data = (nENS + 20) * 5 * 2
        total_n_data_half = (nENS + 20) * 5
        TMs_fs_addition_resh = np.reshape(TMs_fs_addition_re, (total_n_data, 1))
        TMs_fs_addition_resh_1 = np.zeros((int((total_n_data / 2)), 1))
        TMs_fs_addition_resh_2 = np.zeros((int((total_n_data / 2)), 1))
        TMs_fs_addition_resh_1 = TMs_fs_addition_resh[0:int(total_n_data / 2), 0] 
        TMs_fs_addition_resh_2 = TMs_fs_addition_resh[int((total_n_data / 2)):total_n_data, 0] 

        print("Adding the label")
        x = TMs_fs_addition_resh_1; y = TMs_fs_addition_resh_2
        fig, ax = plt.subplots(ncols=2, nrows=2, figsize=(15, 15), dpi=300, layout="constrained")


        xm = [0 ,1] ; ym = [0, 1]


        ax[0, 0].scatter(x[0:int(total_n_data_half / 4)], y[0:int(total_n_data_half / 4)], c = 'b')
        text_list = [f'R{i}' for i in range(0, int(total_n_data_half / 4))]
        ta.allocate(ax[0, 0], x[0:int(total_n_data_half / 4)], y[0:int(total_n_data_half / 4)],
            text_list,
            x_scatter=x[0:int(total_n_data_half / 4)], y_scatter=y[0:int(total_n_data_half / 4)],
            textsize=10)
        ax[0, 0].plot(xm, ym, linestyle='dashed', color = 'black')
        ax[0, 0].set_xlim(np.min(x), np.max(x)); ax[0, 0].set_ylim(np.min(y), np.max(y))



        ax[0, 1].scatter(x[int(total_n_data_half / 4):int(total_n_data_half / 2)], y[int(total_n_data_half / 4):int(total_n_data_half / 2)], c = 'b')
        text_list = [f'R{i}' for i in range(int(total_n_data_half / 4), int(total_n_data_half / 2))]
        ta.allocate(ax[0, 1], x[int(total_n_data_half / 4):int(total_n_data_half / 2)], y[int(total_n_data_half / 4):int(total_n_data_half / 2)],
            text_list,
            x_scatter=x[int(total_n_data_half / 4):int(total_n_data_half / 2)], y_scatter=y[int(total_n_data_half / 4):int(total_n_data_half / 2)],
            textsize=10)
        ax[0, 1].plot(xm, ym, linestyle='dashed', color = 'black')
        ax[0, 1].set_xlim(np.min(x), np.max(x)); ax[0, 1].set_ylim(np.min(y), np.max(y))        



        ax[1, 0].scatter(x[int(total_n_data_half / 2):int(total_n_data_half * 3 / 2)], y[int(total_n_data_half / 2):int(total_n_data_half * 3 / 2)], c = 'b')
        text_list = [f'R{i}' for i in range(int(total_n_data_half / 2), int(total_n_data_half * 3 / 2))]
        ta.allocate(ax[1, 0], x[int(total_n_data_half / 2):int(total_n_data_half * 3 / 2)], y[int(total_n_data_half / 2):int(total_n_data_half * 3 / 2)],
            text_list,
            x_scatter=x[int(total_n_data_half / 2):int(total_n_data_half * 3 / 2)], y_scatter=y[int(total_n_data_half / 2):int(total_n_data_half * 3 / 2)],
            textsize=10)
        ax[1, 0].plot(xm, ym, linestyle='dashed', color = 'black')
        ax[1, 0].set_xlim(np.min(x), np.max(x)); ax[1, 0].set_ylim(np.min(y), np.max(y))



        ax[1, 1].scatter(x[int(total_n_data_half * 3 / 2):int(total_n_data_half)], y[int(total_n_data_half * 3 / 2):int(total_n_data_half)], c = 'b')
        text_list = [f'R{i}' for i in range(int(total_n_data_half * 3 / 2), int(total_n_data_half))]
        ta.allocate(ax[1, 1], x[int(total_n_data_half * 3 / 2):int(total_n_data_half)], y[int(total_n_data_half * 3 / 2):int(total_n_data_half)],
            text_list,
            x_scatter=x[int(total_n_data_half * 3 / 2):int(total_n_data_half)], y_scatter=y[int(total_n_data_half * 3 / 2):int(total_n_data_half)],
            textsize=10)
        ax[1, 1].plot(xm, ym, linestyle='dashed', color = 'black')
        ax[1, 1].set_xlim(np.min(x), np.max(x)); ax[1, 1].set_ylim(np.min(y), np.max(y))
        

        plt.savefig('TMscore_fs-region_' + full_cate  + '_' + pdb1_name +  '_label.png', transparent = True)
