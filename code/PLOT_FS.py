#!/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 22 13:40:00 2024
 
@author: Myeongsang (Samuel) Lee
"""
import os
import sys
import textalloc as ta
import seaborn as sns
from pathlib import Path
import numpy as np
from numpy import genfromtxt
from matplotlib import pyplot as plt
from adjustText import adjust_text
import glob

from cal_tmscore_fs_flmsa import *
from fs_seq_compare import *

class plot_2D_scatter():
    def __init__(self, full_cate, random_cate, pdb1, pdb1_name, pdb2, pdb2_name, nMSA, nENS):
        ##### load TM-scores both full- and ramdon-MSA
        TMs_full =  genfromtxt("TMScore_" + full_cate + "_" + pdb1_name + ".csv", delimiter = ' ' )
        TMs_random =  genfromtxt("TMScore_" + random_cate + "_" + pdb1_name + ".csv", delimiter = ' ' )

        ############ load pLDDT scores both full- and ramdon-MSA
        plddt_full = genfromtxt("plddt_" + full_cate + "_" + pdb1_name + ".csv", delimiter = ' ' )
        plddt_random = genfromtxt("plddt_" + random_cate + "_" + pdb1_name + ".csv", delimiter = ' ' )


        #################################################################
        ########### getting the TM-score values of fold-switching region

        pwd = os.getcwd() + '/'

        fs_full_TMs = genfromtxt("TMScore_fs_" + full_cate + "_" + pdb1_name + ".csv", delimiter = ' ')
        TMs_fs_full = fs_full_TMs
        fs_random_TMs = genfromtxt("TMScore_fs_" + random_cate + "_" + pdb1_name + ".csv", delimiter = ' ')
        TMs_fs_random = fs_random_TMs



        ######### plotting the TM-score values as 2D scatter plot
        print("                                        ")
        print("Size of column: ", TMs_random.shape[-1])
        print("Size of row: ", TMs_random.shape[0])
        print("Dimension: ", TMs_random.ndim)

        print("                                        ")
        print(TMs_random)
        print("                                        ")
        print(TMs_full)


        print("checking plddt")
        print(plddt_full)
        print(plddt_random)

        plddt_random = np.reshape(plddt_random, (7,  (nMSA + 5) * 5))
        TMs_fs_full_resh = np.reshape(TMs_fs_full, ((((nMSA + 5) * 2), 5)))




        plt.figure(0)


        for ii in range(0, int(TMs_random.shape[0] / 2) ):
            plt.scatter(TMs_random[ii * 2, :], TMs_random[(ii * 2 + 1), :], c = plddt_random[ii, :], cmap='rocket_r',  vmin=50, vmax=100, s=35, marker="o")
            

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
        for ii in range(0, int(TMs_random.shape[0] / 2) ):
            plt.scatter(TMs_fs_random[ii * 2, :], TMs_fs_random[(ii * 2 + 1), :], c = plddt_random[ii, :], cmap='plasma', vmin=50, vmax=100, s=35, marker="o")


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

