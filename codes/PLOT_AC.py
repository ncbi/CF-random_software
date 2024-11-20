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


class plot_2D_scatter_AC():
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


        print("checking plddt")
        print(plddt_addition)
        print(plddt_full)


        TMs_full_resh = np.reshape(TMs_full, ((((nMSA + 5) * 2), 5)))

        f1 = np.concatenate((TMs_addition[0:(nENS + 20), :], TMs_full_resh[0:(nMSA + 5), :]), axis=0)
        print(f1)
        f2 = np.concatenate((TMs_addition[(nENS + 20):(nENS + 20) * 2, :], TMs_full_resh[(nMSA + 5):(nMSA + 5) * 2, :]), axis=0)
        print(f2)


        if np.all(f1 > f2) or np.all(f1 < f2):
            print("Prediction is biased")   
            sys.exit()
        else:
            print("Prediction is not biased")


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

