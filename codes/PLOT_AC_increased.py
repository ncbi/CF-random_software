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
    def __init__(self, full_cate, random_cate, addition_cate, pdb1, pdb1_name, pdb2, pdb2_name):
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
        #plt.figure(figsize=(11.5, 12.5))


        print("checking plddt")
        print(plddt_addition)
        print(plddt_full)
        



        TMs_full_resh = np.reshape(TMs_full, (50, 5))

        f1 = np.concatenate((TMs_addition[0:30, :], TMs_full_resh[0:25, :]), axis=0)
        print(f1)
        f2 = np.concatenate((TMs_addition[30:60, :], TMs_full_resh[25:50, :]), axis=0)
        print(f2)


        if np.all(f1 > f2) or np.all(f1 < f2):
            print("Prediction is biased")   
            #sys.exit()
            print("Plotting the TM-score result just in case")
        else:
            print("Prediction is not biased")


        plt.figure(0)

        for ii in range(0, int((TMs_addition.shape[0] / 2) - 1)):
            plt.scatter(TMs_addition[ii, :], TMs_addition[(ii + 20), :], c = plddt_addition[ii, :], cmap='plasma', vmin=50, vmax=100, s=35, marker="o")
            

        clb=plt.colorbar()
        clb.ax.tick_params(labelsize=15)


        #plt.scatter(TMs_full[0, :], TMs_full[1, :], c = 'green', vmin=60, vmax=100, s=35, marker="o")
        plt.scatter(TMs_full[0, :], TMs_full[1, :], c = plddt_full, cmap='plasma', vmin=50, vmax=100, s=35, marker="o")


        x = [ 0 , 1 ]
        y = [ 0 , 1 ]

        plt.ylim(0, 1)
        plt.xlim(0, 1)
    

        #plt.xticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
        #plt.yticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
        #plt.colorbar()
        
        
        #plt.scatter(tmp1_0, tmp2_0, color = 'red', s=35, marker="o")
        plt.plot(x, y, linestyle='dashed', color = 'black')
        
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        
        plt.xlabel('TM-Score similar to fold1(' + pdb1_name + ')', fontsize=15); plt.ylabel('TM-score similar to fold2(' + pdb2_name + ')', fontsize=15)
        plt.savefig('TMscore_' + full_cate  + '_' + pdb1_name +  '.png', transparent = True)







        ################################################################
        ######## additional plot for labeling - checking purpose #######
        ################################################################
        #plt.figure(3)
        #plt.figure(figsize=(30, 30), dpi=300)

        #lbl_add = [] ; lbl_full = [];
        #for lll in range(0, int(((TMs_fs_addition.shape[0]) * 5) / 2)):
        #    #lbl_add.append("Ra_" + str(lll))
        #    lbl_add.append(str(lll))
        #print(lbl_add)




        #TMs_fs_addition_re = []
        #plddt_addition_re = []
        #print(TMs_fs_addition)
        #for ii in range(0, int(((TMs_fs_addition.shape[0])))):
        #    for iii in range(0, 5):
        #        TMs_fs_addition_re.append(TMs_fs_addition[ii, iii])

        #


        #for ii in range(0, int(((TMs_fs_addition.shape[0]) / 2))):
        #    for iii in range(0, 5):
        #        plddt_addition_re.append(plddt_addition[ii, iii])                


        #TMs_fs_addition_resh = np.reshape(TMs_fs_addition_re, (200, 1))
        #TMs_fs_addition_resh_1 = np.zeros((100, 1))
        #TMs_fs_addition_resh_2 = np.zeros((100, 1))
        #TMs_fs_addition_resh_1 = TMs_fs_addition_resh[0:100, 0] 
        #TMs_fs_addition_resh_2 = TMs_fs_addition_resh[100:200, 0] 
        


        #print("Adding the label")

        #x = TMs_fs_addition_resh_1; y = TMs_fs_addition_resh_2
        #fig, ax = plt.subplots(ncols=2, nrows=2, figsize=(15, 15), dpi=300, layout="constrained")


        #xm = [0 ,1] ; ym = [0, 1]


        #ax[0, 0].scatter(x[0:25], y[0:25], c = 'b')
        #text_list = [f'R{i}' for i in range(0, 25)]
        #ta.allocate(ax[0, 0], x[0:25], y[0:25],
        #    text_list,
        #    x_scatter=x[0:25], y_scatter=y[0:25],
        #    textsize=10)
        #ax[0, 0].plot(xm, ym, linestyle='dashed', color = 'black')
        #ax[0, 0].set_xlim(np.min(x), np.max(x)); ax[0, 0].set_ylim(np.min(y), np.max(y))



        #ax[0, 1].scatter(x[25:50], y[25:50], c = 'b')
        #text_list = [f'R{i}' for i in range(25, 50)]
        #ta.allocate(ax[0, 1], x[25:50], y[25:50],
        #    text_list,
        #    x_scatter=x[25:50], y_scatter=y[25:50],
        #    textsize=10)
        #ax[0, 1].plot(xm, ym, linestyle='dashed', color = 'black')
        #ax[0, 1].set_xlim(np.min(x), np.max(x)); ax[0, 1].set_ylim(np.min(y), np.max(y))        



        #ax[1, 0].scatter(x[50:75], y[50:75], c = 'b')
        #text_list = [f'R{i}' for i in range(50, 75)]
        #ta.allocate(ax[1, 0], x[50:75], y[50:75],
        #    text_list,
        #    x_scatter=x[50:75], y_scatter=y[50:75],
        #    textsize=10)
        #ax[1, 0].plot(xm, ym, linestyle='dashed', color = 'black')
        #ax[1, 0].set_xlim(np.min(x), np.max(x)); ax[1, 0].set_ylim(np.min(y), np.max(y))



        #ax[1, 1].scatter(x[75:100], y[75:100], c = 'b')
        #text_list = [f'R{i}' for i in range(75, 100)]
        #ta.allocate(ax[1, 1], x[75:100], y[75:100],
        #    text_list,
        #    x_scatter=x[75:100], y_scatter=y[75:100],
        #    textsize=10)
        #ax[1, 1].plot(xm, ym, linestyle='dashed', color = 'black')
        #ax[1, 1].set_xlim(np.min(x), np.max(x)); ax[1, 1].set_ylim(np.min(y), np.max(y))
        

