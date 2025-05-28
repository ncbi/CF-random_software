#!/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 13:00:00 2024
 
@author: Myeongsang (Samuel) Lee
"""
import os
import sys
import textalloc as ta
from pathlib import Path
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from numpy import genfromtxt
from matplotlib import pyplot as plt
from adjustText import adjust_text
import glob


ref_Mc = sorted(glob.glob("ref_data_Mchaourab*csv"))
CF_add = sorted(glob.glob("TMScore_random*csv"))
CF_full = sorted(glob.glob("TMScore_full*csv"))



#### printing name
ref_name = []
for name in ref_Mc:
    ref_Mc_name = name.replace('.csv', '')
    ref_Mc_name = ref_Mc_name.replace('ref_data_Mchaourab_', '')

    ref_name = np.append(ref_name, ref_Mc_name)



max_ref_Mc_f1 = np.zeros((np.size(ref_Mc) ,1))
max_ref_Mc_f2 = np.zeros((np.size(ref_Mc) ,1))

max_CF_f1 = np.zeros((np.size(ref_Mc) ,1))
max_CF_f2 = np.zeros((np.size(ref_Mc) ,1))


i = 0
for ref in ref_Mc:
    TMs = genfromtxt(ref, delimiter = ' ')

    max_ref_Mc_f1[i] = np.max(TMs[:, 0])
    max_ref_Mc_f2[i] = np.max(TMs[:, 1])
    i += 1



i = 0 
for CFa, CFf in zip(CF_add, CF_full):
    TMs_CFa = genfromtxt(CFa, delimiter = ' ')
    TMs_CFf = genfromtxt(CFf, delimiter = ' ')



    model1_addi_pdb1 = []; model1_addi_pdb2 = []
    print(TMs_CFa)
    for ii in range(0, TMs_CFa.shape[0]):
        if ii % 2 != 0:
            model1_addi_pdb1 = np.append(model1_addi_pdb1, TMs_CFa[ii])
        else:
            model1_addi_pdb2 = np.append(model1_addi_pdb2, TMs_CFa[ii])

    print(model1_addi_pdb1)
    print(model1_addi_pdb2)



    if np.max(TMs_CFf[0, :]) > np.max(model1_addi_pdb2):
        max_CF_f1[i] = np.max(TMs_CFf[0, :])
    else:
        max_CF_f1[i] = np.max(model1_addi_pdb2)


    if np.max(TMs_CFf[1, :]) > np.max(model1_addi_pdb1):
        max_CF_f2[i] = np.max(TMs_CFf[1, :])
    else:
        max_CF_f2[i] = np.max(model1_addi_pdb1)

    i += 1


bar_empty = np.zeros((np.size(ref_Mc) ,1))
#bar_enpty = np.nan
bar_empty[:] = np.nan


#max_all = np.concatenate((max_ref_Mc_f1, max_ref_Mc_f2, max_CF_f1, max_CF_f2), axis=1)
#max_all = np.concatenate((max_ref_Mc_f1, max_CF_f1, bar_empty, max_ref_Mc_f2, max_CF_f2), axis = 1)
max_all_f1 = np.concatenate((max_ref_Mc_f1, max_CF_f1), axis=1)
max_all_f2 = np.concatenate((max_ref_Mc_f2, max_CF_f2), axis=1)


# plotting the heatmap
y_values = np.arange(0.5, 14.5, 1)
x_values = np.arange(0.5, 2.5, 1)
x_label = ['SPEACH-AF', 'CF-random']
plt.figure(figsize=(7,5))
plt.rcParams["font.family"] = "Helvetica"

f,(ax1,ax2,axcb) = plt.subplots(1,3, 
            gridspec_kw={'width_ratios':[1,1, 0.08]})

g1 = sns.heatmap(data=max_all_f1, annot=True, vmin=0.4, vmax=1, cbar=False, cmap="rocket_r", ax = ax1, fmt=".2f")
g1.set_yticks(y_values, ref_name, fontsize = 13)
g1.set_yticklabels(ref_name, rotation = 0)
g1.set_xticks([])
g1.set_xlabel("Fold1", fontsize = 13)
g1.xaxis.set_label_position('top')
g1.set_xticks(x_values, x_label, fontsize = 13)
g1.xaxis.set_tick_params(bottom = False)
g2 = sns.heatmap(data=max_all_f2, annot=True, vmin=0.4, vmax=1, cmap="rocket_r", ax = ax2, cbar_ax=axcb, fmt=".2f")
g2.set_ylabel('')
g2.set_xticks([])
g2.set_yticks([])
g2.set_xlabel("Fold2", fontsize = 13)
g2.xaxis.set_label_position('top')
g2.figure.axes[-1].set_ylabel('TM-score', size=15, rotation = 270, labelpad = 15)
g2.set_xticks(x_values, x_label, fontsize = 13)
g2.xaxis.set_tick_params(bottom = False)
#g2.set_ylabel("Fold2", fontsize = 15, rotation = 270, labelpad = 15)
#g2.yaxis.set_label_position('right')




plt.tight_layout()
plt.savefig('heatmap-max TMscore comparison.png',dpi=600, transparent=True)
plt.savefig('heatmap-max TMscore comparison.svg',dpi=600, transparent=True)
