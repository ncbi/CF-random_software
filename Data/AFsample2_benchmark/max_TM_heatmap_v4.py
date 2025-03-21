#!/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 13:00:00 2024
 
@author: Myeongsang (Samuel) Lee
"""
import os
import sys
import csv
import textalloc as ta
from pathlib import Path
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.font_manager
from numpy import genfromtxt
from matplotlib import pyplot as plt
from adjustText import adjust_text
import glob

CF_add = sorted(glob.glob("TMScore_random*csv"))
CF_full = sorted(glob.glob("TMScore_full*csv"))


data_name = open('list_of_OC23-uniprot_ID.csv', 'r')
myreader = csv.reader(data_name)

uniprot_name = []; pdb1_name = []; pdb2_name = []; ref_TM_open = []; ref_TM_close = []

for row in myreader:
    uniprot_name = np.append(uniprot_name, row[0])
    pdb1_name = np.append(pdb1_name, row[1])
    pdb2_name = np.append(pdb2_name, row[2])
    ref_TM_open = np.append(ref_TM_open, row[5])
    ref_TM_close = np.append(ref_TM_close, row[6])

uniprot_name = uniprot_name[1:]
pdb1_name = pdb1_name[1:]
pdb2_name = pdb2_name[1:]
uniprot_size = np.size(uniprot_name)


uniprot_name = uniprot_name.reshape(uniprot_size, 1)
pdb1_name = pdb1_name.reshape(uniprot_size, 1); pdb2_name = pdb2_name.reshape(uniprot_size, 1)
data_all = np.concatenate((uniprot_name, pdb1_name, pdb2_name), axis = 1)

pdb1_name = pdb1_name.ravel().astype(str)
pdb2_name = pdb2_name.ravel().astype(str)
uniprot_name = uniprot_name.ravel().astype(str)

max_CF_f1 = np.zeros((uniprot_size, 1)); max_CF_f2 = np.zeros((uniprot_size, 1))
model1_f1 = np.zeros((uniprot_size, 1)); model1_f2 = np.zeros((uniprot_size, 1))
model2_f1 = np.zeros((uniprot_size, 1)); model2_f2 = np.zeros((uniprot_size, 1))



i = 0
for model1, model2 in zip(pdb1_name, pdb2_name):
    print(model1)
    tmp1_add = "TMScore_random-MSA_" + str(model1) + ".csv" 
    tmp1_ful = "TMScore_full-MSA_" + str(model1) + ".csv"

    model1_addi = genfromtxt(tmp1_add, delimiter = ' ')
    model1_full = genfromtxt(tmp1_ful, delimiter = ' ')

    ##### Some cases, the number of samples more than 300 mimicking the reference data
    #if np.size(TMs_CFa) < 50:
    #    if np.max(TMs_CFf[0, :]) > np.max(TMs_CFa[0:30, :]):
    #        max_CF_f1[i] = np.max(TMs_CFf[0, :])
    #    else:
    #        max_CF_f1[i] = np.max(TMs_CFa[0:30, :])

    #    if np.max(TMs_CFf[1, :]) > np.max(TMs_CFa[30:60, :]):
    #        max_CF_f2[i] = np.max(TMs_CFf[1, :])
    #    else:
    #        max_CF_f2[i] = np.max(TMs_CFa[30:60, :])

    ##### model1 part   
    
    model1_addi_pdb1 = []; model1_addi_pdb2 = []
    for ii in range(0, model1_addi.shape[0]):
        if ii % 2 != 0:
            model1_addi_pdb1 = np.append(model1_addi_pdb1, model1_addi[ii])
        else:
            model1_addi_pdb2 = np.append(model1_addi_pdb2, model1_addi[ii])


    if np.max(model1_full[0, :]) > np.max(model1_addi_pdb2):
        model1_f1[i] = np.max(model1_full[0, :])
    else:
        model1_f1[i] = np.max(model1_addi_pdb2)

    if np.max(model1_full[1, :]) > np.max(model1_addi_pdb1):
        model1_f2[i] = np.max(model1_full[1, :])
    else:
        model1_f2[i] = np.max(model1_addi_pdb1)


    max_CF_f1[i] = model1_f1[i]; max_CF_f2[i] = model1_f2[i]

    print(max_CF_f1[i], max_CF_f2[i])
    i += 1





#max_all = np.concatenate((max_CF_f1, max_CF_f2), axis=1)
#print(max_all)
ref_TM_open = ref_TM_open[1:]; ref_TM_close = ref_TM_close[1:]
ref_TM_open = ref_TM_open.ravel().astype(float)
ref_TM_close = ref_TM_close.ravel().astype(float)
ref_TM_open = np.array(ref_TM_open); ref_TM_close = np.array(ref_TM_close)
ref_TM_open = ref_TM_open.reshape(uniprot_size, 1); ref_TM_close = ref_TM_close.reshape(uniprot_size, 1)


np.set_printoptions(precision=2)

max_all_f1 = np.concatenate((ref_TM_open, max_CF_f1), axis=1)
max_all_f2 = np.concatenate((ref_TM_close, max_CF_f2), axis=1)



# plotting the heatmap
#y_values = np.arange(0.5, int(np.size((uniprot_name))) + 0.5, 1)
#x_values = np.arange(0.5, 2.5, 1)
#x_label = ['Fold1', 'Fold2']
#
#
#plt.figure(figsize=(7,5))
#plt.rcParams["font.family"] = "Helvetica"
#hm = sns.heatmap(data=max_all,
#                       annot=True, vmin=0.5, vmax=1, cmap="rocket_r", cbar_kws={'label': 'TM-score'})
#hm.figure.axes[-1].set_ylabel('TM-score', size=15, rotation = 270, labelpad = 15)
#hm.figure.axes[-1].tick_params(labelsize=15)
#
#
#x_label = ['Fold1', 'Fold2']
#plt.yticks(y_values, uniprot_name, fontsize= 13)
#plt.yticks(rotation=0)
#plt.xticks(x_values, x_label, fontsize = 15)
#plt.tick_params(
#    axis='x',          # changes apply to the x-axis
#    bottom=False,      # ticks along the bottom edge are off
#    top=False)         # ticks along the top edge are off

# plotting the heatmap
y_values = np.arange(0.5, 23.5, 1)
x_values = np.arange(0.5, 2.5, 1)
x_label = ['AFSample2', 'CF-random']
plt.figure(figsize=(7,5))
plt.rcParams["font.family"] = "Helvetica"

f,(ax1,ax2,axcb) = plt.subplots(1,3,
            gridspec_kw={'width_ratios':[1,1, 0.08]})

g1 = sns.heatmap(data=max_all_f1, annot=True, vmin=0.5, vmax=1.0, cbar=False, cmap="rocket_r", ax = ax1, fmt=".2f")
g1.set_yticks(y_values, uniprot_name, fontsize = 13)
g1.set_yticklabels(uniprot_name, rotation = 0)
g1.set_xticks([])
g1.set_xlabel("Open", fontsize = 13)
g1.xaxis.set_label_position('top')
g1.set_xticks(x_values, x_label, fontsize = 13)
g1.xaxis.set_tick_params(bottom = False)
g2 = sns.heatmap(data=max_all_f2, annot=True, vmin=0.5, vmax=1.0, cmap="rocket_r", ax = ax2, cbar_ax=axcb, fmt=".2f")
g2.set_ylabel('')
g2.set_xticks([])
g2.set_yticks([])
g2.set_xlabel("Close", fontsize = 13)
g2.xaxis.set_label_position('top')
g2.figure.axes[-1].set_ylabel('TM-score', size=15, rotation = 270, labelpad = 15)
g2.set_xticks(x_values, x_label, fontsize = 13)
g2.xaxis.set_tick_params(bottom = False)




plt.tight_layout()
plt.savefig('heatmap-max TMscore comparison.png',dpi=600, transparent=True)

