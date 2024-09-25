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
from numpy import genfromtxt
from matplotlib import pyplot as plt
from adjustText import adjust_text
import glob

CF_add = sorted(glob.glob("TMScore_additional*csv"))
CF_full = sorted(glob.glob("TMScore_full*csv"))


data_name = open('list_of_OC23-uniprot_ID.csv', 'r')
myreader = csv.reader(data_name)

uniprot_name = []; pdb1_name = []; pdb2_name = []

for row in myreader:
    uniprot_name = np.append(uniprot_name, row[0])
    pdb1_name = np.append(pdb1_name, row[1])
    pdb2_name = np.append(pdb2_name, row[2])

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
    tmp1_add = "TMScore_additional-MSA_" + str(model1) + ".csv" 
    tmp1_ful = "TMScore_full-MSA_" + str(model1) + ".csv"

    tmp2_add = "TMScore_additional-MSA_" + str(model2) + ".csv"
    tmp2_ful = "TMScore_full-MSA_" + str(model2) + ".csv"

    model1_addi = genfromtxt(tmp1_add, delimiter = ' ')
    model1_full = genfromtxt(tmp1_ful, delimiter = ' ')

    model2_addi = genfromtxt(tmp2_add, delimiter = ' ')
    model2_full = genfromtxt(tmp2_ful, delimiter = ' ')

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
    if np.max(model1_full[0, :]) > np.max(model1_addi[0:20, :]):
        model1_f1[i] = np.max(model1_full[0, :])
    else:
        model1_f1[i] = np.max(model1_addi[0:20, :])

    if np.max(model1_full[1, :]) > np.max(model1_addi[20:40, :]):
        model1_f2[i] = np.max(model1_full[1, :])
    else:
        model1_f2[i] = np.max(model1_addi[20:40, :])


    #### model2 part
    if np.max(model1_full[0, :]) > np.max(model1_addi[0:20, :]):
        model2_f1[i] = np.max(model1_full[0, :])
    else:
        model2_f1[i] = np.max(model1_addi[0:20, :])

    if np.max(model1_full[1, :]) > np.max(model1_addi[20:40, :]):
        model2_f2[i] = np.max(model1_full[1, :])
    else:
        model2_f2[i] = np.max(model1_addi[20:40, :])




    tmp_model1 = np.concatenate((model1_f1, model1_f2), axis = 1); 
    tmp_model2 = np.concatenate((model2_f1, model2_f2), axis = 1);

    if np.all(tmp_model1 > tmp_model2):
        max_CF_f1[i] = model1_f1[i]; max_CF_f2[i] = model1_f2[i]
    else:
        max_CF_f1[i] = model2_f1[i]; max_CF_f2[i] = model2_f2[i]

    print(max_CF_f1[i], max_CF_f2[i])
    i += 1





max_all = np.concatenate((max_CF_f1, max_CF_f2), axis=1)
print(max_all)


# plotting the heatmap
y_values = np.arange(0.5, int(np.size((uniprot_name))) + 0.5, 1)
x_values = np.arange(0.5, 2.5, 1)
x_label = ['Fold1', 'Fold2']


plt.figure(figsize=(7,12))
hm = sns.heatmap(data=max_all,
                       annot=True, vmin=0.5, vmax=1, cmap="rocket_r", cbar_kws={'label': 'TM-score'})
hm.figure.axes[-1].set_ylabel('TM-score', size=20, rotation = 270, labelpad = 13)
hm.figure.axes[-1].tick_params(labelsize=15)


x_label = ['Fold1', 'Fold2']
plt.yticks(y_values, uniprot_name, fontsize= 13)
plt.yticks(rotation=0)
plt.xticks(x_values, x_label, fontsize = 15)
plt.tick_params(
    axis='x',          # changes apply to the x-axis
    bottom=False,      # ticks along the bottom edge are off
    top=False)         # ticks along the top edge are off


plt.tight_layout()
plt.savefig('heatmap-max TMscore comparison.png',dpi=600)

