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





data_name = open('list_of_OC23-uniprot_ID.csv', 'r')
myreader = csv.reader(data_name)

uniprot_name = []; pdb1_name = []; pdb2_name = []; nsampleO = []; nsampleC = []

for row in myreader:
    uniprot_name = np.append(uniprot_name, row[0])
    pdb1_name = np.append(pdb1_name, row[1])
    pdb2_name = np.append(pdb2_name, row[2])
    nsampleO = np.append(nsampleO, row[3])
    nsampleC = np.append(nsampleC, row[4])
    
uniprot_name = uniprot_name[1:]
print(uniprot_name)
pdb1_name = pdb1_name[1:]
pdb2_name = pdb2_name[1:]
nsampleO = nsampleO[1:]
nsampleC = nsampleC[1:]
uniprot_size = np.size(uniprot_name)


pdb1_name = pdb1_name.ravel().astype(str)
pdb2_name = pdb2_name.ravel().astype(str)
#uniprot_name = uniprot_name.ravel().astype(list)
nsampleO = nsampleO.ravel().astype(int)
nsampleC = nsampleC.ravel().astype(int)


print(nsampleO, nsampleC)


pdb1_name = pdb1_name.reshape(uniprot_size, 1); pdb2_name = pdb2_name.reshape(uniprot_size, 1)
nsampleO = nsampleO.reshape(uniprot_size, 1); nsampleC = nsampleC.reshape(uniprot_size, 1);

nsample_all = np.concatenate((nsampleO, nsampleC), axis = 1)



# plotting the heatmap
plt.figure(figsize=(7,5))
plt.rcParams["font.family"] = "Helvetica"
hm = sns.heatmap(data=nsample_all,
                annot=True, vmin=100, vmax=1000, cmap="rocket_r", fmt="3.0f")
hm.figure.axes[-1].set_ylabel('Number of samples', size=15, rotation = 270, labelpad = 15)

y_values = np.arange(0.5, int(uniprot_size) + 0.5, 1)
x_values = np.arange(0.5, 2.5, 1)
#x_label = ['SPEACH F1', 'SPEACH F2', 'CF Fold1', 'CF Fold2']
x_label = ['AFsample2', 'CF-random']
plt.yticks(y_values, uniprot_name, fontsize= 13)
plt.yticks(rotation=0)
plt.xticks(x_values, x_label, fontsize = 15)
plt.tick_params(
    axis='x',          # changes apply to the x-axis
    bottom=False,      # ticks along the bottom edge are off
    top=False)         # ticks along the top edge are off

plt.tight_layout()
plt.savefig('OC23_heatmap-nsamples.png',dpi=600, transparent=True)
plt.savefig('OC23_heatmap-nsamples.svg',dpi=600, transparent=True)
 
