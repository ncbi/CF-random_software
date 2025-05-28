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
import matplotlib.font_manager
from numpy import genfromtxt
from matplotlib import pyplot as plt
from adjustText import adjust_text
import glob



nsample = genfromtxt("number_of_predictions.csv", delimiter = ' ')

ref_Mc = sorted(glob.glob("ref_data_Mchaourab*csv"))

#### printing name
ref_name = []
for name in ref_Mc:
    ref_Mc_name = name.replace('.csv', '')
    ref_Mc_name = ref_Mc_name.replace('ref_data_Mchaourab_', '')

    ref_name = np.append(ref_name, ref_Mc_name)

print(nsample)

# plotting the heatmap
plt.figure(figsize=(7,5))
plt.rcParams["font.family"] = "Helvetica"
hm = sns.heatmap(data=nsample,
                annot=True, vmin=30, vmax=900, cmap="rocket_r", fmt="3.0f")
hm.figure.axes[-1].set_ylabel('Number of samples', size=15, rotation = 270, labelpad = 15)
y_values = np.arange(0.5, 14.5, 1)
x_values = np.arange(0.5, 2.5, 1)
#x_label = ['SPEACH F1', 'SPEACH F2', 'CF Fold1', 'CF Fold2']
x_label = ['SPEACH-AF', 'CF-random']
plt.xlabel("    ", fontsize = 13)
hm.xaxis.set_label_position("top") 
plt.yticks(y_values, ref_name, fontsize= 13)
plt.yticks(rotation=0)
plt.xticks(x_values, x_label, fontsize = 13)
plt.tick_params(
    axis='x',          # changes apply to the x-axis
    bottom=False,      # ticks along the bottom edge are off
    top=False)         # ticks along the top edge are off


plt.tight_layout()
plt.savefig('heatmap-nsamples.png',dpi=600, transparent=True)

