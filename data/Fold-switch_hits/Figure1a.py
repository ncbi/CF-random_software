#! /Users/porterll/miniconda3/bin/python

import sys
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
import matplotlib as mpl

if __name__ == '__main__':

    names = ['AF-cluster','SPEACH_AF','All-AF (no templates)','CF-random']

    mpl.rcParams['font.family'] = 'sans-serif'
    mpl.rcParams['font.sans-serif'] = ['Helvetica']
    mpl.rcParams["axes.labelsize"] = 14

    num_models=[202930,76000,283730,47000]
    #Added 1 to All-AF because we missed 4HDD before
    num_successes = [18,7,25,32]

    plt.plot(num_models,num_successes,'ko')
    plt.xlabel('#Predictions',size=16)
    plt.ylabel('#Successes',size=16)
    plt.savefig('method_performances.png',dpi=600)
