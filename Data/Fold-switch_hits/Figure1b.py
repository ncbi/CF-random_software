import pandas as pd
from matplotlib import pyplot as plt
import matplotlib as mpl

import seaborn as sns

if __name__ == '__main__':

    df = pd.read_csv('Hits_CF_parameters_fold2.csv')

    counts = df['max_seq'].value_counts()

    #fig,ax = plt.subplots()

    mpl.rcParams['font.family'] = 'sans-serif'
    mpl.rcParams['font.sans-serif'] = ['Helvetica']
    mpl.rcParams["axes.labelsize"] = 14

    g = sns.barplot(counts,color='k')

    g.set_xticklabels(['1:0','1:2','2:4','4:8','8:16','16:32','32:64','64:128','Full'])

    plt.ylabel('Counts',size=16)
    plt.xlabel('Sampling depths',size=16)

    plt.ylim([0,10])

    plt.savefig('barplot.png',dpi=600)

    

    
