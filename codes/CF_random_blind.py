import os, glob, re
from os.path import isfile
import sys
import numpy as np
import math
from tmtools import tm_align
from tmtools.io import get_structure, get_residue_data
from tmtools.testing import get_pdb_path
from matplotlib import pyplot as plt
import matplotlib as mpl
from matplotlib import cm
import pandas as pd


class blind_screening():
    def TMscore_blind(self, pdb1, pdb2, tmscores):
    
         ## loading reference pdb for TM-score
         pwd = os.getcwd() + '/'
    
         ##### pdb1_name part
         pdb1_dir = pwd + pdb1
         model1 = pdb1_dir.replace('.pdb','') 
         r2 = get_structure(get_pdb_path(str(model1)))
         coords2, seq2 = get_residue_data(r2)
            
         pdb2_dir = pwd + pdb2
         model2 = pdb2_dir.replace('.pdb','')
         s = get_structure(get_pdb_path(model2))
         coords1, seq1 = get_residue_data(s)
         res = tm_align(coords1, coords2, seq1, seq2)
         tmscore = round(res.tm_norm_chain1,5) # wrt to model
         tmscores.append(tmscore)
    
    def get_reference_ID(self, IDs, pdb1_name):
    
         max_score = 0
         max_ID = ''
         
         with open(IDs) as log:
              f = log.read().splitlines()
              scores = np.array([float(x.split()[4].split('=')[1]) for x in f if 'recycle=3' in x])
              pdbs = [x.split()[2] for x in f if 'recycle=3' in x]
    
         print(scores)
         print(pdbs)
    
         idx = np.argmax(scores)

         return IDs.split('/')[0] + '/' + IDs.split('/')[1] + '/' + IDs.split('/')[2] + '/0_unrelaxed_'+pdbs[idx]+'.pdb'
    


    def __init__(self, pdb1_name, blind_path):
    
        print(pdb1_name)
        ofname = os.getcwd().split('/')[-1]
        

        pwd = os.getcwd() + '/'
    
        ref = self.get_reference_ID(glob.glob(blind_path + '/*full*/log.txt')[0], pdb1_name)
        print(ref)
    
        tmscores = []
        foldseek_pdb = []
        foldseek_val= []
        pdbs_used = []
    
        pdb_files = glob.glob(blind_path + '/*max*/*pdb')
    
        fs = ref.split('.')[0]+'.foldseek'
    
        with open(fs) as f1:
            info = f1.read().splitlines()[0].split()
            pdb_ref = info[1]
            fs_ref = -math.log(float(info[-2]))
    
        for p in pdb_files:
    
            if not os.path.isfile(p.split('.')[0]+'.foldseek'):
                 continue
             
            self.TMscore_blind(ref,p,tmscores)
            with open(p.split('.')[0]+'.foldseek') as f1:
                info2 = f1.read().splitlines()
                if len(info2) < 1:
                     tmscores = tmscores[:-1]
                     continue
                info = re.split('\t |\s',info2[0])
                foldseek_pdb.append(info[1][:5]+info[1][-1])
                foldseek_val.append(-math.log(float(info[-2]))/fs_ref)
                pdbs_used.append(p)
    
        pdbs_used.append(ref)
        tmscores.append(1.0)
        foldseek_pdb.append(pdb_ref[:5]+pdb_ref[-1])
        foldseek_val.append(1.0)
    
        results = pd.DataFrame({'pdb_file':pdbs_used,
                                'FoldSeek_pdb':foldseek_pdb,
                                'TM-score':tmscores,
                                'FoldSeek_score':foldseek_val})
    
        good_results = results[results['FoldSeek_score'] >= 0.65]
    
        
        #good_results.to_csv(ofname+'.csv')
        good_results.to_csv(pdb1_name + '.csv')
     
        unique_IDs = []
        for i in foldseek_pdb:
             if i not in unique_IDs:
                  unique_IDs.append(i)
                  
        unique_idxs= range(len(unique_IDs))
        unique_info = {unique_IDs[i]:unique_idxs[i]for i in range(len(unique_IDs))}
        fs_nums = np.array([unique_info[x] for x in foldseek_pdb])
        fig, ax = plt.subplots()
    
        colors = mpl.cm.rainbow(np.linspace(0,1,len(unique_IDs)))
    
        leg_labs = []
    
        good_tms = []
        good_vals = []
        good_pdbs = []
        good_nums = []
        good_models = []
        all_pdbs = []
    
        pdb_files.append(ref)
        
        for i,v in enumerate(foldseek_val):
             if v >= 0.65:
                  good_tms.append(tmscores[i])
                  good_vals.append(v)
                  good_nums.append(fs_nums[i])
                  good_models.append(pdb_files[i])
                  all_pdbs.append(foldseek_pdb[i])
                  if foldseek_pdb[i] not in good_pdbs:
                       good_pdbs.append(foldseek_pdb[i])
    
        best_nums = [0]*len(unique_IDs)
        best_pdbs = [0]*len(unique_IDs)
        best_TMs = [0]*len(unique_IDs)
    
        for i,v in enumerate(foldseek_val):
                  if v > best_nums[fs_nums[i]]:
                       best_nums[fs_nums[i]] = v
                       best_pdbs[fs_nums[i]] = pdb_files[i]
                       best_TMs[fs_nums[i]] = tmscores[i]
    
        #best_hits = open(ofname+'_best_hits.txt',"w")
        best_hits = open(pdb1_name + '_best_hits.txt', "w")
    
        diffs = np.array(best_nums)-np.array(best_TMs)
    
        for i in range(len(unique_IDs)):
             if best_nums[i] >= 0.65:
                  best_hits.write('%-100s %-30s %1.2f %1.2f %1.4f\n' %(best_pdbs[i],unique_IDs[i],best_nums[i],best_TMs[i],diffs[i]))
    
        best_hits.close()
    
        hits = pd.DataFrame({'TMs':good_tms,
                             'Vals':good_vals,
                             'IDs':all_pdbs})
    
        for i in good_pdbs:
    
             df = hits[hits['IDs']==i]
             scatter = ax.scatter(df['TMs'],df['Vals'],label=i)
    
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.85, box.height])
    
        h,l = ax.get_legend_handles_labels()
        ax.legend(h,l,loc='center left',bbox_to_anchor=(1,0.5))
        
        plt.xlabel('TM-score')
        plt.ylabel('Normalized foldseek score')
        #plt.savefig(ofname+'.png',dpi=600)
        plt.savefig(pdb1_name + '.png', dpi=600)
    
        plt.clf()

    


    
