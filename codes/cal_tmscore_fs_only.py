#!/bin/env python3
# -*- coding: utf-8 -*-
"""
Compare the predicted models with original PDBs
report TM-scores for ranked 0 to 4
input line is pdb1 pdb2 preds_of_pdb dirname

This version requires tmtools 0.0.2 (Python bindings around the TM-align code for structural alignment of proteins)
check this for local installation
https://pypi.org/project/tmtools/

Usage:

python3.8 compare_strs_fs.py 2k42_A 1cee_B 1cee_B 0_msas_models/

Created on Wed Feb 21 14:51:00 2024
@author: Myeongsang (Samuel) Lee
"""
import re
import Bio
import os
from os import listdir
from os.path import isfile, join
import sys
from pathlib import Path
import numpy as np
from numpy import genfromtxt
import matplotlib.pyplot as plt
import glob
import random
import argparse
# call related modules of tmtools after installation
from tmtools import tm_align
from tmtools.io import get_structure, get_residue_data
from tmtools.testing import get_pdb_path
import Bio.PDB
from Bio.PDB import PDBParser, Structure




pdbParser = PDBParser(QUIET=True)

# convert three letter code to one letter code
aa3to1 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
 'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
 'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
 'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}


class TM_score_fs():
    def get_coords(self, pdbfile, fs_range):
            """
            parameters:
            pdbfile - path to pdbfile
            fs_range - range of residues at the fold-switching region, given as string - "112-162"
            returns:
            numpy array of coords
            string of seqs in 1-letter-code
            """
    
            seq = ""
            struct = pdbParser.get_structure('x',str(pdbfile))
            coords = []
            seq_dict = {}
    
            # for residues within a certain range, using numpy to save the coords
            # and save the sequence as a dict and then sorted list of tuples
            # return the coords and the seq
    
            # convert str to residue range for the fs region
            (start,stop) = fs_range.split("-")
            res_range = range(int(start),int(stop)+1)
    
            for atom in struct.get_atoms():
                    residue = atom.get_parent() # from atom we can get the parent residue
                    res_id = residue.get_id()[1]
                    resname = residue.get_resname()
                    if res_id in res_range and atom.get_name()=="CA":
                            x,y,z = atom.get_coord()
                            coords.append([x,y,z])
                            if res_id not in seq_dict:
                                    seq_dict[res_id]=aa3to1[resname]
    
    
            #print(coords)
            # convert to np array
            coords_np = np.array(coords)
            # sort the seq_dict by keys a.k.a res_ids
            sorted_data = sorted(seq_dict.items())
            for i in sorted_data:
                    seq+=i[1]
    
            return  coords_np,seq



    def get_tmscore(self, coords1, seq1, predfilepath, res_range):
           """
           parameters:
           coords1, seq1 - the numpy array of PDB coords and its seqs
           predfilepath - path for predicted files
           res_range - fs range in predicted models
    
           returns:
           tmscore list
    
           """
           
           tmscores = []
           tmscores_ord = []; tmscores_rev = []
           #modelfiles = sorted(glob.glob(str(predfilepath) + "/*_unrelaxed*pdb"))
           modelfiles = (glob.glob(str(predfilepath) + "/*_unrelaxed*pdb"))
    
           if len(modelfiles)==0:
                   tmscores = [0.0,0.0,0.0,0.0,0.0]
                   return tmscores
    
           for model in modelfiles:
                   modelpath = Path(model)
                   coords2, seq2 = self.get_coords(modelpath,res_range)
                   res = tm_align(coords1, coords2, seq1, seq2)
                   tmscore = round(res.tm_norm_chain1,2) # wrt to model
                   tmscores_ord.append(tmscore)

                   res = tm_align(coords2, coords1, seq2, seq1)
                   tmscore = round(res.tm_norm_chain1,5) # wrt to model
                   tmscores_rev.append(tmscore)

                   if np.max(tmscores_ord) > np.max(tmscores_rev):
                       tmscores = tmscores_ord
                   else:
                       tmscores = tmscores_rev

    
           return tmscores
    


    #def run_for_models(self, FH, pdbfile1, pdbfile2, data_dir,pred_range,res_range1,res_range2):
    def run_for_models(self, pdbfile1, pdbfile2, data_dir,pred_range,res_range1,res_range2):
            """
            compare the original PDB
            with the predicted models, 0 to 5
            
            parameters:
            FH - filehandle for writing
            pdbfile1 - path to original PDB, Fold1
            pdbfile2 - path to alternate PDB, Fold2
            data_dir - path for the predicted strs
            res_range1 - fs range in PDB1 and its models
            res_range2 - fs range in PDB2 and its models
             
            returns:
            nothing
             
            saves the TM-scores in a local file
            """
            #print(res_range1,res_range2)
    
            # get list of subdirectories
            all_sub_dir_paths = glob.glob(str(data_dir)) # returns list of sub directory paths
            tmscores_fs = [] 


            # files found then continue    
            if len(all_sub_dir_paths) == 0:
                pass
            
            for subdir in all_sub_dir_paths:
               preddir = Path(subdir)
               # predicted dir doesn't exist then continue
               if not preddir.exists():
                   pass
               
               # only comparing on one set of predicted models
               # but with both PDBs/Folds
               coords1,seq1 = self.get_coords(pdbfile1,res_range1)
               tmscore_lst1 = self.get_tmscore(coords1,seq1,preddir,pred_range) # wrt pdb1
               tmp_tm_fs = tmscore_lst1 
               tmscores_fs.append(tmp_tm_fs)
               
             
            for subdir in all_sub_dir_paths:
                preddir = Path(subdir)
                
                # predicted dir doesn't exist then continue
                if not preddir.exists():
                    pass
                
                # only comparing on one set of predicted models
                # but with both PDBs/Folds
                coords2,seq2 = self.get_coords(pdbfile2,res_range2)
                tmscore_lst2 = self.get_tmscore(coords2,seq2,preddir,pred_range) # wrt pdb2
                tmp_tm_fs = tmscore_lst2 
                tmscores_fs.append(tmp_tm_fs)
                
            print("         ")
            tmscores_fs = np.array(tmscores_fs)
            print("tmscores_fs")
            self.tmscores_fs = tmscores_fs
            


    def __init__(self, pred_path, pdb1, pdb1_name, pdb2, pdb2_name):
        # get numpy arrays for coords at the fold-switching region
        # also return the seq in 1-letter code for the same

        # input arguments: sys.argv[1] - pdb1, sys.argv[2] - pdb2
        #                  sys.argv[3] - preds, sys.argv[4] - current directory
        
        current_dir = os.getcwd() + '/' 
        #pred_dir = 'additional_sampling/' + pdb1_name  
        #pred_path = current_dir + pred_dir + '/'
        #print(pred_path)
        data_dir = Path(pred_path) # Path to the predicted models


        # the range of the fold-switching region
        range_file = current_dir + 'range_fs_pairs_all.txt'

        # convert this file into a dictionary for reference later
        fs_res = {}
        
        # The range_file file has the fold-switching residue ranges
        # for the original PDB/PDB1, alternate PDB/PDB2
        # Predicted model for PDB1, predicted model for PDB2
        with open(range_file,'r') as Infile:
                next(Infile) # skip header line "# pdb1,pdb2,pred1,pred2"
                for line in Infile:
                        line=line.strip()
                        (n1,n2,p1,p2,m1,m2)=line.split(",")
                        # the value of the dictionary is a tuple
                        # the first element of tuple is the fs range in the original PDB 
                        # followed by the range in the predicted model
                        if n1 not in fs_res:
                                fs_res[n1]=(p1,m1)
                        if n2 not in fs_res:
                                fs_res[n2]=(p2,m2)
        


        print("Running for pair ",pdb1_name, pdb2_name, end="..")
        print("                ")
        print("comparing predictions of ", pdb1_name, end="...")
        print("                ")


        try:
                range_pdb1 = fs_res[pdb1_name] # so if pdb1 is '1nqd_A', fs_res['1nqd_A']=('895-919', '1-33')
                range_pdb2 = fs_res[pdb2_name] # and if pdb2 is '1nqj_B', fs_res['1nqj_B']=('894-919', '1-33')
        except:
                print("check PDBIDs ",pdb1_name, pdb2_name)
                sys.exit(1)
        

        range_pred = range_pdb1[1]
        self.run_for_models(pdb1, pdb2, data_dir, range_pred, range_pdb1[0], range_pdb2[0])


#if __name__ == "__main__":
#
#    import warnings
#    warnings.filterwarnings('ignore')
#
#    parser = argparse.ArgumentParser()
#    parser.add_argument("--pdb1", type=str, help='PDB structure for the target crystal structure (target to be predicted)')
#    parser.add_argument("--pdb2", type=str, help='PDB structure for the alternative crystal structure')
#
#    args = parser.parse_args()
#
#    pdb1 = args.pdb1; pdb2 = args.pdb2
#    pdb1_name = pdb1.replace('.pdb','');  pdb2_name = pdb2.replace('.pdb','')
#
#    TM_score_fs(pdb1, pdb1_name, pdb2, pdb2_name)
#
