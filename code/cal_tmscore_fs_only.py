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

from Bio import SeqIO 
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1

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
    def find_resi_index(self, pdb, target_seq):
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("protein", pdb)                                   

        # Extract residues from the first model                                            
        model = structure[0]                                                               
        chain = next(model.get_chains())
        #chain = model[chain_id]                                                            

        residues = [r for r in chain.get_residues() if r.id[0] == " "]  # exclude heteroatoms
        # Build sequence and residue index list
        sequence = "".join(seq1(r.resname) for r in residues)
        res_indices = [r.id[1] for r in residues]

        # Find the subsequence
        start = sequence.find(target_seq)
        if start == -1:
            print("Subsequence not found in chain.")
            return None

        end = start + len(target_seq)
        matched_indices = res_indices[start:end]

        return matched_indices 

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
            print(tmscores_fs)
            self.tmscores_fs = tmscores_fs
            


    def __init__(self, pred_path, pdb1, pdb1_name, pdb2, pdb2_name, seq):
        # get numpy arrays for coords at the fold-switching region
        # also return the seq in 1-letter code for the same

        # input arguments: sys.argv[1] - pdb1, sys.argv[2] - pdb2
        #                  sys.argv[3] - preds, sys.argv[4] - current directory
        
        current_dir = os.getcwd() + '/' 
        #pred_dir = 'additional_sampling/' + pdb1_name  
        #pred_path = current_dir + pred_dir + '/'
        #print(pred_path)
        data_dir = Path(pred_path) # Path to the predicted models




        print("Running for pair ",pdb1_name, pdb2_name, end="..")
        print("                ")
        print("comparing predictions of ", pdb1_name, end="...")
        print("                ")



        #range_pred = range_pdb1[1]
        pred_path = glob.glob(pred_path + '/*pdb')
        pred = pred_path[0]
        pdb1_index = self.find_resi_index(pdb1, seq)
        pdb2_index = self.find_resi_index(pdb2, seq)
        pred_index = self.find_resi_index(pred, seq)
        print(pdb1_index, pdb2_index, pred_index)
        range_pdb1 = str(pdb1_index[0]) + '-' + str(pdb1_index[-1])
        range_pdb2 = str(pdb2_index[0]) + '-' + str(pdb2_index[-1])
        range_pred = str(pred_index[0]) + '-' + str(pred_index[-1])
        print(range_pdb1, range_pdb2, range_pred)

        self.run_for_models(pdb1, pdb2, data_dir, range_pred, range_pdb1, range_pdb2)


