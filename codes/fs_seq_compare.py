import sys
import os
import re
from os import listdir
from os.path import isfile, join
import pandas as pd
import numpy as np
import Bio.PDB
import matplotlib.pyplot as plt
import glob
import random
import argparse
from Bio import SeqIO
from Bio.PDB.PDBParser import PDBParser
from Bio import *
from Bio.SeqRecord import SeqRecord

from thefuzz import fuzz
from thefuzz import process



class fs_range():
    def first_res_check(self, pdb1, pdb2):
        #self.pdb1 = pdb1; self.pdb2 = pdb2

        ## first residue index check
        structure_1 = PDBParser().get_structure('pdb1', pdb1)
        model_1 = structure_1[0]
        print(model_1)
        
        structure_2 = PDBParser().get_structure('pdb2', pdb2)
        model_2 = structure_2[0]
        print(model_2)
        
        res_index_1 = []
        res_index_2 = []
        
        for chain_1 in model_1:
            for i, residue in enumerate(chain_1.get_residues()):
                #res_id = list(residue.id)
                res_index_1.append(residue.id[1])
                #print(residue.id[1])
        
        for chain_2 in model_2:
            for i, residue in enumerate(chain_2.get_residues()):
                res_index_2.append(residue.id[1])
        
        #print(int(res_index_1[0]))
        #print(int(res_index_2[0]))

        self.pdb1_res_index_1 = int(res_index_1[0])
        self.pdb2_res_index_1 = int(res_index_2[0])




    def pydssp(self, crys_pdb, pred_pdb, number, pdb_name):

        ##### generating the command for pydssp
        number = str(number)
        command = 'pydssp ' + crys_pdb + ' ' + pred_pdb + ' -o output_' + pdb_name + '_' + number + '.log'
        print(command)
        os.system(command)
            



    def res_check(self, pdb1, pdb2, pdb1_name, pdb2_name):
        current_dir = os.getcwd() + '/'
        range_file = current_dir + 'range_fs_pairs_all.txt'
    
        crys_fs_res_1 = {}; crys_fs_res_2 = {}
        pred_fs_res_1 = {}; pred_fs_res_2 = {}

        with open(range_file,'r') as Infile:
            next(Infile) # skip header line "# pdb1,pdb2,pred1,pred2"
            for line in Infile:
                line=line.strip()
                (n1,n2,p1,p2,m1,m2)=line.split(",")
                # the value of the dictionary is a tuple
                # the first element of tuple is the fs range in the original PDB
                # followed by the range in the predicted model
                #if n1 == pdb1_name and n2 == pdb2_name:
                if (n1 == pdb1_name and n2 == pdb2_name) or (n2 == pdb1_name and n1 == pdb2_name):
                    #fs_res_1 = (m1); fs_res_2 = (m2)
                    crys_fs_res_1 = (p1); crys_fs_res_2 = (p2);
                    pred_fs_res_1 = (m1); pred_fs_res_2 = (m2);
    
            #fs_res_1_update = fs_res_1.split("-"); fs_res_2_update = fs_res_2.split("-");
            #print(fs_res_1_update, fs_res_2_update)
    
    
        crys_fs_res_1_update = crys_fs_res_1.split("-"); crys_fs_res_2_update = crys_fs_res_2.split("-");
        print(crys_fs_res_1_update, crys_fs_res_2_update)
        pred_fs_res_1_update = pred_fs_res_1.split("-"); pred_fs_res_2_update = pred_fs_res_2.split("-");
        print(pred_fs_res_1_update, pred_fs_res_2_update)
    
        ##### convert list data to int
        self.crys_fs_res_1_update = [int(i) for i in crys_fs_res_1_update]
        self.crys_fs_res_2_update = [int(i) for i in crys_fs_res_2_update]

        self.pred_fs_res_1_update = [int(i) for i in pred_fs_res_1_update]
        self.pred_fs_res_2_update = [int(i) for i in pred_fs_res_2_update]







    def __init__(self, pdb1, pdb2, pdb1_name, pdb2_name, pred_dir):
        ##### check first residue index of query proteins
        #fs_check = fs_range(pdb1, pdb2)
        self.first_res_check(pdb1, pdb2)
        print("    "); print("checking first residue index")
        print(self.pdb1_res_index_1)
        print(self.pdb2_res_index_1) 
    
        
    
        pred_folder = pred_dir
        #pred_folder = '3hdf_A_predicted_models_full_rand_12'
        pred_path = pred_folder
        print(pred_path)
    
    
        pred_files = (glob.glob(str(pred_path) + "/*_relaxed*pdb"))
    
    
        ##### read range file information
        self.res_check(pdb1, pdb2, pdb1_name, pdb2_name)
        print(self.crys_fs_res_1_update, self.pred_fs_res_1_update)
        print(self.crys_fs_res_2_update, self.pred_fs_res_2_update)
    
        crys1_fs_res_st = self.crys_fs_res_1_update[0]; crys1_fs_res_ed = self.crys_fs_res_1_update[1]
        crys2_fs_res_st = self.crys_fs_res_2_update[0]; crys2_fs_res_ed = self.crys_fs_res_2_update[1]
        pred1_fs_res_st = self.pred_fs_res_1_update[0]; pred1_fs_res_ed = self.pred_fs_res_1_update[1]
        pred2_fs_res_st = self.pred_fs_res_2_update[0]; pred2_fs_res_ed = self.pred_fs_res_2_update[1]
    
    
    
        if int(self.pdb1_res_index_1) > 1:                                                              
            print("Initial residue is not starting from 1")
            self.crys_fs_res_1_update[0] = self.crys_fs_res_1_update[0] - int(self.pdb1_res_index_1)                            
            self.crys_fs_res_1_update[1] = self.crys_fs_res_1_update[1] - int(self.pdb1_res_index_1)                            
            crys1_fs_res_st = self.crys_fs_res_1_update[0]; 
            crys1_fs_res_ed = self.crys_fs_res_1_update[1]
                                                                                                 
        if int(self.pdb2_res_index_1) > 1:                                                              
            print("Initial residue is not starting from 1")
            self.crys_fs_res_2_update[0] = self.crys_fs_res_2_update[0] - int(self.pdb2_res_index_1)                            
            self.crys_fs_res_2_update[1] = self.crys_fs_res_2_update[1] - int(self.pdb2_res_index_1)                            
            crys2_fs_res_st = self.crys_fs_res_2_update[0]
            crys2_fs_res_ed = self.crys_fs_res_2_update[1]
                                                                                                 
        print("checking starting and ending residue number")
        print(""); print("crystal structure")
        print(crys1_fs_res_st, crys1_fs_res_ed)
        print(crys2_fs_res_st, crys2_fs_res_ed)
        print(""); print("predicted structure")
        print(pred1_fs_res_st, pred1_fs_res_ed)
        print(pred2_fs_res_st, pred2_fs_res_ed)


        pred_dir_add = 'additional_sampling/' + pdb1_name + '/'
        pred_dir_suc = 'successed_prediction/' + pdb1_name + '/*/'
        pred_dir_fal = ' failed_prediction/'
    
    
        ##### perform pydssp and calculate secondary structure similarity
        index = 0
        print(np.size(pred_files))
        print("         "); print("calculating with pdb1 ", pdb1_name)
        for model in pred_files:
            print(model)
            self.pydssp(pdb1, model, index, pdb1_name)
            dssp_read_tmp = pd.read_csv('output_' + pdb1_name + '_' + str(index) + '.log', sep=' ', header = None)
            ## seq1 = crystal structure, seq2 = predicted structure
            print(dssp_read_tmp)
            print(dssp_read_tmp[0].iloc[0]); seq1 = dssp_read_tmp[0].iloc[0]
            print(dssp_read_tmp[0].iloc[1]); seq2 = dssp_read_tmp[0].iloc[1]
    
            # crystal protein 1 and predictions
            print("     ")
            print(seq1[crys1_fs_res_st:crys1_fs_res_ed])
            print(seq2[pred2_fs_res_st:pred2_fs_res_ed])
            if fuzz.ratio(seq1[crys1_fs_res_st:crys1_fs_res_ed], seq2[pred2_fs_res_st:pred2_fs_res_ed]) > 85:
                print("fs region is correctly predicted")
                f = open("fs_compare_output_" + pdb1_name + ".log", "w")
                f.write("success")
                f.close()
                break
            elif index == (int(np.size(pred_files)) - 1):
                print("fs region is not correctly predicted")

                #command = 'mv ' + pred_dir_add + pred_dir_fal
                #print(command); os.system(command)
                #command = 'mv ' + pred_dir_suc + pred_dir_fal + pdb1_name + '/'
                #print(command); os.system(command)

                #command = 'rm *' + pdb1_name + '*csv'
                #print(command); os.system(command)
                print("calculating TM-score of fs with alternative pdb")

                index = 0
                print("         "); print("calculating with pdb2 ", pdb2_name)

                for model in pred_files:
                    self.pydssp(pdb2, model, index, pdb1_name)
                    dssp_read_tmp = pd.read_csv('output_' + pdb1_name + '_' + str(index) + '.log', sep=' ', header = None)
                    ## seq1 = crystal structure, seq2 = predicted structure
                    print(dssp_read_tmp[0].iloc[0]); seq1 = dssp_read_tmp[0].iloc[0]
                    print(dssp_read_tmp[0].iloc[1]); seq2 = dssp_read_tmp[0].iloc[1]


                    # crystal protein 1 and predictions
                    print("     ")
                    print(seq1[crys2_fs_res_st:crys2_fs_res_ed])
                    print(seq2[pred2_fs_res_st:pred2_fs_res_ed])
                    if fuzz.ratio(seq1[crys2_fs_res_st:crys2_fs_res_ed], seq2[pred2_fs_res_st:pred2_fs_res_ed]) > 85:
                        print("fs region is correctly predicted")
                        break
                    elif index == (int(np.size(pred_files)) - 1):
                        print("fs region is not correctly predicted")

                        f = open("fs_compare_output_" + pdb1_name + ".log", "w")
                        f.write("fail")
                        f.close()

                        #command = 'mv ' + pred_dir_add + pred_dir_fal
                        #print(command); os.system(command)
                        #command = 'mv ' + pred_dir_suc + pred_dir_fal + pdb1_name + '/'
                        #print(command); os.system(command)

                    else:
                        index += 1


            else:
                index += 1
    
            # index += 1
    
    
        #index = 0
        #print("         "); print("calculating with pdb2 ", pdb2_name)
        #for model in pred_files:
        #    self.pydssp(pdb2, model, index)
        #    dssp_read_tmp = pd.read_csv('output_' + str(index) + '.log', sep=' ', header = None)
        #    ## seq1 = crystal structure, seq2 = predicted structure
        #    print(dssp_read_tmp[0].iloc[0]); seq1 = dssp_read_tmp[0].iloc[0]
        #    print(dssp_read_tmp[0].iloc[1]); seq2 = dssp_read_tmp[0].iloc[1]
    
        #    # crystal protein 1 and predictions
        #    print("     ")
        #    print(seq1[crys2_fs_res_st:crys2_fs_res_ed])
        #    print(seq2[pred2_fs_res_st:pred2_fs_res_ed])
        #    if fuzz.ratio(seq1[crys2_fs_res_st:crys2_fs_res_ed], seq2[pred2_fs_res_st:pred2_fs_res_ed]) > 85:
        #        print("fs region is correctly predicted")
        #        break
        #    elif index == (int(np.size(pred_files)) - 1):
        #        print("fs region is not correctly predicted")

        #        command = 'mv ' + pred_dir_add + pred_dir_fal
        #        print(command); os.system(command)
        #        command = 'mv ' + pred_dir_suc + pred_dir_fal + pdb1_name + '/'
        #        print(command); os.system(command)

        #        #command = 'rm *' + pdb1_name + '*csv'
        #        #print(command); os.system(command)



        #    else:
        #        index += 1


