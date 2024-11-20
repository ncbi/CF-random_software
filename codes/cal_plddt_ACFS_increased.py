"""
find the average pLDDT score
"""
import sys,re
import glob
import json
import numpy as np
from pathlib import Path

#define pattern for regular expression

# 0_000_scores_rank_001_alphafold2_ptm_model_4_seed_000.json
pattrn = re.compile(r'.*?_scores_rank_(?P<rank>\d+)_alphafold2.*')

# default if pattern doesn't work
rank = "000"

def read_plddt (jsonfile):
    """
    read the json file
    return the plddt scores
    as numpy array
    """
    with open(jsonfile) as json_file:
        data = json.load(json_file)
    
    plddt_scores = np.array(data['plddt'],dtype='float64')

    return plddt_scores

def fract_good (score):
    """
    return percentage
    of residue with
    plddt score > 70
    """
    vals_greater_70 = (score > 70).sum()
    percent_good = round((vals_greater_70 / score.size)*100,2)
    avg_plddt = round(np.average(score),2)
    #return percent_good,avg_plddt
    return avg_plddt




class plddt_cal_inc():
    def __init__(self, sub_list, category, pdb_name):
        # if files found then continue     
        if len(sub_list) == 0:
            sys.exit(1)
            
        # create a data dictionary
        out_dict_all = {}
        
        values_all = []
        cnt = 0

        if category =='full-MSA':
        #if category == 'additional-MSA':
            print("working...")            
            print(sub_list)
            for subdir in sub_list:
                print(subdir)
                if Path(subdir).is_dir():
                    subdir_name = Path(subdir).name
                    jsonfiles = glob.glob(str(subdir) + "/*_scores*json")

                    for jsonfile in jsonfiles:
                        plddt_score = read_plddt(jsonfile)
                        values = fract_good(plddt_score)
                        values_all = np.append(values_all, values)
                        jsonfilepath = Path(jsonfile)
                        jsonfilename = jsonfilepath.stem
                        match = pattrn.match(jsonfilename)
                        if match:
                            rank = match.group('rank')

                        key_pair = subdir_name + ":" + rank
                        # for all
                        if key_pair not in out_dict_all:
                            out_dict_all[key_pair]=values

                        cnt = cnt + 1
            cnt = int(cnt / 5)



        elif category == 'additional-MSA':
            print("working...")
            print(sub_list)
            for subdir in sub_list:
                print(subdir)
                if Path(subdir).is_dir():
                    subdir_name = Path(subdir).name
                    jsonfiles = glob.glob(str(subdir) + "/*_scores*json")

                    for jsonfile in jsonfiles:
                        plddt_score = read_plddt(jsonfile)
                        values = fract_good(plddt_score)
                        values_all = np.append(values_all, values)
                        jsonfilepath = Path(jsonfile)
                        jsonfilename = jsonfilepath.stem
                        match = pattrn.match(jsonfilename)
                        if match:
                            rank = match.group('rank')

                        key_pair = subdir_name + ":" + rank
                        # for all
                        if key_pair not in out_dict_all:
                            out_dict_all[key_pair]=values

                        cnt = cnt + 1

            

        else:
            for subdir in sub_list:
            #for subdir in all_sub_dir_paths:
                # make sure subdir exists
                if Path(subdir).is_dir():
                    subdir_name = Path(subdir).name
                    # get the list of json files
                    jsonfiles = glob.glob(str(subdir) + "/*_scores*json")
                    for jsonfile in jsonfiles:
                        plddt_score = read_plddt(jsonfile)
                        values = fract_good(plddt_score)
                        values_all = np.append(values_all, values)
                        jsonfilepath = Path(jsonfile)
                        jsonfilename = jsonfilepath.stem
                        match = pattrn.match(jsonfilename)
                        if match:
                            rank = match.group('rank')   
                        
                        key_pair = subdir_name + ":" + rank
                        # for all
                        if key_pair not in out_dict_all:
                            out_dict_all[key_pair]=values

                cnt = cnt + 1


        print(cnt)
        print(values_all)
        if category =='full-MSA':
            values_all_resh = values_all.reshape(25, 5)
        elif category == 'additional-MSA':
            values_all_resh = values_all.reshape(30, 5)
        else:
            values_all_resh = values_all.reshape(175, 5)
        print("                ")
        print("Calculated pLDDT")
        print(values_all_resh)
        np.savetxt('plddt_' + category + '_' + pdb_name +'.csv', values_all_resh, fmt='%2.3f')
