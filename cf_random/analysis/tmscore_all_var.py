#!/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 21 14:51:00 2024

@author: Myeongsang (Samuel) Lee
"""

import glob
import os
import sys

import numpy as np

# call related modules of tmtools after installation
from tmtools import (
    tm_align,
)
from tmtools.io import (
    get_residue_data,
    get_structure,
)
from tmtools.testing import (
    get_pdb_path,
)

from ..utils.convert_multi_single import (
    convert_m2s,
)


class TMScore:
    def __init__(self, pred_dir, pdb1, pdb1_name, pdb2, pdb2_name, model_type):

        ## loading reference pdb for TM-score
        pwd = os.getcwd() + "/"
        tmscores = []
        tmscores_ord = []
        tmscores_rev = []

        # files_list = sorted(glob.glob(str(pred_dir) + "/*_unrelaxed*pdb"))
        if model_type != "alphafold2_multimer_v3":
            files_list = glob.glob(str(pred_dir) + "/*_unrelaxed*pdb")
            print(files_list)
        else:
            # Convert multimer prediction to single-chain PDBs if needed
            check_files_list = glob.glob(str(pred_dir) + "/rmTER*_unrelaxed*pdb")
            print(check_files_list)
            if not check_files_list:
                convert_m2s(pred_dir, pdb1_name, pdb2_name)
                files_list = glob.glob(str(pred_dir) + "/rmTER*_unrelaxed*pdb")
                print(files_list)
            else:
                files_list = glob.glob(str(pred_dir) + "/rmTER*_unrelaxed*pdb")
                print(files_list)

        # Reference: pdb1
        pdb1_dir = pwd + pdb1_name
        r2 = get_structure(get_pdb_path(str(pdb1_dir)))
        coords2, seq2 = get_residue_data(r2)

        if len(files_list) == 0:
            tmscores = [0.0, 0.0, 0.0, 0.0, 0.0]
            return tmscores

        for model in files_list:
            # modelpath = Path(model)
            # model  = str(modelpath.parent) + "/" + modelpath.stem
            model = model.replace(".pdb", "")
            # model = model.replace('_converted.pdb','_converted')
            model = pwd + model
            s = get_structure(get_pdb_path(model))
            coords1, seq1 = get_residue_data(s)
            res = tm_align(coords1, coords2, seq1, seq2)
            tmscore = round(res.tm_norm_chain1, 5)  # wrt to model
            tmscores_ord.append(tmscore)

            res = tm_align(coords2, coords1, seq2, seq1)
            tmscore = round(res.tm_norm_chain1, 5)  # wrt to model
            tmscores_rev.append(tmscore)

        # Reference: pdb2
        pdb2_dir = pwd + pdb2_name
        r3 = get_structure(get_pdb_path(str(pdb2_dir)))
        coords2, seq2 = get_residue_data(r3)

        for model in files_list:
            # modelpath = Path(model)
            # model  = str(modelpath.parent) + "/" + modelpath.stem
            model = model.replace(".pdb", "")
            # model = model.replace('_converted.pdb','_converted')
            model = pwd + model
            s = get_structure(get_pdb_path(model))
            coords1, seq1 = get_residue_data(s)
            res = tm_align(coords1, coords2, seq1, seq2)
            tmscore = round(res.tm_norm_chain1, 5)  # wrt to model
            tmscores_ord.append(tmscore)

            res = tm_align(coords2, coords1, seq2, seq1)
            tmscore = round(res.tm_norm_chain1, 5)  # wrt to model
            tmscores_rev.append(tmscore)

        print("normal")
        print(tmscores_ord)
        print("reverse")
        print(tmscores_rev)
        if np.max(tmscores_ord) > np.max(tmscores_rev):
            tmscores = tmscores_ord
        else:
            tmscores = tmscores_rev

        print(tmscores)
        self.tmscores = tmscores

    def select_size(self, TMscores_random_alter, pdb1_name, pdb2_name, alt_name, num_seeds):

        TMscores_random_reshape = TMscores_random_alter.reshape(14, num_seeds * 5)
        TMscores_random_locat = np.zeros((7, num_seeds * 5))

        #### finding locatnative pdb_name

        if alt_name == pdb2_name:
            # for i in 1, 3, 5, 7, 9, 11, 13 in TM_scores:
            tmp_cnt = 0
            for i in range(1, 14, 2):
                print(TMscores_random_reshape[i, :])
                TMscores_random_locat[tmp_cnt, :] = TMscores_random_reshape[i, :]
                tmp_cnt = tmp_cnt + 1
        else:
            # for i in 0, 2, 4, 6, 8, 10, 12 in TM_scores:
            tmp_cnt = 0
            for i in range(0, 13, 2):
                print(TMscores_random_reshape[i, :])
                TMscores_random_locat[tmp_cnt, :] = TMscores_random_reshape[i, :]
                tmp_cnt = tmp_cnt + 1

        TMscore_data = TMscores_random_locat
        TMscore_data = TMscores_random_locat.reshape(7, num_seeds * 5)
        TMscore_data_sum = np.zeros((7, 1))

        for ii in range(0, int(TMscore_data.shape[0])):
            TMscore_data_sum[ii] = np.sum(TMscore_data[ii])

        location = np.argmax(np.max(TMscore_data_sum, axis=1))

        print("Selecting...")

        TMscore_data = TMscores_random_alter
        TMscore_data = TMscores_random_alter.reshape(14, num_seeds * 5)

        location_org = location

        if alt_name == pdb2_name:
            location = (location * 2) + 1
        else:
            location = location * 2

        if alt_name == pdb2_name and np.any(TMscore_data[location, :] >= 0.5):
            print(TMscore_data[location, :])
            selection = int((location - 1) / 2)
            self.selection = selection

        elif alt_name == pdb1_name and np.any(TMscore_data[location, :] >= 0.5):
            print(TMscore_data[location, :])
            selection = int(location / 2)
            self.selection = selection

        else:
            print("Predictions are bad")
            print("Predictions of whole structure are bad")
            rm_folder_cmd = "rm -rf successed_prediction/" + self.pdb1_name + "/"
            print(rm_folder_cmd)
            os.system(rm_folder_cmd)
            sys.exit()


class TMScoreCalAllVar:
    def __init__(self, pdb1, pdb1_name, pdb2, pdb2_name, nMSA, option, model_type):
        num_seeds = 5 + nMSA
        pwd = os.getcwd() + "/"

        if model_type != "alphafold2_multimer_v3":
            # Evaluate TM-scores from predictions using full-length MSAs (whole structure)
            pred_dir = (
                "predictions_all/" + pdb1_name + "/" + pdb1_name + "_predicted_models_full_rand_*"
            )
            MSA_full_TMscore = TMScore(pred_dir, pdb1, pdb1_name, pdb2, pdb2_name, model_type)
            full_TMscore = np.array(MSA_full_TMscore.tmscores)
            full_TMscore = full_TMscore.reshape(2, num_seeds * 5)

            # Determine reference vs alternative structure based on full-MSA TM-scores
            if np.any(full_TMscore[0, :] > 0.5) or np.any(full_TMscore[1, :] > 0.5):
                if np.average(full_TMscore[0, :]) > np.average(full_TMscore[1, :]):
                    ref_name = pdb1_name
                    alt_name = pdb2_name
                else:
                    ref_name = pdb2_name
                    alt_name = pdb1_name
            elif np.all(full_TMscore[0, :] < 0.5) and np.all(full_TMscore[1, :] < 0.5):
                # If prediction is failed, move the folder to "failed_prediction""
                gen_dir = "failed_prediction/" + pdb1_name
                if not os.path.exists(gen_dir):
                    os.mkdir(gen_dir)

                mv_folder_cmd = (
                    "mv "
                    + pdb1_name
                    + "_predicted_models_full_rand_"
                    + str(random_seed_full_MSA)
                    + " failed_prediction/"
                    + pdb1_name
                )
                print(mv_folder_cmd)
                os.system(mv_folder_cmd)
                print("All predictions with deep MSA are failed")
                sys.exit()
            else:
                if np.average(full_TMscore[0, :]) > np.average(full_TMscore[1, :]):
                    ref_name = pdb1_name
                    alt_name = pdb2_name
                else:
                    ref_name = pdb2_name
                    alt_name = pdb1_name

            print("Reference structure: ", ref_name)
            print("Alternative structure: ", alt_name)

            # save TM-score from full-length MSA
            np.savetxt("TMScore_full-MSA_" + pdb1_name + ".csv", full_TMscore, fmt="%2.3f")
            print("Full-MSA prediction is tightly aligned to crystal structure")
            print("               ")

            # Evaluate TM-scores from shallow/random MSAs
            max_msa = 1
            ext_msa = 2
            TMscores_random = []

            for multi in (1, 2, 2, 2, 2, 2, 2):
                max_msa = max_msa * multi
                ext_msa = ext_msa * multi

                pred_dir = (
                    "predictions_all/"
                    + pdb1_name
                    + "/"
                    + pdb1_name
                    + "_predicted_models_rand_*"
                    + "_max_"
                    + str(max_msa)
                    + "_ext_"
                    + str(ext_msa)
                    + "/"
                )
                print(pred_dir)
                # TM-score for whole structure
                MSA_shallow_TMscore = TMScore(
                    pred_dir, pdb1, pdb1_name, pdb2, pdb2_name, model_type
                )
                TMscores_random = np.append(TMscores_random, MSA_shallow_TMscore.tmscores)

            fin_pred_dir = (
                "predictions_all/"
                + pdb1_name
                + "/"
                + pdb1_name
                + "_predicted_models_rand_*"
                + "_max_*"
            )
            TMscores_random_reshape = TMscores_random.reshape(14, num_seeds * 5)
            TMscores_random_alter = np.zeros((7, num_seeds * 5))

            # Extract TM-scores for the alternative conformation
            if alt_name == pdb2_name:
                # for i in 1, 3, 5, 7, 9, 11, 13 in TM_scores:
                tmp_cnt = 0
                for i in range(1, 14, 2):
                    print(TMscores_random_reshape[i, :])
                    TMscores_random_alter[tmp_cnt, :] = TMscores_random_reshape[i, :]
                    tmp_cnt = tmp_cnt + 1
            else:
                # for i in 0, 2, 4, 6, 8, 10, 12 in TM_scores:
                tmp_cnt = 0
                for i in range(0, 13, 2):
                    print(TMscores_random_reshape[i, :])
                    TMscores_random_alter[tmp_cnt, :] = TMscores_random_reshape[i, :]
                    tmp_cnt = tmp_cnt + 1

            # Evaluate varied MSA sizes to find an optimal shallow MSA
            if np.all(TMscores_random_alter) < 0.5:
                print("All predictions are failed")
                sys.exit()

            else:
                print()
                print("Finding optimal size of random MSA...")
                MSA_shallow_TMscore.select_size(
                    TMscores_random_reshape, pdb1_name, pdb2_name, alt_name, num_seeds
                )

                size_selection = MSA_shallow_TMscore.selection
                print(size_selection)
                self.size_selection = size_selection
                ## save all TM-scores from random MSA (1-2, 2-4, 4-8.... in order)
                np.savetxt(
                    "TMScore_random-MSA_" + pdb1_name + ".csv", TMscores_random_reshape, fmt="%2.3f"
                )

        elif model_type == "alphafold2_multimer_v3":
            print("Currently working on")
            MSA_multi = PredictionAll_multimer(
                pdb1_name, pdb2_name, search_dir, nMSA, model_type, search_multi_dir
            )
            self.size_selection = MSA_multi.size_selection
