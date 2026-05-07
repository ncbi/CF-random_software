#!/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 21 14:51:00 2024

@author: Myeongsang (Samuel) Lee
"""

import glob
import os
import random
import sys
from typing import List, Tuple, Union

import numpy as np

from ..utils.convert_multi_single import convert_m2s

from .pred_cal_tmscore_multimer import CF_MSA_max, CF_MSA_var


class TM_score:
    """Compute TM-scores of predicted models against reference PDBs.

    The class finds predicted model files in `pred_dir`, optionally converts
    multimer outputs to single-chain format, and computes TM-align scores
    against `pdb1` and `pdb2`.
    """

    def __init__(
        self,
        pred_dir: str,
        pdb1: str,
        pdb1_name: str,
        pdb2: str,
        pdb2_name: str,
        model_type: str,
    ) -> None:
        """Initialize and compute TM-scores.

        Args:
            pred_dir: Path to predicted model directory.
            pdb1: Path or id for reference structure 1.
            pdb1_name: Name/ID for reference 1 (used for multimer conversion).
            pdb2: Path or id for reference structure 2.
            pdb2_name: Name/ID for reference 2.
            model_type: Model type string; affects multimer handling.
        """
        import glob
        from tmtools import tm_align
        from tmtools.io import get_residue_data, get_structure
        from tmtools.testing import get_pdb_path

        pwd = os.getcwd() + "/"

        # Gather predicted model files
        if model_type != "alphafold2_multimer_v3":
            files_list = glob.glob(str(pred_dir) + "/*_unrelaxed*pdb")
        else:
            check_files_list = glob.glob(str(pred_dir) + "/rmTER*_unrelaxed*pdb")
            if not check_files_list:
                convert_m2s(pred_dir, pdb1_name, pdb2_name)
            files_list = glob.glob(str(pred_dir) + "/rmTER*_unrelaxed*pdb")

        if len(files_list) == 0:
            # No predicted models found; store a default score list
            self.tmscores = [0.0, 0.0, 0.0, 0.0, 0.0]
            return

        # Load reference structures
        pdb1_dir = pwd + pdb1_name
        r2 = get_structure(get_pdb_path(str(pdb1_dir)))
        coords2_ref, seq2_ref = get_residue_data(r2)

        pdb2_dir = pwd + pdb2_name
        r3 = get_structure(get_pdb_path(str(pdb2_dir)))
        coords3_ref, seq3_ref = get_residue_data(r3)

        # Helper to compute scores against a reference
        def _scores_against(ref_coords, ref_seq) -> Tuple[List[float], List[float]]:
            ord_list: List[float] = []
            rev_list: List[float] = []
            for model in files_list:
                m = model.replace(".pdb", "")
                modelpath = pwd + m
                s = get_structure(get_pdb_path(modelpath))
                coords1, seq1 = get_residue_data(s)
                res = tm_align(coords1, ref_coords, seq1, ref_seq)
                ord_list.append(round(res.tm_norm_chain1, 5))

                res = tm_align(ref_coords, coords1, ref_seq, seq1)
                rev_list.append(round(res.tm_norm_chain1, 5))
            return ord_list, rev_list

        tmscores_ord1, tmscores_rev1 = _scores_against(coords2_ref, seq2_ref)
        tmscores_ord2, tmscores_rev2 = _scores_against(coords3_ref, seq3_ref)

        # Combine results and select orientation with higher max
        tmscores_ord = tmscores_ord1 + tmscores_ord2
        tmscores_rev = tmscores_rev1 + tmscores_rev2

        if np.max(tmscores_ord) > np.max(tmscores_rev):
            tmscores = tmscores_ord
        else:
            tmscores = tmscores_rev

        self.tmscores = tmscores


class CF_MSA_max:
    """Run ColabFold batch with maximum MSA settings.

    This is a thin wrapper that constructs and executes the `colabfold_batch`
    command for deep/full MSA generation.
    """

    def __init__(
        self,
        search_dir: str,
        output_dir: str,
        pdb_name: str,
        rseed: Union[int, str],
        num_seeds: int,
        model_type: str,
    ) -> None:
        command = (
            "colabfold_batch --num-seeds "
            + str(num_seeds)
            + " --model-type "
            + str(model_type)
            + " --random-seed "
            + str(rseed)
            + search_dir
            + output_dir
        )
        print(command)
        os.system(command)


class CF_MSA_var:
    """Run ColabFold with varying (shallow) MSA sizes.

    This wrapper launches a series of `colabfold_batch` runs with increasing
    MSA depth and extra-sequence parameters.
    """

    def __init__(
        self,
        pdb1: str,
        pdb1_name: str,
        pdb2: str,
        pdb2_name: str,
        search_dir: str,
        output_dir: str,
        rseed: Union[int, str],
        num_seeds: int,
        model_type: str,
    ) -> None:
        # shallow MSA section
        max_msa = 1
        ext_msa = 2
        # keep original seed representation
        random_seed = np.array(rseed)

        self.pdb1_name = pdb1_name

        for multi in (1, 2, 2, 2, 2, 2, 2):
            max_msa = max_msa * multi
            ext_msa = ext_msa * multi

            # Build and execute ColabFold command for this MSA size
            command = (
                "colabfold_batch --num-seeds "
                + str(num_seeds)
                + " --model-type "
                + str(model_type)
                + " --max-seq "
                + str(max_msa)
                + " --max-extra-seq "
                + str(ext_msa)
                + search_dir
                + output_dir
                + str(rseed)
                + "_max_"
                + str(max_msa)
                + "_ext_"
                + str(ext_msa)
            )
            print(command)
            os.system(command)

    def select_size(
        self,
        TMscores_random_alter: np.ndarray,
        pdb1_name: str,
        pdb2_name: str,
        alt_name: str,
        num_seeds: int,
    ) -> None:
        """Select the optimal MSA size based on TM-scores.

        Args:
            TMscores_random_alter: Flattened TM-score array for varied MSA runs.
            pdb1_name: Name of reference structure 1.
            pdb2_name: Name of reference structure 2.
            alt_name: Name of the alternative structure chosen earlier.
            num_seeds: Number of seeds used per run.
        """
        TMscores_random_reshape = TMscores_random_alter.reshape(14, num_seeds * 5)
        TMscores_random_locat = np.zeros((7, num_seeds * 5))

        # Extract location-specific rows depending on which alt_name was chosen
        if alt_name == pdb2_name:
            rows = list(range(1, 14, 2))
        else:
            rows = list(range(0, 14, 2))

        for tmp_cnt, i in enumerate(rows):
            TMscores_random_locat[tmp_cnt, :] = TMscores_random_reshape[i, :]

        TMscore_data = TMscores_random_locat.reshape(7, num_seeds * 5)
        TMscore_data_sum = TMscore_data.sum(axis=1)

        location = int(np.argmax(TMscore_data_sum))

        print("Selecting...")

        TMscore_data_full = TMscores_random_alter.reshape(14, num_seeds * 5)
        if alt_name == pdb2_name:
            location = (location * 2) + 1
        else:
            location = location * 2

        if alt_name == pdb2_name and np.any(TMscore_data_full[location, :] >= 0.5):
            selection = int((location - 1) / 2)
            self.selection = selection
        elif alt_name == pdb1_name and np.any(TMscore_data_full[location, :] >= 0.5):
            selection = int(location / 2)
            self.selection = selection
        else:
            print("Predictions are bad")
            print("Predictions of whole structure are bad")
            rm_folder_cmd = f"rm -rf successed_prediction/{self.pdb1_name}/"
            print(rm_folder_cmd)
            os.system(rm_folder_cmd)
            sys.exit()


class prediction_all_AC:
    """Orchestrate the full prediction workflow for AC (alternate conformation) analysis.

    This class runs deep (full-MSA) and varied shallow MSA ColabFold predictions,
    evaluates TM-scores, selects a reference/alternative structure, and
    determines the optimal shallow-MSA size.
    """

    def __init__(
        self,
        pdb1: str,
        pdb1_name: str,
        pdb2: str,
        pdb2_name: str,
        search_dir: str,
        nMSA: int,
        model_type: str,
        search_multi_dir: str,
    ) -> None:
        num_seeds = 5 + nMSA

        if model_type != "alphafold2_multimer_v3":

            ##### Perform prediction with full-length MSA
            pre_random_seed = np.random.randint(0, 16, 1)
            random_seed_full_MSA = "".join(map(str, pre_random_seed))
            output_dir = (
                " " + pdb1_name + "_predicted_models_full_rand_" + str(random_seed_full_MSA)
            )
            MSA_full = CF_MSA_max(
                search_dir, output_dir, pdb1_name, random_seed_full_MSA, num_seeds, model_type
            )

            ##### Perform prediction with random shallow MSA
            ##### check out varied-MSA with (msa-max: 1, 2, 4, 8, 16, 32, 64) (msa-extra: 2, 4, 8, 16, 32, 64, 128)
            output_dir = " " + pdb1_name + "_predicted_models_rand_"
            random_seed = random.sample(range(100), 1)
            random_seed = "".join(map(str, random_seed))
            MSA_var = CF_MSA_var(
                pdb1,
                pdb1_name,
                pdb2,
                pdb2_name,
                search_dir,
                output_dir,
                random_seed,
                num_seeds,
                model_type,
            )

            ####################################################################
            ##### check-out TM-scores of prediction with full-length-MSA (whole)
            pred_dir = pdb1_name + "_predicted_models_full_rand_" + str(random_seed_full_MSA) + "/"
            print(pred_dir)
            MSA_full_TMscore = TM_score(pred_dir, pdb1, pdb1_name, pdb2, pdb2_name, model_type)
            full_TMscore = np.array(MSA_full_TMscore.tmscores)
            full_TMscore = full_TMscore.reshape(2, num_seeds * 5)

            ##### check-out the 1st prediction results are good or not
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

            # Directory section and save to successed_prediction folder
            gen_dir = "successed_prediction/" + pdb1_name

            if not os.path.exists(gen_dir):
                os.mkdir(gen_dir)

            mv_folder_cmd = "mv " + pred_dir + " successed_prediction/" + pdb1_name
            print(mv_folder_cmd)
            os.system(mv_folder_cmd)
            print("Full-MSA prediction is tightly aligned to crystal structure")
            print("               ")

            ################################################################
            ##### chech-out TM-scores of prediction with shallow random MSAs
            max_msa = 1
            ext_msa = 2
            TMscores_random = []

            for multi in (1, 2, 2, 2, 2, 2, 2):
                max_msa = max_msa * multi
                ext_msa = ext_msa * multi

                pred_dir = (
                    pdb1_name
                    + "_predicted_models_rand_"
                    + str(random_seed)
                    + "_max_"
                    + str(max_msa)
                    + "_ext_"
                    + str(ext_msa)
                    + "/"
                )
                print(pred_dir)
                MSA_shallow_TMscore = TM_score(
                    pred_dir, pdb1, pdb1_name, pdb2, pdb2_name, model_type
                )
                TMscores_random = np.append(TMscores_random, MSA_shallow_TMscore.tmscores)

            fin_pred_dir = pdb1_name + "_predicted_models_rand_" + str(random_seed) + "_max_*"
            TMscores_random_reshape = TMscores_random.reshape(14, num_seeds * 5)
            TMscores_random_alter = np.zeros((7, num_seeds * 5))

            #### finding alternative pdb_name
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

            ##### check out varied-MSA with (msa-max: 1, 2, 4, 8, 16, 32, 64) (msa-extra: 2, 4, 8, 16, 32, 64, 128)
            if np.all(TMscores_random_alter) < 0.5:
                print("All predictions are failed")
                mv_command = "mv " + fin_pred_dir + " failed_prediction/" + pdb1_name
                print(mv_command)
                os.system(mv_command)
                sys.exit()

            else:
                print("     ")
                print("Finding optimal size of ramdon MSA...")
                MSA_var.select_size(
                    TMscores_random_reshape, pdb1_name, pdb2_name, alt_name, num_seeds
                )

                size_selection = MSA_var.selection
                print(size_selection)
                self.size_selection = size_selection
                ## save all TM-scores from random MSA (1-2, 2-4, 4-8.... in order)
                np.savetxt(
                    "TMScore_random-MSA_" + pdb1_name + ".csv", TMscores_random_reshape, fmt="%2.3f"
                )

                mv_command = "mv " + fin_pred_dir + " successed_prediction/" + pdb1_name
                print(mv_command)
                os.system(mv_command)

        elif model_type == "alphafold2_multimer_v3":
            print("Currently working on")
            MSA_multi = prediction_all_multimer(
                pdb1_name, pdb2_name, search_dir, nMSA, model_type, search_multi_dir
            )
            self.size_selection = MSA_multi.size_selection
