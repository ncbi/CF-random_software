import argparse
import glob
import os
import sys
import warnings
from pathlib import (
    Path,
)

import numpy as np

warnings.filterwarnings("ignore")

from colabfold.download import (
    download_alphafold_params,
)

from ..analysis.cal_plddt_ACFS import (
    PlddtCal,
)
from ..analysis.tmscore_all_var import (
    TMScoreCalAllVar,
)
from ..analysis.tmscore_all_var_fs import (
    TMScoreCalAllVarFS,
)
from ..plotting.plot_ac import (
    Plot2DScatterAC,
)
from ..plotting.plot_fc import (
    Plot2DScatter,
)
from ..prediction.prediction_all_var import (
    PredictionAll,
)
from ..utils.search_foldseek_cluster import (
    BlindScreening,
)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--pdb1",
        type=str,
        help="PDB structure for the target crystal structure (target to be predicted)",
    )
    parser.add_argument(
        "--pdb2",
        type=str,
        help="PDB structure for the alternative crystal structure",
    )
    parser.add_argument("--fname", type=str, help="MSA folder name after colabsearch")
    parser.add_argument("--fmname", type=str, help="multimer MSA folder name after colabsearch")
    parser.add_argument("--pname", type=str, help="job name for predicting blind mode")
    parser.add_argument(
        "--nMSA",
        type=int,
        default=0,
        help="number of additional MSA seeds to run (added to default 5)",
    )
    parser.add_argument(
        "--nENS",
        type=int,
        default=0,
        help="number of ensemble samples to generate (integer)",
    )
    parser.add_argument(
        "--option",
        type=str,
        help=(
            "select prediction mode: AC (alternative conformation), "
            "FS (fold-switching), inAC (increased sampling for alternative conformation), "
            "or blind (no crystal structures)"
        ),
    )
    parser.add_argument(
        "--type",
        type=str,
        help="select model-type of Colabfold: ptm, monomer, or multimer",
    )
    args = parser.parse_args()

    download_alphafold_params("alphafold2_ptm", Path("."))

    blind = "predictions_all"
    success = "predictions_all"
    fail = "failed_prediction"
    multi = "multimer_prediction"
    pwd = os.getcwd() + "/"

    pdb1 = None
    pdb2 = None
    pdb1_name = None
    pdb2_name = None

    if args.option == "blind":
        if args.pname is not None:
            pdb1_name = args.pname
        elif args.fname is not None:
            pdb1_name = args.fname.replace("/", "")
        else:
            print("Error: blind mode requires --pname or --fname")
            sys.exit(1)
        print("work name:", pdb1_name)
    elif args.pdb1 is not None and args.pdb2 is not None:
        pdb1 = args.pdb1
        pdb2 = args.pdb2
        pdb1_name = pdb1.replace(".pdb", "")
        pdb2_name = pdb2.replace(".pdb", "")
        print(pdb1_name, pdb2_name)
    else:
        print("Error: non-blind modes require --pdb1 and --pdb2")
        sys.exit(1)

    nMSA = args.nMSA
    nENS = args.nENS

    # --- Resolve search directories ---
    if args.fname is None:
        print("Error: --fname (MSA folder) is required for all modes")
        sys.exit(1)

    search_dir = args.fname
    search_multi_dir = args.fmname if args.fmname is not None else None

    if args.type is None or args.type == "ptm":
        model_type = "alphafold2_ptm"
    elif args.type == "monomer":
        model_type = "alphafold2"
    elif args.type == "multimer":
        model_type = "alphafold2_multimer_v3"
        target_dir = multi if args.option == "blind" else multi
        if not os.path.exists(target_dir):
            os.mkdir(target_dir)
        if args.option != "blind" and pdb1 is not None:
            TER_count = 0
            with open(pdb1, "r") as f:
                for line in f:
                    TER_count += line.split().count("TER")
            print(TER_count, "chain(s) in this multimer file.")
    else:
        print("Error: unrecognised --type option. Choose from: ptm, monomer, multimer")
        sys.exit(1)

    success = "predictions_all/" + pdb1_name + "/"

    if args.option == "AC":
        print("Predicting alternative conformations")

        succ_dir_count = 0
        if not os.path.exists(success):
            os.mkdir(success)
        else:
            for _root, cur_dir, _files in os.walk(pwd + success):
                succ_dir_count += len(cur_dir)

        if os.path.exists(success) and succ_dir_count > 0 and succ_dir_count < 8:
            print("Folder exists but is incomplete — cleaning subfolders")
            os.system("rm -rf " + success)

        prediction_option = args.option
        if os.path.exists(success) and succ_dir_count >= 8:
            print("Predictions including full- and random-MSA were already done")
            cal_TMscore = TMScoreCalAllVar(
                pdb1, pdb1_name, pdb2, pdb2_name, nMSA, prediction_option, model_type
            )
        else:
            PredictionAll(pdb1_name, search_dir, search_multi_dir, nMSA, model_type)
            cal_TMscore = TMScoreCalAllVar(
                pdb1, pdb1_name, pdb2, pdb2_name, nMSA, prediction_option, model_type
            )

        shallow_MSA_size = np.append([], cal_TMscore.size_selection)
        print("               ")
        print("Specific size of shallow random MSA is similar to full-MSA")
        print(shallow_MSA_size)
        np.savetxt("selected_MSA-size_" + pdb1_name + ".csv", shallow_MSA_size)

        if model_type == "alphafold2_multimer_v3":
            base = pwd + multi + "/" + pdb1_name
            list_org_samplings = glob.glob(base + "/*full_rand*/")
            list_ran_samplings = glob.glob(base + "/*max*/")
        else:
            list_org_samplings = glob.glob(pwd + success + "*full_rand*/")
            list_ran_samplings = glob.glob(pwd + success + "*max*/")

        full = "full-MSA"
        random = "random-MSA"
        PlddtCal(list_org_samplings, full, pdb1_name, nMSA, nENS, model_type)
        PlddtCal(list_ran_samplings, random, pdb1_name, nMSA, nENS, model_type)
        Plot2DScatterAC(full, random, pdb1, pdb1_name, pdb2, pdb2_name, nMSA, nENS, model_type)

    elif args.option == "FS":
        print("Predicting fold-switching models")

        succ_dir_count = 0
        if not os.path.exists(success):
            os.mkdir(success)
        else:
            for _root, cur_dir, _files in os.walk(pwd + success):
                succ_dir_count += len(cur_dir)

        if os.path.exists(success) and succ_dir_count > 0 and succ_dir_count < 8:
            print("Folder exists but is incomplete — cleaning subfolders")
            os.system("rm -rf " + success)

        prediction_option = args.option
        if os.path.exists(success) and succ_dir_count >= 8:
            print("Predictions including full- and random-MSA were already done")
            cal_TMscore = TMScoreCalAllVarFS(
                pdb1, pdb1_name, pdb2, pdb2_name, nMSA, prediction_option, model_type
            )
            shallow_MSA_size = np.append([], cal_TMscore.size_selection)
        else:
            PredictionAll(pdb1_name, search_dir, search_multi_dir, nMSA, model_type)
            if args.type != "multimer":
                cal_TMscore = TMScoreCalAllVarFS(
                    pdb1, pdb1_name, pdb2, pdb2_name, nMSA, prediction_option, model_type
                )
                shallow_MSA_size = np.append([], cal_TMscore.size_selection)
            else:
                shallow_MSA_size = np.array([])

        print("               ")
        print("Specific size of shallow random MSA is similar to full-MSA")
        print(shallow_MSA_size)
        np.savetxt("selected_MSA-size_" + pdb1_name + ".csv", shallow_MSA_size)

        if model_type == "alphafold2_multimer_v3":
            base = pwd + multi + "/" + pdb1_name
            list_org_samplings = glob.glob(base + "/*full_rand*/")
            list_ran_samplings = glob.glob(base + "/*max*/")
        else:
            list_org_samplings = glob.glob(pwd + success + "*full_rand*/")
            list_ran_samplings = glob.glob(pwd + success + "*max*/")

        full = "full-MSA"
        random = "random-MSA"
        PlddtCal(list_org_samplings, full, pdb1_name, nMSA, nENS, model_type)
        PlddtCal(list_ran_samplings, random, pdb1_name, nMSA, nENS, model_type)

        if model_type == "alphafold2_multimer_v3":
            Plot2DScatterAC(full, random, pdb1, pdb1_name, pdb2, pdb2_name, nMSA, nENS, model_type)
        else:
            Plot2DScatter(full, random, pdb1, pdb1_name, pdb2, pdb2_name, nMSA, nENS)

    elif args.option == "blind":
        print("Predicting fold-switching proteins without crystal structures")

        if not os.path.exists(blind):
            os.mkdir(blind)

        blind_pdb_dir = blind + "/" + pdb1_name
        blind_dir_count = 0
        if os.path.exists(blind_pdb_dir):
            for _root, cur_dir, _files in os.walk(pwd + blind_pdb_dir + "/"):
                blind_dir_count += len(cur_dir)
            if blind_dir_count > 0 and blind_dir_count < 8:
                print("Folder exists but is incomplete — cleaning subfolders")
                os.system("rm -rf " + blind_pdb_dir + "/")

        blind_pred_path = "predictions_all/" + pdb1_name
        print(blind_pred_path)

        if os.path.exists(blind_pdb_dir) and blind_dir_count >= 8:
            print("Predictions including full- and random-MSA were already done")
            fseek_file_count = 0
            for _root, _dirs, files in os.walk(pwd + blind_pdb_dir + "/"):
                fseek_file_count += len(files)
            print(fseek_file_count)
            BlindScreening(pdb1_name, blind_pred_path)
        else:
            PredictionAll(pdb1_name, search_dir, search_multi_dir, nMSA, model_type)
            print("               ")
            print("Finished running predictions using full- and shallow random-MSAs")
            print("               ")
            print("Running Foldseek to find related crystal structures")
            BlindScreening(pdb1_name, blind_pred_path)

    else:
        print("Error: unrecognised --option. Choose from: AC, FS, blind, inAC")
        sys.exit(1)


if __name__ == "__main__":
    main()
