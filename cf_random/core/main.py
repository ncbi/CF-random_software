#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Main entry point for CF-random protein structure prediction pipeline.

This module orchestrates the prediction and analysis workflow for alternative
conformation (AC), fold-switching (FS), and blind mode predictions using
ColabFold and AlphaFold models.
"""

import argparse
import glob
import logging
import os
import shutil
import warnings
from pathlib import (
    Path,
)
from typing import (
    Optional,
    Union,
)

import numpy as np
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
from ..plotting.plot_fs import (
    Plot2DScatter,
)
from ..prediction.prediction_all_var import (
    PredictionAll,
)
from ..utils.search_foldseek_cluster import (
    BlindScreening,
)

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)

# Suppress warnings from dependencies
warnings.filterwarnings("ignore")

# Constants
MODEL_TYPES = {
    "ptm": "alphafold2_ptm",
    "monomer": "alphafold2",
    "multimer": "alphafold2_multimer_v3",
}
MULTI = "multimer_prediction"
BLIND = "predictions_all"
SUCCESS = "predictions_all"

FOLDSEEK_DONE_FILE_COUNT = 640


def parse_arguments() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="CF-random protein structure prediction pipeline",
        formatter_class=argparse.RawTextHelpFormatter,
    )
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
        type=str,
        help="number of additional MSA seeds to run (added to default 5)",
    )
    parser.add_argument(
        "--nENS",
        type=str,
        help="number of ensemble samples to generate (integer)",
    )
    parser.add_argument(
        "--option",
        type=str,
        required=True,
        help=(
            "select prediction mode: AC (alternative conformation), "
            "FS (fold-switching), inAC (increased AC sampling), "
            "or blind (no crystal structures)"
        ),
    )
    parser.add_argument(
        "--type",
        type=str,
        choices=["ptm", "monomer", "multimer"],
        help="select model-type of Colabfold: ptm, monomer, or multimer",
    )
    return parser.parse_args()


def resolve_pdb1_name(args: argparse.Namespace) -> str:
    """Resolve the working name for blind mode."""
    if args.pdb1 is None and args.pdb2 is None:
        return args.pname
    elif args.pdb1 is None and args.pname is None:
        return args.fname.replace("/", "")
    else:
        return args.fname.replace("/", "")


def resolve_nMSA_nENS(args: argparse.Namespace):
    """Resolve nMSA and nENS from optional string arguments."""
    nMSA_raw = args.nMSA
    nENS_raw = args.nENS

    if nMSA_raw is None and nENS_raw is None:
        return 0, 0
    elif nMSA_raw is not None and nENS_raw is not None:
        return int(nMSA_raw), int(nENS_raw)
    elif nMSA_raw is None and nENS_raw is not None:
        return 0, int(nENS_raw)
    elif nMSA_raw is not None and nENS_raw is None:
        return int(nMSA_raw), 0
    else:
        raise ValueError("Please provide a valid combination of --nMSA and --nENS")


def resolve_search_dirs(args: argparse.Namespace):
    """Resolve search_dir and search_multi_dir."""
    if args.fname is None and args.fmname is None:
        raise ValueError("--fname (MSA folder) is required for all modes")
    elif args.fname is None and args.fmname is not None:
        raise ValueError("--fname (monomer MSA folder) is required alongside --fmname")
    elif args.fname is not None and args.fmname is None:
        return args.fname, 0
    else:
        return args.fname, " " + args.fmname


def determine_model_type(args: argparse.Namespace, pdb1: Optional[str]) -> str:
    """Determine the model type and set up multimer directory if needed."""
    if args.type is None or args.type == "ptm":
        return MODEL_TYPES["ptm"]
    elif args.type == "monomer":
        return MODEL_TYPES["monomer"]
    elif args.type == "multimer" and args.option == "blind":
        model_type = MODEL_TYPES["multimer"]
        Path(MULTI).mkdir(parents=True, exist_ok=True)
        return model_type
    elif args.type == "multimer":
        ter_count = count_chains(pdb1)
        logger.info("%d chain(s) in this multimer file.", ter_count)
        model_type = MODEL_TYPES["multimer"]
        Path(MULTI).mkdir(parents=True, exist_ok=True)
        return model_type
    else:
        raise ValueError(
            f"Unrecognized model type: {args.type!r}. Choose from: ptm, monomer, multimer"
        )


def count_chains(pdb_file: str) -> int:
    """Count the number of chains in a PDB file."""
    ter_count = 0
    with open(pdb_file, "r") as f:
        for line in f:
            ter_count += line.split().count("TER")
    return ter_count


def main() -> None:
    """Main entry point for the CF-random pipeline."""
    args = parse_arguments()

    download_alphafold_params("alphafold2_ptm", Path("."))

    pwd = os.getcwd() + "/"
    pdb1: Optional[str] = None
    pdb2: Optional[str] = None
    pdb1_name: Optional[str] = None
    pdb2_name: Optional[str] = None

    # Resolve working names
    if args.option == "blind":
        pdb1_name = resolve_pdb1_name(args)
        logger.info("Work name: %s", pdb1_name)
    elif args.pdb1 is not None and args.pdb2 is not None:
        pdb1 = args.pdb1
        pdb2 = args.pdb2
        pdb1_name = pdb1.replace(".pdb", "")
        pdb2_name = pdb2.replace(".pdb", "")
        logger.info("PDB names: %s, %s", pdb1_name, pdb2_name)

    nMSA, nENS = resolve_nMSA_nENS(args)
    search_dir, search_multi_dir = resolve_search_dirs(args)
    model_type = determine_model_type(args, pdb1)

    search_dir = args.fname
    success_dir = f"{SUCCESS}/{pdb1_name}/"

    if args.option == "AC":
        run_alternative_conformation_workflow(
            pdb1,
            pdb1_name,
            pdb2,
            pdb2_name,
            nMSA,
            nENS,
            model_type,
            search_dir,
            search_multi_dir,
            success_dir,
            pwd,
        )
    elif args.option == "FS":
        run_fold_switching_workflow(
            pdb1,
            pdb1_name,
            pdb2,
            pdb2_name,
            nMSA,
            nENS,
            model_type,
            search_dir,
            search_multi_dir,
            success_dir,
            pwd,
        )
    elif args.option == "blind":
        run_blind_workflow(pdb1_name, search_dir, search_multi_dir, nMSA, model_type)
    else:
        raise ValueError(f"Unrecognized option: {args.option!r}. Choose from: AC, FS, inAC, blind")


def run_alternative_conformation_workflow(
    pdb1: str,
    pdb1_name: str,
    pdb2: str,
    pdb2_name: str,
    nMSA: int,
    nENS: int,
    model_type: str,
    search_dir: str,
    search_multi_dir: Union[int, str],
    success_dir: str,
    pwd: str,
) -> None:
    """Run the alternative conformation prediction workflow."""
    logger.info("Predicting alternative conformations")

    if not os.path.exists(success_dir):
        Path(success_dir).mkdir(parents=True, exist_ok=True)
        succ_dir_count = 0
    else:
        succ_dir_count = 0
        for root_dir, cur_dir, files in os.walk(pwd + success_dir + "/"):
            succ_dir_count += len(cur_dir)

    if os.path.exists(success_dir) and 0 < succ_dir_count < 8:
        logger.info("Folder exists but is incomplete — cleaning subfolders")
        shutil.rmtree(success_dir)

    if os.path.exists(success_dir) and succ_dir_count >= 8:
        logger.info("Predictions including full and random-MSA were already completed.")
        calculate_tm_score = TMScoreCalAllVar(
            pdb1, pdb1_name, pdb2, pdb2_name, nMSA, "AC", model_type, search_dir, search_multi_dir
        )
    else:
        PredictionAll(pdb1_name, search_dir, search_multi_dir, nMSA, model_type)
        calculate_tm_score = TMScoreCalAllVar(
            pdb1, pdb1_name, pdb2, pdb2_name, nMSA, "AC", model_type, search_dir, search_multi_dir
        )

    shallow_MSA_size = np.append([], calculate_tm_score.size_selection)
    logger.info("Specific size of shallow random MSA is similar to full-MSA: %s", shallow_MSA_size)
    np.savetxt("selected_MSA-size_" + pdb1_name + ".csv", shallow_MSA_size)

    if model_type == MODEL_TYPES["multimer"]:
        base = pwd + MULTI + "/" + pdb1_name
        list_org_samplings = glob.glob(base + "/*full_rand*/")
        list_ran_samplings = glob.glob(base + "/*max*/")
    else:
        list_org_samplings = glob.glob(pwd + success_dir + "*full_rand*/")
        list_ran_samplings = glob.glob(pwd + success_dir + "*max*/")

    full = "full-MSA"
    random = "random-MSA"
    PlddtCal(list_org_samplings, full, pdb1_name, nMSA, nENS, model_type)
    PlddtCal(list_ran_samplings, random, pdb1_name, nMSA, nENS, model_type)
    Plot2DScatterAC(full, random, pdb1, pdb1_name, pdb2, pdb2_name, nMSA, nENS, model_type)


def run_fold_switching_workflow(
    pdb1: str,
    pdb1_name: str,
    pdb2: str,
    pdb2_name: str,
    nMSA: int,
    nENS: int,
    model_type: str,
    search_dir: str,
    search_multi_dir: Union[int, str],
    success_dir: str,
    pwd: str,
) -> None:
    """Run the fold-switching prediction workflow."""
    logger.info("Predicting fold-switching models.")

    if not os.path.exists(success_dir):
        Path(success_dir).mkdir(parents=True, exist_ok=True)
        succ_dir_count = 0
    else:
        succ_dir_count = 0
        for root_dir, cur_dir, files in os.walk(pwd + success_dir + "/"):
            succ_dir_count += len(cur_dir)

    if os.path.exists(success_dir) and 0 < succ_dir_count < 8:
        logger.info("Folder exists but is incomplete — cleaning subfolders")
        shutil.rmtree(success_dir)

    shallow_MSA_size = np.array([])

    if os.path.exists(success_dir) and succ_dir_count >= 8:
        logger.info("Predictions including full and random-MSA were already completed.")
        calculate_tm_score = TMScoreCalAllVarFS(
            pdb1, pdb1_name, pdb2, pdb2_name, nMSA, "FS", model_type, search_dir, search_multi_dir
        )
        shallow_MSA_size = np.append(shallow_MSA_size, calculate_tm_score.size_selection)
    else:
        PredictionAll(pdb1_name, search_dir, search_multi_dir, nMSA, model_type)
        if model_type != MODEL_TYPES["multimer"]:
            calculate_tm_score = TMScoreCalAllVarFS(
                pdb1,
                pdb1_name,
                pdb2,
                pdb2_name,
                nMSA,
                "FS",
                model_type,
                search_dir,
                search_multi_dir,
            )
            shallow_MSA_size = np.append(shallow_MSA_size, calculate_tm_score.size_selection)

    logger.info("Specific size of shallow random MSA is similar to full-MSA: %s", shallow_MSA_size)
    np.savetxt("selected_MSA-size_" + pdb1_name + ".csv", shallow_MSA_size)

    if model_type == MODEL_TYPES["multimer"]:
        base = pwd + MULTI + "/" + pdb1_name
        list_org_samplings = glob.glob(base + "/*full_rand*/")
        list_ran_samplings = glob.glob(base + "/*max*/")
    else:
        list_org_samplings = glob.glob(pwd + success_dir + "/*full_rand*/")
        list_ran_samplings = glob.glob(pwd + success_dir + "/*max*/")

    full = "full-MSA"
    random = "random-MSA"
    PlddtCal(list_org_samplings, full, pdb1_name, nMSA, nENS, model_type)
    PlddtCal(list_ran_samplings, random, pdb1_name, nMSA, nENS, model_type)

    if model_type == MODEL_TYPES["multimer"]:
        Plot2DScatterAC(full, random, pdb1, pdb1_name, pdb2, pdb2_name, nMSA, nENS, model_type)
    else:
        Plot2DScatter(full, random, pdb1, pdb1_name, pdb2, pdb2_name, nMSA, nENS)


def run_blind_workflow(
    pdb1_name: str,
    search_dir: str,
    search_multi_dir: Union[int, str],
    nMSA: int,
    model_type: str,
) -> None:
    """Run the blind prediction workflow."""
    logger.info("Predicting fold-switching proteins without crystal structures.")

    if not os.path.exists(BLIND):
        Path(BLIND).mkdir(parents=True, exist_ok=True)

    blind_pdb_dir = BLIND + "/" + pdb1_name
    blind_pred_path = f"{BLIND}/{pdb1_name}"
    logger.info("Blind prediction path: %s", blind_pred_path)

    blind_dir_count = 0
    if os.path.exists(blind_pdb_dir):
        for root_dir, cur_dir, files in os.walk(blind_pdb_dir + "/"):
            blind_dir_count += len(cur_dir)

    if os.path.exists(blind_pdb_dir) and 0 < blind_dir_count < 8:
        logger.info("Folder exists but is incomplete — cleaning subfolders")
        shutil.rmtree(blind_pdb_dir)

    if os.path.exists(blind_pdb_dir) and blind_dir_count >= 8:
        logger.info("Predictions including full and random-MSA were already completed.")

        # Count total files to determine whether Foldseek has already been run
        fseek_file_count = 0
        for root_dir, cur_dir, files in os.walk(blind_pdb_dir + "/"):
            fseek_file_count += len(files)
        logger.info("Number of files in blind prediction path: %d", fseek_file_count)

        if fseek_file_count >= FOLDSEEK_DONE_FILE_COUNT:
            logger.info("Foldseek search was already done")

        # Run blind screening regardless — Foldseek searches skip existing results
        BlindScreening(pdb1_name, blind_pred_path)
    else:
        PredictionAll(pdb1_name, search_dir, search_multi_dir, nMSA, model_type)
        logger.info("Finished running predictions using full and shallow random-MSAs")
        logger.info("Running Foldseek to find related crystal structures")
        BlindScreening(pdb1_name, blind_pred_path)


if __name__ == "__main__":
    main()
