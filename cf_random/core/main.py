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

from ..analysis.cal_plddt_ac_fs import (
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

warnings.filterwarnings("ignore")

# Constants
MODEL_TYPES = {
    "ptm": "alphafold2_ptm",
    "monomer": "alphafold2",
    "multimer": "alphafold2_multimer_v3",
}
BLIND = "predictions_all"
SUCCESS = "predictions_all"

ALTERNATIVE_CONFORMATION = "AC"
FOLD_SWITCHING = "FS"
BLIND_MODE = "blind"
FULL_MSA = "full-MSA"
RANDOM_MSA = "random-MSA"

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
        "--num_msa",
        type=str,
        help="number of additional MSA seeds to run (added to default 5)",
    )
    parser.add_argument(
        "--num_ens",
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
        "--model-type",
        type=str,
        choices=["ptm", "monomer", "multimer"],
        help="select model-type of Colabfold: ptm, monomer, or multimer",
    )
    return parser.parse_args()


def resolve_pdb_name(args: argparse.Namespace) -> str:
    """Resolve the working name for blind mode."""
    if args.pdb1 is None and args.pdb2 is None:
        return args.pname

    if args.pdb1 is None and args.pname is None:
        return args.fname.replace("/", "")
    else:
        return args.fname.replace("/", "")


def resolve_num_msa_num_ens(args: argparse.Namespace):
    """Resolve num_msa and num_ens from optional string arguments."""
    num_msa_raw = args.num_msa
    num_ens_raw = args.num_ens

    if num_msa_raw is None and num_ens_raw is None:
        return 0, 0

    if num_msa_raw is not None and num_ens_raw is not None:
        return int(num_msa_raw), int(num_ens_raw)

    if num_msa_raw is None and num_ens_raw is not None:
        return 0, int(num_ens_raw)

    if num_msa_raw is not None and num_ens_raw is None:
        return int(num_msa_raw), 0
    else:
        raise ValueError("Please provide a valid combination of --num_msa and --num_ens")


def resolve_search_dirs(args: argparse.Namespace):
    """Resolve search_dir and search_multi_dir."""
    if args.fname is None and args.fmname is None:
        raise ValueError("--fname (MSA folder) is required for all modes")

    if args.fname is None and args.fmname is not None:
        raise ValueError("--fname (monomer MSA folder) is required alongside --fmname")

    if args.fname is not None and args.fmname is None:
        return args.fname, None
    else:
        return args.fname, " " + args.fmname


def determine_model_type(args: argparse.Namespace, pdb1: Optional[str]) -> str:
    """Determine the model type and set up multimer directory if needed."""
    if args.model_type is None or args.model_type == "ptm":
        return MODEL_TYPES["ptm"]
    elif args.model_type == "monomer":
        return MODEL_TYPES["monomer"]
    elif args.model_type == "multimer" and args.option == BLIND_MODE:
        model_type = MODEL_TYPES["multimer"]
        return model_type
    elif args.model_type == "multimer":
        ter_count = count_chains(pdb1)
        logger.info("%d chain(s) in this multimer file.", ter_count)
        model_type = MODEL_TYPES["multimer"]
        return model_type
    else:
        raise ValueError(
            f"Unrecognized model type: {args.model_type!r}. Choose from: ptm, monomer, multimer"
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

    if args.model_type is None:
        args.model_type = "ptm"

    download_alphafold_params(MODEL_TYPES[args.model_type], Path("."))

    pwd = os.getcwd() + "/"
    pdb1: Optional[str] = None
    pdb2: Optional[str] = None
    pdb1_name: Optional[str] = None
    pdb2_name: Optional[str] = None

    # Resolve working names
    if args.option == BLIND_MODE:
        pdb1_name = resolve_pdb_name(args)
        logger.info("Work name: %s", pdb1_name)

    if args.pdb1 and args.pdb2:
        pdb1 = args.pdb1
        pdb2 = args.pdb2
        pdb1_name = pdb1.replace(".pdb", "")
        pdb2_name = pdb2.replace(".pdb", "")
        logger.info("PDB names: %s, %s", pdb1_name, pdb2_name)

    num_msa, num_ens = resolve_num_msa_num_ens(args)
    search_dir, search_multi_dir = resolve_search_dirs(args)
    model_type = determine_model_type(args, pdb1)

    search_dir = args.fname
    success_dir = f"{SUCCESS}/{pdb1_name}/"

    logger.info(
        "Running CF-random pipeline with updated options: %s",
        {
            "pdb1": pdb1,
            "pdb1_name": pdb1_name,
            "pdb2": pdb2,
            "pdb2_name": pdb2_name,
            "num_msa": num_msa,
            "num_ens": num_ens,
            "model_type": model_type,
            "search_dir": search_dir,
            "search_multi_dir": search_multi_dir,
            "success_dir": success_dir,
        },
    )

    if args.option == ALTERNATIVE_CONFORMATION:
        run_alternative_conformation_workflow(
            pdb1=pdb1,
            pdb1_name=pdb1_name,
            pdb2=pdb2,
            pdb2_name=pdb2_name,
            num_msa=num_msa,
            num_ens=num_ens,
            model_type=model_type,
            search_dir=search_dir,
            search_multi_dir=search_multi_dir,
            success_dir=success_dir,
            pwd=pwd,
        )
    elif args.option == FOLD_SWITCHING:
        run_fold_switching_workflow(
            pdb1=pdb1,
            pdb1_name=pdb1_name,
            pdb2=pdb2,
            pdb2_name=pdb2_name,
            num_msa=num_msa,
            num_ens=num_ens,
            model_type=model_type,
            search_dir=search_dir,
            search_multi_dir=search_multi_dir,
            success_dir=success_dir,
            pwd=pwd,
        )
    elif args.option == BLIND_MODE:
        run_blind_workflow(
            pdb1_name=pdb1_name,
            search_dir=search_dir,
            search_multi_dir=search_multi_dir,
            num_msa=num_msa,
            model_type=model_type,
        )
    else:
        raise ValueError(f"Unrecognized option: {args.option!r}. Choose from: AC, FS, inAC, blind")


def run_alternative_conformation_workflow(
    pdb1: str,
    pdb1_name: str,
    pdb2: str,
    pdb2_name: str,
    num_msa: int,
    num_ens: int,
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
        for _, cur_dir, _ in os.walk(pwd + success_dir + "/"):
            succ_dir_count += len(cur_dir)

    if os.path.exists(success_dir) and 0 < succ_dir_count < 8:
        logger.info("Folder exists but is incomplete — cleaning subfolders")
        shutil.rmtree(success_dir)

    if os.path.exists(success_dir) and succ_dir_count >= 8:
        logger.info("Predictions including full and random-MSA were already completed.")
        calculate_tm_score = TMScoreCalAllVar(
            pdb1=pdb1,
            pdb1_name=pdb1_name,
            pdb2=pdb2,
            pdb2_name=pdb2_name,
            num_msa=num_msa,
            option=ALTERNATIVE_CONFORMATION,
            model_type=model_type,
            search_dir=search_dir,
            search_multi_dir=search_multi_dir,
        )
    else:
        PredictionAll(
            pdb1_name=pdb1_name,
            search_dir=search_dir,
            search_multi_dir=search_multi_dir,
            num_msa=num_msa,
            model_type=model_type,
        )
        calculate_tm_score = TMScoreCalAllVar(
            pdb1=pdb1,
            pdb1_name=pdb1_name,
            pdb2=pdb2,
            pdb2_name=pdb2_name,
            num_msa=num_msa,
            option=ALTERNATIVE_CONFORMATION,
            model_type=model_type,
            search_dir=search_dir,
            search_multi_dir=search_multi_dir,
        )

    shallow_msa_size = np.append([], calculate_tm_score.size_selection)
    logger.info("Specific size of shallow random MSA is similar to full-MSA: %s", shallow_msa_size)
    np.savetxt("selected_MSA-size_" + pdb1_name + ".csv", shallow_msa_size)

    base = os.path.join(pwd, success_dir)
    list_org_samplings = glob.glob(os.path.join(base, "*full_rand*"))
    list_ran_samplings = glob.glob(os.path.join(base, "*max*"))

    logger.info("Searching for pLDDT folders in: %s", base)
    logger.info("Found %d folders for pLDDT calculation", len(list_org_samplings) + len(list_ran_samplings))

    PlddtCal(
        sub_list=list_org_samplings,
        category=FULL_MSA,
        pdb_name=pdb1_name,
        num_msa=num_msa,
        num_ens=num_ens,
        model_type=model_type,
    )
    PlddtCal(
        sub_list=list_ran_samplings,
        category=RANDOM_MSA,
        pdb_name=pdb1_name,
        num_msa=num_msa,
        num_ens=num_ens,
        model_type=model_type,
    )
    Plot2DScatterAC(
        full_category=FULL_MSA,
        random_category=RANDOM_MSA,
        pdb1=pdb1,
        pdb1_name=pdb1_name,
        pdb2=pdb2,
        pdb2_name=pdb2_name,
        num_msa=num_msa,
        num_ens=num_ens,
        model_type=model_type,
    )


def run_fold_switching_workflow(
    pdb1: str,
    pdb1_name: str,
    pdb2: str,
    pdb2_name: str,
    num_msa: int,
    num_ens: int,
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
        for _, cur_dir, _ in os.walk(pwd + success_dir + "/"):
            succ_dir_count += len(cur_dir)

    if os.path.exists(success_dir) and 0 < succ_dir_count < 8:
        logger.info("Folder exists but is incomplete — cleaning subfolders")
        shutil.rmtree(success_dir)

    shallow_msa_size = np.array([])

    logging.info("Success directory and count: %s, %d", success_dir, succ_dir_count)

    if os.path.exists(success_dir) and succ_dir_count >= 8:
        logger.info("Predictions including full and random-MSA were already completed.")
        calculate_tm_score = TMScoreCalAllVarFS(
            pdb1=pdb1,
            pdb1_name=pdb1_name,
            pdb2=pdb2,
            pdb2_name=pdb2_name,
            num_msa=num_msa,
            option=FOLD_SWITCHING,
            model_type=model_type,
            search_dir=search_dir,
            search_multi_dir=search_multi_dir,
        )
        shallow_msa_size = np.append(shallow_msa_size, calculate_tm_score.size_selection)
    else:
        PredictionAll(
            pdb1_name=pdb1_name,
            search_dir=search_dir,
            search_multi_dir=search_multi_dir,
            num_msa=num_msa,
            model_type=model_type,
        )
        calculate_tm_score = TMScoreCalAllVarFS(
            pdb1=pdb1,
            pdb1_name=pdb1_name,
            pdb2=pdb2,
            pdb2_name=pdb2_name,
            num_msa=num_msa,
            option=FOLD_SWITCHING,
            model_type=model_type,
            search_dir=search_dir,
            search_multi_dir=search_multi_dir,
        )
        shallow_msa_size = np.append(shallow_msa_size, calculate_tm_score.size_selection)

        logger.info(
            "Specific size of shallow random MSA is similar to full-MSA: %s", shallow_msa_size
        )
        np.savetxt("selected_MSA-size_" + pdb1_name + ".csv", shallow_msa_size)

    base = os.path.join(pwd, success_dir)
    list_org_samplings = glob.glob(os.path.join(base, "*full_rand*"))
    list_ran_samplings = glob.glob(os.path.join(base, "*max*"))

    logger.info("Searching for pLDDT folders in: %s", base)
    logger.info("Found %d folders for pLDDT calculation", len(list_org_samplings) + len(list_ran_samplings))

    PlddtCal(
        sub_list=list_org_samplings,
        category=FULL_MSA,
        pdb_name=pdb1_name,
        num_msa=num_msa,
        num_ens=num_ens,
        model_type=model_type,
    )
    PlddtCal(
        sub_list=list_ran_samplings,
        category=RANDOM_MSA,
        pdb_name=pdb1_name,
        num_msa=num_msa,
        num_ens=num_ens,
        model_type=model_type,
    )

    if model_type == MODEL_TYPES["multimer"]:
        Plot2DScatterAC(
            full_category=FULL_MSA,
            random_category=RANDOM_MSA,
            pdb1=pdb1,
            pdb1_name=pdb1_name,
            pdb2=pdb2,
            pdb2_name=pdb2_name,
            num_msa=num_msa,
            num_ens=num_ens,
            model_type=model_type,
        )
    else:
        Plot2DScatter(
            full_category=FULL_MSA,
            random_category=RANDOM_MSA,
            pdb1=pdb1,
            pdb1_name=pdb1_name,
            pdb2=pdb2,
            pdb2_name=pdb2_name,
            num_msa=num_msa,
            num_ens=num_ens,
        )


def run_blind_workflow(
    pdb1_name: str,
    search_dir: str,
    search_multi_dir: Union[int, str],
    num_msa: int,
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
        for _, cur_dir, _ in os.walk(blind_pdb_dir + "/"):
            blind_dir_count += len(cur_dir)

    if os.path.exists(blind_pdb_dir) and 0 < blind_dir_count < 8:
        logger.info("Folder exists but is incomplete — cleaning subfolders")
        shutil.rmtree(blind_pdb_dir)

    if os.path.exists(blind_pdb_dir) and blind_dir_count >= 8:
        logger.info("Predictions including full and random-MSA were already completed.")

        # Count total files to determine whether Foldseek has already been run
        fseek_file_count = 0
        for _, cur_dir, files in os.walk(blind_pdb_dir + "/"):
            fseek_file_count += len(files)
        logger.info("Number of files in blind prediction path: %d", fseek_file_count)

        if fseek_file_count >= FOLDSEEK_DONE_FILE_COUNT:
            logger.info("Foldseek search was already done")

        # Run blind screening regardless — Foldseek searches skip existing results
        BlindScreening(pdb1_name=pdb1_name, blind_path=blind_pred_path)
    else:
        PredictionAll(
            pdb1_name=pdb1_name,
            search_dir=search_dir,
            search_multi_dir=search_multi_dir,
            num_msa=num_msa,
            model_type=model_type,
        )
        logger.info("Finished running predictions using full and shallow random-MSAs")
        logger.info("Running Foldseek to find related crystal structures")
        BlindScreening(pdb1_name=pdb1_name, blind_path=blind_pred_path)


if __name__ == "__main__":
    main()
