"""Module for comparing fold-switching regions in predicted models against reference structures."""

import glob
import logging
import os

import numpy as np
import pandas as pd
from Bio.PDB.PDBParser import PDBParser
from thefuzz import fuzz

logger = logging.getLogger(__name__)


class FSRange:
    """Class for comparing fold-switching regions in predicted models against reference structures."""

    def __init__(self, pdb1: str, pdb2: str, pdb1_name: str, pdb2_name: str, pred_dir: str):
        self.first_res_check(pdb1, pdb2)
        logger.debug(
            "First residue indices — pdb1: %s, pdb2: %s",
            self.pdb1_res_index_1,
            self.pdb2_res_index_1,
        )

        pred_files = glob.glob(str(pred_dir) + "/*_unrelaxed*pdb")
        logger.debug("Prediction directory: %s (%d files found)", pred_dir, len(pred_files))

        self.res_check(pdb1, pdb2, pdb1_name, pdb2_name)
        logger.debug(
            "Crystal FS residues: %s, %s — Predicted FS residues: %s, %s",
            self.crys_fs_res_1_update,
            self.crys_fs_res_2_update,
            self.pred_fs_res_1_update,
            self.pred_fs_res_2_update,
        )

        crys1_fs_res_st = self.crys_fs_res_1_update[0]
        crys1_fs_res_ed = self.crys_fs_res_1_update[1]
        crys2_fs_res_st = self.crys_fs_res_2_update[0]
        crys2_fs_res_ed = self.crys_fs_res_2_update[1]
        pred1_fs_res_st = self.pred_fs_res_1_update[0]
        pred1_fs_res_ed = self.pred_fs_res_1_update[1]
        pred2_fs_res_st = self.pred_fs_res_2_update[0]
        pred2_fs_res_ed = self.pred_fs_res_2_update[1]

        if int(self.pdb1_res_index_1) > 1:
            logger.debug(
                "pdb1 residue numbering does not start at 1 (starts at %s); adjusting FS range",
                self.pdb1_res_index_1,
            )
            self.crys_fs_res_1_update[0] -= int(self.pdb1_res_index_1)
            self.crys_fs_res_1_update[1] -= int(self.pdb1_res_index_1)
            crys1_fs_res_st = self.crys_fs_res_1_update[0]
            crys1_fs_res_ed = self.crys_fs_res_1_update[1]

        if int(self.pdb2_res_index_1) > 1:
            logger.debug(
                "pdb2 residue numbering does not start at 1 (starts at %s); adjusting FS range",
                self.pdb2_res_index_1,
            )
            self.crys_fs_res_2_update[0] -= int(self.pdb2_res_index_1)
            self.crys_fs_res_2_update[1] -= int(self.pdb2_res_index_1)
            crys2_fs_res_st = self.crys_fs_res_2_update[0]
            crys2_fs_res_ed = self.crys_fs_res_2_update[1]

        logger.debug(
            "Adjusted FS ranges — crystal: pdb1=[%s,%s] pdb2=[%s,%s]; predicted: pdb1=[%s,%s] pdb2=[%s,%s]",
            crys1_fs_res_st, crys1_fs_res_ed,
            crys2_fs_res_st, crys2_fs_res_ed,
            pred1_fs_res_st, pred1_fs_res_ed,
            pred2_fs_res_st, pred2_fs_res_ed,
        )

        # Compare secondary structure of predicted models against pdb1
        index = 0
        logger.info("Comparing FS region secondary structure against %s (%d models)", pdb1_name, np.size(pred_files))
        for model in pred_files:
            logger.debug("Processing model: %s", model)
            self.pydssp(pdb1, model, index, pdb1_name)
            dssp_read_tmp = pd.read_csv(
                f"output_{pdb1_name}_{index}.log", sep=" ", header=None
            )
            seq1 = dssp_read_tmp[0].iloc[0]
            seq2 = dssp_read_tmp[0].iloc[1]

            logger.debug(
                "Crystal FS region: %s | Predicted FS region: %s",
                seq1[crys1_fs_res_st:crys1_fs_res_ed],
                seq2[pred2_fs_res_st:pred2_fs_res_ed],
            )

            if (
                fuzz.ratio(
                    seq1[crys1_fs_res_st:crys1_fs_res_ed],
                    seq2[pred2_fs_res_st:pred2_fs_res_ed],
                )
                > 85
            ):
                logger.info("FS region correctly predicted (matched pdb1: %s)", pdb1_name)
                with open(f"fs_compare_output_{pdb1_name}.log", "w", encoding="utf-8") as f:
                    f.write("success")
                break

            elif index == (int(np.size(pred_files)) - 1):
                logger.info(
                    "FS region not matched against %s; retrying against %s", pdb1_name, pdb2_name
                )
                index = 0

                for model in pred_files:
                    self.pydssp(pdb2, model, index, pdb1_name)
                    dssp_read_tmp = pd.read_csv(
                        f"output_{pdb1_name}_{index}.log", sep=" ", header=None
                    )
                    seq1 = dssp_read_tmp[0].iloc[0]
                    seq2 = dssp_read_tmp[0].iloc[1]

                    logger.debug(
                        "Crystal FS region: %s | Predicted FS region: %s",
                        seq1[crys2_fs_res_st:crys2_fs_res_ed],
                        seq2[pred2_fs_res_st:pred2_fs_res_ed],
                    )

                    if (
                        fuzz.ratio(
                            seq1[crys2_fs_res_st:crys2_fs_res_ed],
                            seq2[pred2_fs_res_st:pred2_fs_res_ed],
                        )
                        > 85
                    ):
                        logger.info(
                            "FS region correctly predicted (matched pdb2: %s)", pdb2_name
                        )
                        break
                    elif index == (int(np.size(pred_files)) - 1):
                        logger.warning(
                            "FS region not correctly predicted for %s against either reference",
                            pdb1_name,
                        )
                        with open(
                            f"fs_compare_output_{pdb1_name}.log", "w", encoding="utf-8"
                        ) as f:
                            f.write("fail")
                    else:
                        index += 1

            else:
                index += 1

    def first_res_check(self, pdb1: str, pdb2: str) -> None:
        """Check the first residue index of both PDB files to ensure correct residue numbering."""
        structure_1 = PDBParser().get_structure("pdb1", pdb1)
        structure_2 = PDBParser().get_structure("pdb2", pdb2)

        res_index_1 = [r.id[1] for c in structure_1[0] for r in c.get_residues()]
        res_index_2 = [r.id[1] for c in structure_2[0] for r in c.get_residues()]

        self.pdb1_res_index_1 = int(res_index_1[0])
        self.pdb2_res_index_1 = int(res_index_2[0])
        logger.debug(
            "Parsed first residue indices — pdb1: %d, pdb2: %d",
            self.pdb1_res_index_1,
            self.pdb2_res_index_1,
        )

    def pydssp(self, crys_pdb: str, pred_pdb: str, number: int, pdb_name: str) -> None:
        """Generate and execute the pydssp command to compare secondary structures."""
        command = f"pydssp {crys_pdb} {pred_pdb} -o output_{pdb_name}_{number}.log"
        logger.debug("Executing: %s", command)
        os.system(command)

    def res_check(self, pdb1: str, pdb2: str, pdb1_name: str, pdb2_name: str) -> None:
        """Read fold-switching residue ranges from file and store for the given PDB pair."""
        range_file = os.path.join(os.getcwd(), "range_fs_pairs_all.txt")

        crys_fs_res_1 = crys_fs_res_2 = pred_fs_res_1 = pred_fs_res_2 = ""

        with open(range_file, "r", encoding="utf-8") as file:
            next(file)  # skip header
            for line in file:
                line = line.strip()
                n1, n2, p1, p2, m1, m2 = line.split(",")
                if (n1 == pdb1_name and n2 == pdb2_name) or (n2 == pdb1_name and n1 == pdb2_name):
                    crys_fs_res_1, crys_fs_res_2 = p1, p2
                    pred_fs_res_1, pred_fs_res_2 = m1, m2

        self.crys_fs_res_1_update = [int(i) for i in crys_fs_res_1.split("-")]
        self.crys_fs_res_2_update = [int(i) for i in crys_fs_res_2.split("-")]
        self.pred_fs_res_1_update = [int(i) for i in pred_fs_res_1.split("-")]
        self.pred_fs_res_2_update = [int(i) for i in pred_fs_res_2.split("-")]

        logger.debug(
            "Loaded FS ranges — crystal: %s / %s, predicted: %s / %s",
            self.crys_fs_res_1_update,
            self.crys_fs_res_2_update,
            self.pred_fs_res_1_update,
            self.pred_fs_res_2_update,
        )