#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""pLDDT score calculation for AlphaFold predictions."""

import json
import logging
import re
from pathlib import (
    Path,
)
from typing import (
    Dict,
    List,
)

import numpy as np

logger = logging.getLogger(__name__)

# Regex pattern for parsing AlphaFold JSON filenames
RANK_PATTERN = re.compile(r".*?_scores_rank_(?P<rank>\d+)_alphafold2.*")


def read_plddt(jsonfile: str) -> np.ndarray:
    """Reads pLDDT scores from an AlphaFold prediction JSON file.

    Args:
        jsonfile: Path to the JSON file containing prediction scores.

    Returns:
        Array of pLDDT scores as float64.
    """
    with open(jsonfile) as json_file:
        data = json.load(json_file)

    plddt_scores = np.array(data["plddt"], dtype=np.float64)
    return plddt_scores


def calculate_average_plddt(score: np.ndarray) -> float:
    """Calculates the average pLDDT score from an array of scores.

    Args:
        score: Array of pLDDT scores.

    Returns:
        Average pLDDT score rounded to 2 decimal places.
    """
    avg_plddt = round(np.average(score), 2)
    return avg_plddt


class PlddtCal:
    """Calculates average pLDDT scores for protein models across different categories.

    This class processes AlphaFold prediction JSON files to extract pLDDT scores
    and computes average scores for different MSA categories and model types.
    """

    def __init__(
        self,
        sub_list: List[str],
        category: str,
        pdb_name: str,
        nMSA: int,
        nENS: int,
        model_type: str,
    ) -> None:
        """Initializes pLDDT calculation for given subdirectories and parameters.

        Args:
            sub_list: List of subdirectory paths to process.
            category: MSA category ('full-MSA', 'additional-MSA', 'random-MSA').
            pdb_name: Name of the PDB structure.
            nMSA: Number of MSA sequences.
            nENS: Number of ensemble models.
            model_type: Type of AlphaFold model.
        """
        if not sub_list:
            raise ValueError("No subdirectories provided for pLDDT calculation")

        logger.info("Processing pLDDT scores for %d subdirectories", len(sub_list))
        logger.debug("Subdirectories: %s", sub_list)

        values_all, out_dict_all, cnt = self._process_subdirs(sub_list)

        if category == "full-MSA":
            cnt = int(cnt / 5)

        logger.debug("Processed %d files", cnt)

        # Reshape based on category and model_type
        if category == "full-MSA":
            values_all_resh = values_all.reshape(nMSA + 5, 5)
        elif category == "additional-MSA" and model_type == "alphafold2_multimer_v3":
            values_all_resh = values_all.reshape((nENS + 20), 5)
        elif category == "additional-MSA":
            values_all_resh = values_all.reshape((nENS + 20), 5)
        elif category == "random-MSA" and model_type != "alphafold2_multimer_v3":
            values_all_resh = values_all.reshape((nMSA + 5) * 7, 5)
        elif category == "random-MSA":
            values_all_resh = values_all.reshape((nMSA + 5) * 7, 5)
        else:
            raise ValueError(f"Unknown category/model_type combination: {category}/{model_type}")

        logger.info("Calculated pLDDT scores")
        output_file = f"plddt_{category}_{pdb_name}.csv"
        np.savetxt(output_file, values_all_resh, fmt="%2.3f")
        logger.info("Saved pLDDT results to %s", output_file)

    def _process_subdirs(self, sub_list: List[str]) -> tuple[np.ndarray, Dict[str, float], int]:
        """Processes subdirectories to extract pLDDT scores.

        Args:
            sub_list: List of subdirectory paths.

        Returns:
            Tuple of (values_all, out_dict_all, cnt) where values_all is numpy array
            of scores, out_dict_all is dict of key-value pairs, cnt is count.
        """
        out_dict_all: Dict[str, float] = {}
        values_all = np.array([])
        cnt = 0

        for subdir in sub_list:
            subdir_path = Path(subdir)
            if not subdir_path.is_dir():
                logger.warning("Skipping non-directory: %s", subdir)
                continue

            subdir_name = subdir_path.name
            jsonfiles = list(subdir_path.glob("*_scores*json"))

            for jsonfile in jsonfiles:
                plddt_score = read_plddt(str(jsonfile))
                values = calculate_average_plddt(plddt_score)
                values_all = np.append(values_all, values)

                jsonfilename = jsonfile.stem
                match = RANK_PATTERN.match(jsonfilename)
                rank = match.group("rank") if match else "000"

                key_pair = f"{subdir_name}:{rank}"
                if key_pair not in out_dict_all:
                    out_dict_all[key_pair] = values

                cnt += 1

        return values_all, out_dict_all, cnt
