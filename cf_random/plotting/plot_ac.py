#!/bin/env python3
# -*- coding: utf-8 -*-
"""Plotting utilities for alternative TM-score comparisons.

Provides plotting helpers used by the command-line tools.
"""

import logging

import numpy as np
from matplotlib import pyplot as plt
from numpy import genfromtxt

logger = logging.getLogger(__name__)


class Plot2DScatterAC:
    """Create 2D scatter plots of TM-scores for alternative conformations."""

    def __init__(
        self,
        full_cate: str,
        random_cate: str,
        pdb1: str,
        pdb1_name: str,
        pdb2: str,
        pdb2_name: str,
        num_msa: int,
        num_ens: int,
        model_type: str,
    ):
        # Load TM-scores both full- and random-MSA
        tmscores_full = genfromtxt(
            f"TMScore_{full_cate}_{pdb1_name}.csv", delimiter=" "
        )
        tmscores_random = genfromtxt(
            f"TMScore_{random_cate}_{pdb1_name}.csv", delimiter=" "
        )

        # Load pLDDT scores both full- and random-MSA
        plddt_full = genfromtxt(f"plddt_{full_cate}_{pdb1_name}.csv", delimiter=" ")
        plddt_random = genfromtxt(f"plddt_{random_cate}_{pdb1_name}.csv", delimiter=" ")

        logger.debug(
            "TM-score array shape: rows=%d, cols=%d, ndim=%d",
            tmscores_random.shape[0],
            tmscores_random.shape[-1],
            tmscores_random.ndim,
        )

        plddt_random = np.reshape(plddt_random, (7, (num_msa + 5) * 5))
        tmscore_full_resh = np.reshape(tmscores_full, (((num_msa + 5) * 2), 5))

        self._plot(
            tmscores_random,
            tmscore_full_resh,
            plddt_random,
            plddt_full,
            num_msa,
            pdb1_name,
            pdb2_name,
            full_cate,
        )

    def _plot(
        self,
        tmscores_random: np.ndarray,
        tmscore_full_resh: np.ndarray,
        plddt_random: np.ndarray,
        plddt_full: np.ndarray,
        num_msa: int,
        pdb1_name: str,
        pdb2_name: str,
        full_cate: str,
    ) -> None:
        """Render and save the 2D TM-score scatter plot."""
        plt.figure(0)
        for i in range(0, int(tmscores_random.shape[0] / 2)):
            plt.scatter(
                tmscores_random[i * 2, :],
                tmscores_random[(i * 2 + 1), :],
                c=plddt_random[i, :],
                cmap="rocket_r",
                vmin=50,
                vmax=100,
                s=35,
                marker="o",
            )

        clb = plt.colorbar()
        clb.ax.tick_params(labelsize=15)

        plt.scatter(
            tmscore_full_resh[0 : (num_msa + 5), :],
            tmscore_full_resh[(num_msa + 5) : (num_msa + 5) * 2, :],
            c=plddt_full,
            cmap="rocket_r",
            vmin=50,
            vmax=100,
            s=35,
            marker="o",
        )

        plt.ylim(0, 1)
        plt.xlim(0, 1)
        plt.plot([0, 1], [0, 1], linestyle="dashed", color="black")
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        plt.xlabel(f"TM-Score similar to fold1({pdb1_name})", fontsize=15)
        plt.ylabel(f"TM-score similar to fold2({pdb2_name})", fontsize=15)

        output_file = f"TMscore_{full_cate}_{pdb1_name}.png"
        plt.savefig(output_file, transparent=True)
        logger.info("Saved AC scatter plot to %s", output_file)