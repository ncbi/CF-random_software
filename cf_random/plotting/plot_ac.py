#!/bin/env python3
# -*- coding: utf-8 -*-
"""Plotting utilities for alternative TM-score comparisons.

Provides plotting helpers used by the command-line tools.
"""

import numpy as np
from matplotlib import (
    pyplot as plt,
)
from numpy import (
    genfromtxt,
)


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
        tmscores_full = genfromtxt("TMScore_" + full_cate + "_" + pdb1_name + ".csv", delimiter=" ")
        tmscores_random = genfromtxt(
            "TMScore_" + random_cate + "_" + pdb1_name + ".csv", delimiter=" "
        )

        # Load pLDDT scores both full- and random-MSA
        plddt_full = genfromtxt("plddt_" + full_cate + "_" + pdb1_name + ".csv", delimiter=" ")
        plddt_random = genfromtxt("plddt_" + random_cate + "_" + pdb1_name + ".csv", delimiter=" ")

        # Plotting the TM-score values as 2D scatter plot
        print("Size of column: ", tmscores_random.shape[-1])
        print("Size of row: ", tmscores_random.shape[0])
        print("Dimension: ", tmscores_random.ndim)

        print(f"Random MSA TM-scores: {tmscores_random}")
        print(f"Full MSA TM-scores: {tmscores_full}")

        print("Checking plddt")
        print(f"Full MSA pLDDT scores: {plddt_full}")
        print(f"Random MSA pLDDT scores: {plddt_random}")

        plddt_random = np.reshape(plddt_random, (7, (num_msa + 5) * 5))
        print(f"Reshaped Random MSA pLDDT scores: {plddt_random}")

        if model_type != "alphafold2_multimer_v3":
            tmscore_full_resh = np.reshape(tmscores_full, ((((num_msa + 5) * 2), 5)))
        else:
            tmscore_full_resh = np.reshape(tmscores_full, (((num_msa + 5) * 2), 5))

        if model_type != "alphafold2_multimer_v3":

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

            x = [0, 1]
            y = [0, 1]

            plt.ylim(0, 1)
            plt.xlim(0, 1)

            plt.plot(x, y, linestyle="dashed", color="black")

            plt.xticks(fontsize=15)
            plt.yticks(fontsize=15)

            plt.xlabel("TM-Score similar to fold1(" + pdb1_name + ")", fontsize=15)
            plt.ylabel("TM-score similar to fold2(" + pdb2_name + ")", fontsize=15)
            plt.savefig("TMscore_" + full_cate + "_" + pdb1_name + ".png", transparent=True)

        else:
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

            x = [0, 1]
            y = [0, 1]

            plt.ylim(0, 1)
            plt.xlim(0, 1)

            plt.plot(x, y, linestyle="dashed", color="black")

            plt.xticks(fontsize=15)
            plt.yticks(fontsize=15)

            plt.xlabel("TM-Score similar to fold1(" + pdb1_name + ")", fontsize=15)
            plt.ylabel("TM-score similar to fold2(" + pdb2_name + ")", fontsize=15)
            plt.savefig("TMscore_" + full_cate + "_" + pdb1_name + ".png", transparent=True)
