#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""2D scatter plotting utilities for fold-switching analysis.

Provides visualization of TM-score distributions and quality metrics
for predicted structures across varied MSA configurations.
"""

import logging

import numpy as np
from matplotlib import (
    pyplot as plt,
)
from numpy import (
    genfromtxt,
)

logger = logging.getLogger(__name__)

# Visualization constants
DEFAULT_FIGURE_SIZE = (10, 8)
DEFAULT_DPI = 300
COLORBAR_LABELSIZE = 15
AXIS_LABELSIZE = 15
TICK_LABELSIZE = 15
PLDT_MIN = 50
PLDT_MAX = 100
SCATTER_SIZE = 35


class Plot2DScatter:
    """Create 2D scatter plots of TM-scores and pLDDT metrics.

    Visualizes the relationship between TM-scores for two conformations
    and color-codes points by predicted local distance difference test (pLDDT)
    score, which indicates model confidence.

    Attributes:
        pdb1_name (str): Name of first reference structure.
        pdb2_name (str): Name of second reference structure.
        nMSA (int): Number of additional MSA samples used.
    """

    def __init__(
        self,
        full_cate: str,
        random_cate: str,
        pdb1: str,
        pdb1_name: str,
        pdb2: str,
        pdb2_name: str,
        nMSA: int,
        nENS: int,
    ) -> None:
        """Initialize and generate scatter plots.

        Args:
            full_cate: Category name for full MSA results.
            random_cate: Category name for random/variable MSA results.
            pdb1: Path to first reference structure (unused, kept for API compatibility).
            pdb1_name: Name of first reference structure (used for filenames).
            pdb2: Path to second reference structure (unused, kept for API compatibility).
            pdb2_name: Name of second reference structure (for axis labels).
            nMSA: Number of additional MSA samples.
            nENS: Number of ensemble samples (unused, kept for API compatibility).

        Raises:
            FileNotFoundError: If required CSV files not found.
            ValueError: If data shapes are incompatible.
        """
        self.pdb1_name = pdb1_name
        self.pdb2_name = pdb2_name
        self.nMSA = nMSA
        self.full_cate = full_cate
        self.random_cate = random_cate

        try:
            self._load_data()
            self._create_plots()
        except Exception as e:
            logger.error(f"Failed to create plots: {e}")
            raise

    def _load_data(self) -> None:
        """Load TM-score and pLDDT data from CSV files.

        Raises:
            FileNotFoundError: If required CSV files cannot be found.
        """
        try:
            # Load TM-scores
            self.tmscore_full = genfromtxt(
                f"TMScore_{self.full_cate}_{self.pdb1_name}.csv", delimiter=" "
            )
            self.tmscore_random = genfromtxt(
                f"TMScore_{self.random_cate}_{self.pdb1_name}.csv", delimiter=" "
            )

            # Load pLDDT scores
            self.plddt_full = genfromtxt(
                f"plddt_{self.full_cate}_{self.pdb1_name}.csv", delimiter=" "
            )
            self.plddt_random = genfromtxt(
                f"plddt_{self.random_cate}_{self.pdb1_name}.csv", delimiter=" "
            )

            # Load fold-switching region TM-scores
            self.tmscore_fs_full = genfromtxt(
                f"TMScore_fs_{self.full_cate}_{self.pdb1_name}.csv", delimiter=" "
            )
            self.tmscore_fs_random = genfromtxt(
                f"TMScore_fs_{self.random_cate}_{self.pdb1_name}.csv", delimiter=" "
            )

            logger.info("Successfully loaded all data files")

        except FileNotFoundError as e:
            logger.error(f"Required CSV file not found: {e}")
            raise

    def _create_plots(self) -> None:
        """Create and save scatter plot visualizations."""
        logger.info("Creating scatter plots...")

        # Reshape pLDDT data for proper color mapping
        self.plddt_random = np.reshape(self.plddt_random, (7, (self.nMSA + 5) * 5))

        # Create whole structure plot
        self._plot_whole_structure()

        # Create fold-switching region plot
        self._plot_foldswitching_region()

    def _plot_whole_structure(self) -> None:
        """Create scatter plot for whole protein structure TM-scores."""
        plt.figure(figsize=DEFAULT_FIGURE_SIZE, dpi=DEFAULT_DPI)

        # Plot variable MSA results
        num_msa_depths = int(self.tmscore_random.shape[0] / 2)
        for ii in range(num_msa_depths):
            plt.scatter(
                self.tmscore_random[ii * 2, :],
                self.tmscore_random[(ii * 2 + 1), :],
                c=self.plddt_random[ii, :],
                cmap="rocket_r",
                vmin=PLDT_MIN,
                vmax=PLDT_MAX,
                s=SCATTER_SIZE,
                marker="o",
                alpha=0.7,
            )

        # Add colorbar
        cbar = plt.colorbar()
        cbar.ax.tick_params(labelsize=COLORBAR_LABELSIZE)
        cbar.set_label("pLDDT Score", fontsize=AXIS_LABELSIZE)

        # Overlay full MSA results
        plt.scatter(
            self.tmscore_full[0, :],
            self.tmscore_full[1, :],
            c=self.plddt_full,
            cmap="plasma",
            vmin=PLDT_MIN,
            vmax=PLDT_MAX,
            s=SCATTER_SIZE,
            marker="D",
            alpha=0.8,
            label="Full MSA",
        )

        # Add diagonal reference line
        plt.plot([0, 1], [0, 1], linestyle="dashed", color="black", linewidth=2)

        # Configure axes
        plt.xlim(0, 1)
        plt.ylim(0, 1)
        plt.xticks(fontsize=TICK_LABELSIZE)
        plt.yticks(fontsize=TICK_LABELSIZE)
        plt.xlabel(f"TM-Score similarity to {self.pdb1_name}", fontsize=AXIS_LABELSIZE)
        plt.ylabel(f"TM-score similarity to {self.pdb2_name}", fontsize=AXIS_LABELSIZE)
        plt.legend(fontsize=TICK_LABELSIZE)
        plt.grid(True, alpha=0.3)

        # Save figure
        output_file = f"TMscore_{self.full_cate}_{self.pdb1_name}.png"
        plt.savefig(output_file, dpi=DEFAULT_DPI, bbox_inches="tight")
        logger.info(f"Saved whole structure plot to {output_file}")
        plt.close()

    def _plot_foldswitching_region(self) -> None:
        """Create scatter plot for fold-switching region TM-scores."""
        plt.figure(figsize=DEFAULT_FIGURE_SIZE, dpi=DEFAULT_DPI)

        # Plot variable MSA results
        num_msa_depths = int(self.tmscore_fs_random.shape[0] / 2)
        for ii in range(num_msa_depths):
            plt.scatter(
                self.tmscore_fs_random[ii * 2, :],
                self.tmscore_fs_random[(ii * 2 + 1), :],
                c=self.plddt_random[ii, :],
                cmap="plasma",
                vmin=PLDT_MIN,
                vmax=PLDT_MAX,
                s=SCATTER_SIZE,
                marker="o",
                alpha=0.7,
            )

        # Add colorbar
        cbar = plt.colorbar()
        cbar.ax.tick_params(labelsize=COLORBAR_LABELSIZE)
        cbar.set_label("pLDDT Score", fontsize=AXIS_LABELSIZE)

        # Overlay full MSA results
        plt.scatter(
            self.tmscore_fs_full[0, :],
            self.tmscore_fs_full[1, :],
            c=self.plddt_full,
            cmap="plasma",
            vmin=PLDT_MIN,
            vmax=PLDT_MAX,
            s=SCATTER_SIZE,
            marker="D",
            alpha=0.8,
            label="Full MSA",
        )

        # Add diagonal reference line
        plt.plot([0, 1], [0, 1], linestyle="dashed", color="black", linewidth=2)

        # Configure axes
        plt.xlim(0, 1)
        plt.ylim(0, 1)
        plt.xticks(fontsize=TICK_LABELSIZE)
        plt.yticks(fontsize=TICK_LABELSIZE)
        plt.xlabel(f"TM-Score similarity to {self.pdb1_name}", fontsize=AXIS_LABELSIZE)
        plt.ylabel(f"TM-score similarity to {self.pdb2_name}", fontsize=AXIS_LABELSIZE)
        plt.legend(fontsize=TICK_LABELSIZE)
        plt.grid(True, alpha=0.3)

        # Save figure
        output_file = f"TMscore_fs-region_{self.full_cate}_{self.pdb1_name}.png"
        plt.savefig(output_file, dpi=DEFAULT_DPI, bbox_inches="tight")
        logger.info(f"Saved fold-switching region plot to {output_file}")
        plt.close()
