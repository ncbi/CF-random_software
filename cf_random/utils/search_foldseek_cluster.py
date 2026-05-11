#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Blind structure screening using Foldseek and clustering analysis.

This module performs unsupervised clustering of predicted structures using
Foldseek bit-score similarity matrices, PCA dimensionality reduction, HDBSCAN
clustering, and k-medoids representative selection.

Requires:
    - Foldseek: For structure similarity computation
    - MDAnalysis: For structure analysis and DSSP filtering
    - scikit-learn: For clustering algorithms
    - PyMOL (optional): For .pse session file output
"""

import csv
import logging
import re
import shutil
import subprocess
from pathlib import (
    Path,
)
from typing import (
    Dict,
    List,
    Tuple,
)

import matplotlib.pyplot as plt
import MDAnalysis as mda
import numpy as np
from MDAnalysis.analysis.dssp import (
    DSSP,
)
from scipy import (
    stats,
)
from scipy.spatial import (
    distance,
)
from sklearn.cluster import (
    HDBSCAN,
)
from sklearn.decomposition import (
    PCA,
)
from sklearn.metrics import (
    silhouette_score,
)
from sklearn.preprocessing import (
    minmax_scale,
)

logger = logging.getLogger(__name__)

try:
    import pymol
except ImportError:
    pymol = None
    logger.warning("PyMOL not available; .pse session file will not be saved")

# Configuration constants
HDBSCAN_K_RANGE = range(2, 51)
HDBSCAN_MIN_SAMPLES = 1
KMEDOIDS_DEFAULT_K = 3
KMEDOIDS_MIN_CLUSTER_SIZE = 4
KMEDOIDS_MAX_ITER = 100
RANDOM_SEED = 42
DISTANCE_METRIC = "euclidean"
ZSCORE_OUTLIER_THRESHOLD = 3.0
FOLDSEEK_SENSITIVITY = "9.5"
FOLDSEEK_FORMAT_OUTPUT = "query,target,alntmscore,qaln,taln,alnlen,evalue,bits"
FOLDSEEK_BIT_SCORE_BUG_VALUE = -2147483648
PCA_N_COMPONENTS = 4


class BlindScreening:
    """Perform blind structural screening and clustering analysis.

    1. Stage PDB files into a flat directory and build a Foldseek DB.
    2. Run per-file exhaustive Foldseek easy-search against the DB.
    3. Parse bit scores into a pairwise correlation matrix.
    4. Remove unfolded outliers via DSSP z-score filtering.
    5. Normalise the matrix, run PCA, HDBSCAN, and k-medoids.
    6. Write CSV summaries, a cluster PNG, and (optionally) a PyMOL .pse.
    """

    def __init__(self, pdb1_name: str, blind_path: str) -> None:
        """Initialise and run the full blind screening pipeline."""
        self.pdb1_name = pdb1_name
        self.blind_path = Path(blind_path)

        if not self.blind_path.exists():
            raise FileNotFoundError(f"Blind screening path not found: {blind_path}")

        logger.info("Starting blind screening for %s", pdb1_name)

        self._stage_pdb_files()
        self._build_foldseek_database()
        self._run_foldseek_searches()
        self._perform_clustering_analysis()

        logger.info("Blind screening completed successfully")

    @staticmethod
    def cluster_structures(x: np.ndarray) -> np.ndarray:
        """Find optimal HDBSCAN clustering labels for reduced structure features.

        Iterates over a range of min_cluster_size values and selects the value
        that maximises the silhouette score.

        Args:
            x: PCA-reduced feature matrix (n_samples, n_components).

        Returns:
            np.ndarray: Cluster labels for each sample (-1 = noise).
        """
        sil_scores: List[float] = []
        for k in HDBSCAN_K_RANGE:
            clustering = HDBSCAN(min_cluster_size=k, min_samples=HDBSCAN_MIN_SAMPLES)
            clustering.fit(x)
            n_unique = len(set(clustering.labels_))
            if 1 < n_unique < len(x):
                score = silhouette_score(x, clustering.labels_, metric=DISTANCE_METRIC)
                sil_scores.append(score)
            else:
                sil_scores.append(-1.0)

        opt_k = HDBSCAN_K_RANGE[int(np.argmax(sil_scores))]
        logger.info("Optimal HDBSCAN min_cluster_size: %s", opt_k)

        final = HDBSCAN(min_cluster_size=opt_k, min_samples=HDBSCAN_MIN_SAMPLES)
        final.fit(x)
        return final.labels_

    @staticmethod
    def k_medoids(
        x: np.ndarray,
        cluster_label: int,
        labels: np.ndarray,
        k: int = KMEDOIDS_DEFAULT_K,
        max_iter: int = KMEDOIDS_MAX_ITER,
    ) -> Tuple[np.ndarray, float]:
        """PAM-style k-medoids to find representative structures in a cluster.

        Args:
            x:             Full PCA-reduced feature matrix (all samples).
            cluster_label: The HDBSCAN label whose members are to be processed.
            labels:        Full label array from HDBSCAN (length = n_samples).
            k:             Number of medoids to return.
            max_iter:      Maximum PAM swap iterations.

        Returns:
            Tuple of (medoid_indices, total_cost) where medoid_indices indexes
            into the *full* x array (not only the cluster subset).
        """
        np.random.seed(RANDOM_SEED)

        temp = x.copy()
        mask = np.zeros(x.shape, dtype=bool)
        mask[np.argwhere(labels == cluster_label)] = True

        # Count members of this cluster
        unique_vals, counts = np.unique(mask[:, 0], return_counts=True)
        true_idx = [i for i, v in enumerate(unique_vals) if v]
        if not true_idx:
            return np.array([], dtype=int), float("nan")
        cluster_count = counts[true_idx[0]]

        if cluster_count < KMEDOIDS_MIN_CLUSTER_SIZE:
            logger.debug(
                "Cluster %s has only %s members; returning all indices",
                cluster_label,
                cluster_count,
            )
            return np.ravel(np.argwhere(mask[:, 0])), float("nan")

        # Mask out non-cluster members so they never win a distance comparison
        temp[~mask] = 9999.0

        n = temp.shape[0]
        medoids = np.random.choice(n, k, replace=False)
        c_dis = distance.cdist(temp, temp[medoids], metric=DISTANCE_METRIC)
        tot_cost = float(np.sum(np.min(c_dis, axis=1)))

        for _ in range(max_iter):
            improved = False
            for m_idx in range(k):
                for candidate in range(n):
                    if candidate in medoids:
                        continue
                    new_medoids = medoids.copy()
                    new_medoids[m_idx] = candidate
                    new_dis = distance.cdist(temp, temp[new_medoids], metric=DISTANCE_METRIC)
                    new_cost = float(np.sum(np.min(new_dis, axis=1)))
                    if new_cost < tot_cost:
                        medoids = new_medoids
                        tot_cost = new_cost
                        improved = True
                        break
                if improved:
                    break
            if not improved:
                break

        return medoids, tot_cost

    def _stage_pdb_files(self) -> None:
        """Copy all PDB files into a flat staging directory for Foldseek.

        Replicates the original's flattening logic: replace '/' with '-' in the
        full path, then strip the first 17 characters to derive the DB filename.
        The mapping from source Path -> flat label is stored in
        ``self.pdb_label_map`` so that ``_build_correlation_matrix`` can match
        Foldseek target fields back to matrix rows.
        """
        self.db_directory = self.blind_path / "pdbs_for_db"
        self.db_directory.mkdir(exist_ok=True)

        raw_pdb_files = sorted(self.blind_path.rglob("*.pdb"))
        raw_pdb_files = [f for f in raw_pdb_files if self.db_directory not in f.parents]

        if not raw_pdb_files:
            raise FileNotFoundError(f"No PDB files found under {self.blind_path}")

        # pdb_label_map: source Path -> flat label used as Foldseek DB key
        self.pdb_label_map: Dict[Path, str] = {}

        self.pdb_files: List[Path] = []
        logger.info("Staging %d PDB files for Foldseek database", len(raw_pdb_files))

        for src in raw_pdb_files:
            dest_name = str(src).replace("/", "-")[17:]
            dest = self.db_directory / dest_name
            if not dest.exists():
                shutil.copyfile(src, dest)
            self.pdb_files.append(src)
            # Label is the dest_name without the .pdb extension
            self.pdb_label_map[src] = dest_name.replace(".pdb", "")

        logger.info("Staged %d files to %s", len(self.pdb_files), self.db_directory)

    def _build_foldseek_database(self) -> None:
        """Create a Foldseek database from the staged PDB directory."""
        self.foldseek_db = self.db_directory / "DB"

        if self.foldseek_db.exists():
            logger.info(
                "Existing Foldseek database found at %s; skipping creation", self.foldseek_db
            )
            return

        cmd = ["foldseek", "createdb", str(self.db_directory), str(self.foldseek_db)]
        logger.info("Building Foldseek database...")
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"foldseek createdb failed:\n{result.stderr.strip()}")
        logger.info("Foldseek database created at %s", self.foldseek_db)

    def _run_foldseek_searches(self) -> None:
        """Run per-file exhaustive Foldseek easy-search against the database."""
        tmp_dir = self.blind_path / "tmp"
        tmp_dir.mkdir(exist_ok=True)

        for pdb_file in self.pdb_files:
            result_file = pdb_file.with_suffix("").parent / (pdb_file.stem + "-self.foldseek")
            if result_file.exists():
                logger.debug("Foldseek result exists; skipping %s", pdb_file.name)
                continue

            cmd = [
                "foldseek",
                "easy-search",
                str(pdb_file),
                str(self.foldseek_db),
                str(result_file),
                str(tmp_dir),
                "--format-mode",
                "0",
                "--format-output",
                FOLDSEEK_FORMAT_OUTPUT,
                "--exhaustive-search",
                "1",
                "-s",
                FOLDSEEK_SENSITIVITY,
            ]
            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode != 0:
                logger.warning(
                    "foldseek easy-search failed for %s:\n%s",
                    pdb_file.name,
                    result.stderr.strip(),
                )
            else:
                logger.info("Foldseek search completed for %s", pdb_file.name)

    def _perform_clustering_analysis(self) -> None:
        """Build similarity matrix, filter outliers, cluster, and save outputs."""
        foldseek_files = sorted(self.blind_path.rglob("*-self.foldseek"))
        if not foldseek_files:
            raise FileNotFoundError("No .foldseek result files found")

        # _filter_unfolded returns files sorted; labels are derived after sort
        foldseek_files, pdb_labels = self._filter_unfolded(foldseek_files)

        corr_mtx = self._build_correlation_matrix(foldseek_files, pdb_labels)

        norm = minmax_scale(corr_mtx, axis=1)
        norm = (norm + norm.T) / 2.0

        pca_coords = PCA(n_components=PCA_N_COMPONENTS).fit_transform(norm)
        labels = self.cluster_structures(pca_coords)

        files_of_interest: List[Tuple[Path, int]] = []
        pca_of_interest: List[np.ndarray] = []

        for cluster_label in np.unique(labels):
            medoid_indices, _ = self.k_medoids(pca_coords, cluster_label, labels)
            for idx in medoid_indices:
                files_of_interest.append((foldseek_files[idx], int(cluster_label)))
                pca_of_interest.append(pca_coords[idx])

        self._save_cluster_plot(pca_coords, labels)
        self._save_structures_of_interest(files_of_interest, pca_of_interest)
        self._save_all_structures(foldseek_files, labels, pca_coords)
        self._save_pymol_session(foldseek_files, files_of_interest)

    def _filter_unfolded(self, foldseek_files: List[Path]) -> Tuple[List[Path], List[str]]:
        """Remove unfolded predictions using DSSP loop-content z-scores.

        A structure is flagged as an outlier if its loop ('-') residue count
        has a z-score > ZSCORE_OUTLIER_THRESHOLD across all structures.

        Returns:
            Filtered (sorted) foldseek file list and matching flat PDB labels.
        """
        files_count: List[np.ndarray] = []

        for ff in foldseek_files:
            pdb_path = Path(str(ff).replace("-self.foldseek", ".pdb"))
            u = mda.Universe(str(pdb_path))
            s = DSSP(u).run().results.dssp[0]
            dssp_types, counts = np.unique(s, return_counts=True)

            # Ensure all three categories are present
            for missing, pos in [("-", 0), ("E", 1), ("H", 2)]:
                if missing not in dssp_types:
                    dssp_types = np.insert(dssp_types, pos, missing)
                    counts = np.insert(counts, pos, 0)

            files_count.append(counts)

        files_count_arr = np.array(files_count)
        z_scores = stats.zscore(files_count_arr[:, 0])  # column 0 = '-' (loops)

        filtered: List[Path] = []
        for i, ff in enumerate(foldseek_files):
            if z_scores[i] > ZSCORE_OUTLIER_THRESHOLD:
                logger.info(
                    "Removed unfolded structure from analysis: %s",
                    str(ff).replace("-self.foldseek", ".pdb"),
                )
            else:
                filtered.append(ff)

        filtered = sorted(filtered)

        # Derive flat labels using the same slash->dash + [17:] rule as staging,
        # applied to the .pdb path so they match what Foldseek indexed.
        pdb_labels = [
            str(Path(str(ff).replace("-self.foldseek", ".pdb")))
            .replace("/", "-")[17:]
            .replace(".pdb", "")
            for ff in filtered
        ]

        return filtered, pdb_labels

    def _build_correlation_matrix(
        self, foldseek_files: List[Path], pdb_labels: List[str]
    ) -> np.ndarray:
        """Parse bit scores from .foldseek files into a pairwise matrix."""
        corr_mtx: List[List[float]] = []

        for ff in foldseek_files:
            row_dict: Dict[str, float] = {label: 0.0 for label in pdb_labels}

            with ff.open("r") as fh:
                for line in fh:
                    parts = line.rstrip().split("\t")
                    if len(parts) < 8:
                        continue
                    target = parts[1]
                    try:
                        bit_score = int(parts[-1])
                    except ValueError:
                        continue
                    if bit_score == FOLDSEEK_BIT_SCORE_BUG_VALUE:
                        bit_score = 0
                    if target in row_dict:
                        row_dict[target] = float(bit_score)

            corr_mtx.append([row_dict[label] for label in pdb_labels])

        return np.array(corr_mtx)

    def _save_cluster_plot(self, pca_coords: np.ndarray, labels: np.ndarray) -> None:
        """Save a 2D PCA scatter plot coloured by HDBSCAN cluster."""
        plot_path = self.blind_path / f"{self.pdb1_name}-cluster.png"
        plt.figure(figsize=(8, 6))
        plt.scatter(pca_coords[:, 0], pca_coords[:, 1], c=labels, cmap="viridis", s=45)
        plt.xlabel("PC 1")
        plt.ylabel("PC 2")
        plt.title(f"Blind screening cluster map — {self.pdb1_name}")
        plt.tight_layout()
        plt.savefig(plot_path)
        plt.clf()
        logger.info("Saved cluster plot to %s", plot_path)

    def _save_structures_of_interest(
        self,
        files_of_interest: List[Tuple[Path, int]],
        pca_of_interest: List[np.ndarray],
    ) -> None:
        """Write the structures-of-interest CSV (group, file, pca_1, pca_2)."""
        out_path = self.blind_path / f"{self.pdb1_name}-structures_of_interest.csv"
        with out_path.open("w", newline="") as fh:
            writer = csv.writer(fh)
            writer.writerow(["group", "file", "pca_1", "pca_2"])
            for (ff, cluster_label), pca_pt in zip(files_of_interest, pca_of_interest):
                writer.writerow([cluster_label, str(ff), pca_pt[0], pca_pt[1]])
        logger.info("Saved structures of interest to %s", out_path)

    def _save_all_structures(
        self,
        foldseek_files: List[Path],
        labels: np.ndarray,
        pca_coords: np.ndarray,
    ) -> None:
        """Write cluster assignments for every structure to structures_all.csv."""
        out_path = self.blind_path / "structures_all.csv"
        with out_path.open("w", newline="") as fh:
            writer = csv.writer(fh)
            writer.writerow(["group", "file", "pca_1", "pca_2"])
            for i, ff in enumerate(foldseek_files):
                writer.writerow([int(labels[i]), str(ff), pca_coords[i, 0], pca_coords[i, 1]])
        logger.info("Saved all-structures summary to %s", out_path)

    def _save_pymol_session(
        self,
        foldseek_files: List[Path],
        files_of_interest: List[Tuple[Path, int]],
    ) -> None:
        """Align structures of interest and save a PyMOL .pse session file.

        Uses foldseek_files[0] as the alignment reference.

        Skipped silently if PyMOL is not installed.
        """
        if pymol is None:
            logger.warning("PyMOL unavailable; skipping .pse export")
            return

        pse_path = self.blind_path / f"{self.pdb1_name}-structures_of_interest.pse"
        viridis = plt.get_cmap("viridis", len(files_of_interest))
        largest_label = max(files_of_interest, key=lambda x: x[1])[1]

        # Load the first of all filtered files as the alignment reference
        dominant_pdb = str(foldseek_files[0]).replace("-self.foldseek", ".pdb")
        pymol.cmd.load(dominant_pdb, "Dominant")

        for idx, (ff, cluster_label) in enumerate(files_of_interest):
            pdb_path = str(ff).replace("-self.foldseek", ".pdb")

            tokens = re.findall(r"(full)|(max\w+)|(rank_\d+)", str(ff))
            short = "_".join(t for group in tokens for t in group if t)
            obj_name = f"{idx}_{short}" if short else f"struct_{idx}"

            if largest_label == -1:
                colour_val = 0.0
            else:
                colour_val = (cluster_label + 1) / (largest_label + 1)
            rgb = viridis(colour_val)[:3]

            pymol.cmd.load(pdb_path, obj_name)
            pymol.cmd.align(obj_name, "Dominant")
            colour_name = f"col_{cluster_label}"
            pymol.cmd.set_color(colour_name, list(rgb))
            pymol.cmd.color(colour_name, obj_name)

        pymol.cmd.save(str(pse_path), "pse")
        pymol.cmd.delete("all")
        pymol.cmd.reinitialize()
        logger.info("Saved PyMOL session to %s", pse_path)
