#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Blind structure screening using Foldseek and clustering analysis.

This module performs unsupervised clustering of predicted structures using
Foldseek similarity scores, PCA dimensionality reduction, and HDBSCAN
clustering followed by k-medoids representative selection.

Requires:
    - Foldseek: For structure similarity computation
    - MDAnalysis: For structure analysis
    - scikit-learn: For clustering algorithms
"""

import csv
import logging
import subprocess
from pathlib import (
    Path,
)
from typing import (
    List,
    Tuple,
)

import matplotlib.pyplot as plt
import MDAnalysis as mda
import numpy as np
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

logger = logging.getLogger(__name__)

try:
    import pymol
except ImportError:
    pymol = None
    logger.warning("PyMOL not available; visualization features will be disabled")

# Configuration constants
HDBSCAN_K_RANGE = range(2, 51)
HDBSCAN_MIN_SAMPLES = 1
KMEDOIDS_DEFAULT_K = 3
KMEDOIDS_MIN_CLUSTER_SIZE = 4
RANDOM_SEED = 42
DISTANCE_METRIC = "euclidean"


class BlindScreening:
    """Perform blind structural screening and clustering analysis.

    Automatically identifies representative structures from a set of
    predictions using unsupervised clustering.
    """

    @staticmethod
    def cluster_structures(X: np.ndarray) -> np.ndarray:
        """Find optimal HDBSCAN clustering labels for reduced structure features.

        Args:
            X: PCA-reduced feature matrix (n_samples, n_components).

        Returns:
            np.ndarray: Cluster labels for each sample.

        Raises:
            ValueError: If input array is invalid.
        """
        if X.size == 0:
            raise ValueError("Input array cannot be empty")

        best_k = None
        best_score = float("-inf")

        for k in HDBSCAN_K_RANGE:
            try:
                clustering = HDBSCAN(min_cluster_size=k, min_samples=HDBSCAN_MIN_SAMPLES)
                labels = clustering.fit_predict(X)
                unique_labels = set(labels)
                if len(unique_labels) <= 1 or len(unique_labels) >= len(X):
                    continue

                score = silhouette_score(X, labels, metric=DISTANCE_METRIC)
                if score > best_score:
                    best_score = score
                    best_k = k
            except Exception:
                continue

        if best_k is None:
            best_k = min(HDBSCAN_K_RANGE)
            logger.warning(
                "No valid HDBSCAN configuration found; using default min_cluster_size=%s", best_k
            )

        final_clustering = HDBSCAN(min_cluster_size=best_k, min_samples=HDBSCAN_MIN_SAMPLES)
        return final_clustering.fit_predict(X)

    @staticmethod
    def k_medoids(
        distances: np.ndarray,
        cluster_label: int,
        all_labels: np.ndarray,
        k: int = KMEDOIDS_DEFAULT_K,
    ) -> Tuple[np.ndarray, float]:
        """Select representative medoids for a cluster.

        Args:
            distances: Pairwise distance matrix for all samples.
            cluster_label: The cluster identifier to process.
            all_labels: Full label array from HDBSCAN.
            k: Number of medoids to return.

        Returns:
            Tuple[np.ndarray, float]: Indices of medoids and total medoid cost.
        """
        cluster_indices = np.where(all_labels == cluster_label)[0]
        if cluster_indices.size == 0:
            return np.array([], dtype=int), float("nan")

        if cluster_indices.size < KMEDOIDS_MIN_CLUSTER_SIZE:
            logger.debug(
                "Cluster %s size %s is too small for k-medoids; returning all samples",
                cluster_label,
                cluster_indices.size,
            )
            return cluster_indices, float("nan")

        cluster_distances = distances[np.ix_(cluster_indices, cluster_indices)]
        costs = cluster_distances.sum(axis=1)
        best_indices = np.argsort(costs)[: min(k, cluster_indices.size)]
        medoids = cluster_indices[best_indices]
        return medoids, float(costs[best_indices[0]])

    @staticmethod
    def _extract_structure_features(pdb_files: List[Path]) -> np.ndarray:
        """Extract low-dimensional structural features from PDB files."""
        rows = []
        for pdb_file in pdb_files:
            universe = mda.Universe(str(pdb_file))
            backbone = universe.select_atoms("name CA")
            if len(backbone) == 0:
                backbone = universe.select_atoms("all")

            coords = backbone.positions
            centroid = coords.mean(axis=0) if coords.size else np.zeros(3, dtype=float)
            distances = (
                np.linalg.norm(coords - centroid, axis=1)
                if coords.size
                else np.zeros(1, dtype=float)
            )

            rows.append(
                [
                    float(len(coords)),
                    float(coords.shape[0]),
                    float(distances.mean() if distances.size else 0.0),
                    float(distances.std() if distances.size else 0.0),
                    float(np.linalg.norm(centroid)),
                ]
            )

        return np.vstack(rows) if rows else np.empty((0, 5), dtype=float)

    @staticmethod
    def _save_cluster_plot(X: np.ndarray, labels: np.ndarray, output_path: Path) -> None:
        """Save a 2D scatter plot for clusters."""
        fig, ax = plt.subplots(figsize=(8, 6))
        unique_labels = sorted(set(labels))
        for label in unique_labels:
            mask = labels == label
            label_text = "noise" if label == -1 else f"cluster_{label}"
            ax.scatter(
                X[mask, 0],
                X[mask, 1],
                label=label_text,
                alpha=0.7,
                edgecolors="w",
                s=60,
            )

        ax.set_title(f"Blind screening cluster map for {output_path.stem}")
        ax.set_xlabel("PC 1")
        ax.set_ylabel("PC 2")
        ax.legend(loc="best", fontsize="small")
        fig.tight_layout()
        fig.savefig(output_path)
        plt.close(fig)

    @staticmethod
    def _save_cluster_summary(
        output_path: Path,
        pdb_files: List[Path],
        labels: np.ndarray,
        medoid_indices: List[int],
    ) -> None:
        """Persist cluster assignments and medoid summaries to CSV."""
        with output_path.open("w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(["cluster", "file_path", "is_medoid"])
            medoid_set = set(medoid_indices)
            for index, pdb_file in enumerate(pdb_files):
                writer.writerow(
                    [
                        int(labels[index]),
                        str(pdb_file.relative_to(output_path.parent)),
                        int(index in medoid_set),
                    ]
                )

    @staticmethod
    def _save_best_hits(output_path: Path, medoids: dict) -> None:
        """Persist medoid structure best hits to a plain text file."""
        with output_path.open("w") as handle:
            for cluster_label, cluster_medoids in sorted(medoids.items()):
                handle.write(f"Cluster {cluster_label}:\n")
                for medoid_path in cluster_medoids:
                    handle.write(f"  {medoid_path}\n")
                handle.write("\n")

    def __init__(self, pdb1_name: str, blind_path: str) -> None:
        """Initialize blind screening analysis."""
        self.pdb1_name = pdb1_name
        self.blind_path = Path(blind_path)

        if not self.blind_path.exists():
            raise FileNotFoundError(f"Blind screening path not found: {blind_path}")

        logger.info("Starting blind screening for %s", pdb1_name)

        self._setup_foldseek_database()
        self._compute_similarities()
        self._perform_clustering_analysis()
        logger.info("Blind screening completed successfully")

    def _setup_foldseek_database(self) -> None:
        """Create Foldseek database from PDB files."""
        pdb_files = sorted(self.blind_path.glob("*.pdb"))
        if not pdb_files:
            raise FileNotFoundError(f"No PDB files found in {self.blind_path}")

        db_path = self.blind_path / "blind_db"
        cmd = ["foldseek", "createdb"] + [str(pdb_file) for pdb_file in pdb_files] + [str(db_path)]

        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"Foldseek createdb failed: {result.stderr.strip()}")

        logger.info("Created Foldseek database at %s", db_path)

    def _compute_similarities(self) -> None:
        """Compute pairwise structure similarities using Foldseek."""
        db_path = self.blind_path / "blind_db"
        result_path = self.blind_path / "blind_result"
        temp_dir = self.blind_path / "blind_tmp"
        temp_dir.mkdir(parents=True, exist_ok=True)

        cmd = ["foldseek", "search", str(db_path), str(db_path), str(result_path), str(temp_dir)]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"Foldseek search failed: {result.stderr.strip()}")

        logger.info("Computed structure similarities")

    def _perform_clustering_analysis(self) -> None:
        """Perform PCA, HDBSCAN clustering, and k-medoids selection."""
        pdb_files = sorted(self.blind_path.rglob("*.pdb"))
        if not pdb_files:
            raise FileNotFoundError("No PDB files found for clustering")

        features = self._extract_structure_features(pdb_files)
        if features.shape[0] < 2:
            raise ValueError("At least two structures are required for clustering")

        reduced = PCA(n_components=2).fit_transform(features)
        labels = self.cluster_structures(reduced)

        distance_matrix = distance.cdist(reduced, reduced, metric=DISTANCE_METRIC)
        medoid_map = {}
        medoid_indices = []
        for cluster_label in sorted(set(labels)):
            if cluster_label == -1:
                continue

            medoids, cost = self.k_medoids(distance_matrix, cluster_label, labels)
            medoid_indices.extend(medoids.tolist())
            medoid_map[cluster_label] = [
                str(pdb_files[i].relative_to(self.blind_path)) for i in medoids
            ]
            logger.info("Cluster %s medoids: %s", cluster_label, medoid_map[cluster_label])

        summary_file = self.blind_path / f"{self.pdb1_name}_structures_of_interest.csv"
        self._save_cluster_summary(summary_file, pdb_files, labels, medoid_indices)

        best_hits_file = self.blind_path / f"{self.pdb1_name}_best_hits.txt"
        self._save_best_hits(best_hits_file, medoid_map)

        plot_file = self.blind_path / f"{self.pdb1_name}_blind_clusters.png"
        self._save_cluster_plot(reduced, labels, plot_file)

        logger.info("Saved clustering summary to %s", summary_file)
        logger.info("Saved best hits to %s", best_hits_file)
        logger.info("Saved cluster plot to %s", plot_file)
