# CF-random

CF-random predicts alternative protein conformations and fold-switching behavior using AlphaFold2 with varied MSA depths.

To run CF-random in a Colab notebook:

<a target="_blank" href="https://colab.research.google.com/drive/16pD2tUMkUx1gwDxZXcSr9WOosYp0ZU6j?authuser=0">
  <img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open CF-random Colab"/>
</a>

---

## Installation
CF-random depends on ColabFold (structure prediction) and Foldseek (structure search). A convenience script sets up the conda environment and required packages:

```bash
chmod +x install.sh
./install.sh
```

---

## Usage

CF-random supports three modes:

- **FS** — fold-switching prediction with two reference PDBs
- **AC** — alternative conformation prediction with two reference PDBs
- **blind** — alternative conformation or fold-switching without reference PDBs

### General notes

- FS and AC modes require two reference PDB files (fold1 and fold2).
- Each target needs its own MSA directory; do not reuse output folders across runs.
- All required PDB files, MSA files, and Python scripts must be in the same directory.
- PDB files should contain a single chain. Multi-chain PDBs may be converted automatically in some workflows, but providing single-chain PDBs avoids issues.
- FS mode requires a `range_fs_pairs_all.txt` file describing the fold-switching region (see [Fold-switching mode](#1-fold-switching-mode-fs) below for format details). ColabFold uses 1-based residue indexing; ensure ranges match your PDB/sequence.
- `--num_msa` and `--num_ens` apply to all modes except blind, which does not support `--num_ens`.
- Activate the conda environment before running:

```bash
conda activate cf-random
```

### Arguments

| Argument | Required | Description |
|---|---|---|
| `--option` | Yes | Run mode: `AC`, `FS`, `inAC`, or `blind` |
| `--fname` | Yes | MSA folder name (output of ColabSearch) |
| `--pdb1` | FS / AC | Target crystal structure PDB |
| `--pdb2` | FS / AC | Alternative crystal structure PDB |
| `--fmname` | No | Multimer MSA folder name (output of ColabSearch) |
| `--pname` | blind | Job name for blind mode output naming |
| `--num_msa` | No | Number of additional MSA seeds to run, added to the default 5 |
| `--num_ens` | No | Number of ensemble samples to generate |
| `--model_type` | No | ColabFold model type: `ptm`, `monomer`, or `multimer` |

---

## Examples

### 1. Fold-switching mode (FS)

Predicts RfaH using two reference structures and an MSA.

```bash
cf-random --fname 2oug_C-search/ --pdb1 2oug_C.pdb --pdb2 6c6s_D.pdb --option FS
```

**Input files:**
- `2oug_C.pdb` — dominant reference
- `6c6s_D.pdb` — alternative reference
- `2oug_C-search/0.a3m` — MSA
- `range_fs_pairs_all.txt` — fold-switching region definition

**`range_fs_pairs_all.txt` format:**

Each line defines the fold-switching region for a pair of reference structures:
pdb1, pdb2, XXX-XXX, XXX-XXX, XXX-XXX, XXX-XXX

Fields are: PDB1 name, PDB2 name, residue range of reference 1, residue range of reference 2, residue range of prediction 1, residue range of prediction 2. ColabFold generates residue indices starting from 1; ensure your ranges match the residue numbering in your PDB files.

*Generates 200 structures; takes under 30 minutes on an A100 GPU.*

**Output files** (written to `predictions_all/`):

| File | Description |
|---|---|
| `TMscore_full-MSA_2oug_C.png` | TM-score scatter plot, whole structure |
| `TMscore_fs-region_full-MSA_2oug_C.png` | TM-score scatter plot, fold-switching region |
| `TMscore_fs-region_full-MSA_2oug_C_label.png` | Same plot with prediction rank labels |
| `TMs_plDDT_full_all_2oug_C.csv` | TM-scores and pLDDT for deep MSA predictions |
| `TMs_plDDT_rand_all_2oug_C.csv` | TM-scores and pLDDT for random MSA predictions |
| `selected_MSA-size_2oug_C.csv` | Selected MSA depth (when CF-random identifies the optimal depth) |

MSA depth encoding: `0=1:2`, `1=2:4`, `2=4:8`, `3=8:16`, `4=16:32`, `5=32:64`, `6=64:128`

---

### 2. Alternative conformation mode (AC)

Predicts Lactococcal OppA using two reference structures and an MSA.

```bash
cf-random --fname P71447-search/ --pdb1 2wfa_A.pdb --pdb2 2wf5_A.pdb --option AC --num_msa 5
```

**Input files:**
- `2wfa_A.pdb` — dominant reference
- `2wf5_A.pdb` — alternative reference
- `P71447-search/P71447_converted.a3m` — MSA

*Generates 200 structures; takes under 70 minutes on an A100 GPU (~250 residue protein).*

**Output files** (written to `predictions_all/`):

| File | Description |
|---|---|
| `TMscore_full-MSA_2wfa_A.png` | TM-score scatter plot, whole structure |
| `TMs_plDDT_full_all_2wfa_A.csv` | TM-scores and pLDDT for deep MSA predictions |
| `TMs_plDDT_rand_all_2wfa_A.csv` | TM-scores and pLDDT for random MSA predictions |
| `selected_MSA-size_2wfa_A.csv` | Selected MSA depth (when CF-random identifies the optimal depth) |

MSA depth encoding: `0=1:2`, `1=2:4`, `2=4:8`, `3=8:16`, `4=16:32`, `5=32:64`, `6=64:128`

---

### 3. Blind mode

Predicts alternative conformations or fold-switching without reference PDBs. Uses Foldseek to cluster predicted structures and identify structures of interest.

```bash
cf-random --pname Mad2_test --fname 2vfx_L-search/ --option blind
```

> **Note:** Before running blind mode, make a symbolic link to the Foldseek PDB libraries in the directory where you run the command.

**Input files:**
- `2vfx_L-search/0.a3m` — MSA

*Generates 200 structures + 200 Foldseek result files; takes under 70 minutes on an A100 GPU.*

**Output files** (written to `blind_prediction/Mad2_test/`):

| File | Description |
|---|---|
| `Mad2_test-cluster.png` | PCA + HDBSCAN cluster plot |
| `Mad2_test-structures_of_interest.csv` | Representative structures per cluster |
| `Mad2_test-structures_of_interest.pse` | PyMOL session with aligned representatives |

---

## Citation

Lee, M., Schafer, J.W., Prabakaran, J. et al. Large-scale predictions of alternative protein conformations by AlphaFold2-based sequence association. *Nat Commun* **16**, 5622 (2025). https://doi.org/10.1038/s41467-025-60759-5

## License

See [LICENSE.md](LICENSE.md).