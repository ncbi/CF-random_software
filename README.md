# Data and code for CF-random

General installation and usage guidance of CF-random for predicting the alternative conformation and fold-switching proteins.

To run CF-random in a Colab notebook, please use following [link](https://colab.research.google.com/drive/16pD2tUMkUx1gwDxZXcSr9WOosYp0ZU6j?authuser=0).

<a target="_blank" href="https://colab.research.google.com/drive/16pD2tUMkUx1gwDxZXcSr9WOosYp0ZU6j?authuser=0">
 <img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open CF-random Colab"/>
</a>

# Installation

CF-random depends on ColabFold (structure prediction) and Foldseek (structure search). We provide a convenience script `install.sh` to set up the conda environment and common packages.

Quick start (recommended):

Run the installer (creates/activates a `cf-random` conda env if needed):

```bash
bash install.sh
```


# Usage
* CF-random supports three main modes: fold-switching (FS), alternative conformation (AC), and blind mode.
* For FS and AC modes you should provide two reference PDBs (fold1 and fold2). Blind mode does not require reference PDBs.
* Predictions use MSAs; provide an MSA directory per target. Do not reuse the same output folder for multiple runs.
* PDB files should contain a single chain; multi-chain PDBs may be converted automatically in some workflows, but providing single-chain PDBs avoids issues.

```
 --fname ####    |  folder name having a multiple sequence alignment (MSA)
 --pname ####    |  project name for running blind mode (only for blind mode)
 --pdb1  ####    |  dominant reference model used to calculate TM-score with predicted models
 --pdb2  ####    |  alternative reference model used to calculate TM-score with predicted models
 --nMSA  ####    |  the number of additional samples for predicting the structure with MSAs, default = 0
 --type  ####    |  can choose the model type of Colabfold. e.g. ptm, monomer, and multimer
 --options ###   |  AC: predicting alternative conformations of protein with references, FS: predicting the fold-switching protein with references, and blind: predicting the alternative conformations or fold-switching proteins without reference PDB files.
```
* Output: TM-score CSV/PNG, plDDT values, and selected MSA info. Successful predictions are saved under `successed_prediction/<pdb1_name>/`.
* For FS runs you must provide `range_fs_pairs_all.txt` describing the FS region ranges. ColabFold uses 1-based residue indexing; ensure ranges match your PDB/sequence.
* --nMSA can be applied for all options, but --nESN cannot be used for blind mode.
* In blind mode, predicted files are deposited under blind_prediction/pdb1_name . CF-random with blind mode produces the comparison result with Foldseek. <br>
* ### For running the foldseek in blind mode, Foldseek parameter files and running Python scripts should be in same directory. <br>

* Before running the CF-random, ensure that the CF-random conda environment is activated:<br>
```
conda activate CF-random
```
<br>

# Examples
We provide some examples how users can run the CF-random with different modes.<br>
First two modes such as fold-switching and alternative conformation are default modes of CF-random and the last one is a blind mode.
## 1. For CF-random with fold-switching mode. <br>
For this example, RfaH would be predicted with two reference structures (i.e., 2oug_C.pdb and 6c6s_D.pdb) and a MSA file.
```
python main.py --fname 2oug_C-search/ --pdb1 2oug_C.pdb --pdb2 6c6s_D.pdb --option FS
```
### Used input files: <br>
* PDB1: 2oug_C.pdb <br>
* PDB2: 6c6s_D.pdb <br>
* MSA: 2oug_C-search/0.a3m (MSA file should be in a folder) <br>
* range_fs_pairs_all.txt (This file is required for reading the fold-switching region in refernece and predicted structures. Users should check the region before running this mode.) <br>

*This takes <30 Minutes to run on an A100 GPU (generates 200 structures total).* <br>

### Generated output files: <br>
_Predicted files from deep and random MSAs are deposited in 'predictions_all' directory._ <br>
_If CF-random fails to find the selected random MSA, all generated files will be in 'predictions_all' directory._ <br>
* TM-score plot of whole structure: TMscore_fs-region_full-MSA_2oug_C.png <br>
* TM-score plot of fold-switching region: TMscore_full-MSA_2oug_C.png <br>
* TM-score plot of fold-switching region with label of prediction rank: TMscore_fs-region_full-MSA_2oug_C_label.png <br> 
* TM-scores and plDDT scores of predictions with deep MSA: TMs_plDDT_full_all_2oug_C.csv <br>
* TM-scores and plDDT scores of predictions with random MSAs: TMs_plDDT_rand_all_2oug_C.csv <br>
* Selection of random MSA: selected_MSA-size_2oug_C.csv (When CF-random finds the MSA depth)
  - MSA depth information (e.g. # = max-seq:max-seq-extra) (0 = 1:2, 1 = 2:4, 2 = 4:8, 3 = 8:16, 4 = 16:32, 5 = 32:64, 6 = 64:128) <br>


## 2. For CF-random with alternative conformation mode. <br>
For this mode, Lactococcal OppA would be predicted with two reference structures (i.e., 3drh.pdb and 3drf.pdb) and an MSA file. <br>
```
python main.py --fname 5olw_A-search --pdb1 5olw_A.pdb --pdb2 5olx_A.pdb --option AC --nMSA 5
```
### Used input files: <br>
* PDB1: 5olw_A.pdb <br>
* PDB2: 5olx_A.pdb <br>
* MSA: 5olw_A-search/0.a3m (MSA file should be in a folder) <br>

*This takes <70 Minutes to run on an A100 GPU (generates 200 structures total; protein is large: ~250 residues).* <br>

### Generated output files: <br>
_Predicted files from deep and random MSAs are deposited in 'predictions_all' directory._ <br>
_If CF-random fails to find the selected random MSA, all generated files will be in 'predictions_all' directory._ <br>
* TM-score plot of whole structure: TMscore_full-MSA_5olw_A.png <br>
* TM-scores and plDDT scores of predictions with deep MSA: TMs_plDDT_full_all_5olw_A.csv <br>
* TM-scores and plDDT scores of predictions with random MSAs: TMs_plDDT_rand_all_5olw_A.csv <br>
* Selection of random MSA: selected_MSA-size_3drh_A.csv (When CF-random finds the MSA depth)
  - MSA depth information (e.g. # = max-seq:max-seq-extra) (0 = 1:2, 1 = 2:4, 2 = 4:8, 3 = 8:16, 4 = 16:32, 5 = 32:64, 6 = 64:128) <br>

## 3. For CF-random with blind mode covering both fold-switching and alternative conformation. <br>
```
python main.py --pname Mad2_test --fname 2vfx_L-search/ --option blind
```

*Before running this code, make a symbolic link to the foldseek pdb libraries in the directory where you run the command above.*

### Used input files: <br>
MSA: 2vfx_L-search/0.a3m (MSA file should be in a folder) <br>


### Generated output files: <br>
_Predicted files from deep and random MSAs are deposited in 'blind_prediction' directory._ <br>
_If user uses the option '--pname', the name of output files would be entered '--pname'._<br>
* List of prediction files: Mad2-structures_of_interest.csv
* The best hit list of alternative conformations: Mad2-structures_of_interest.csv
* Cluster analysis result as an image file: Mad2-cluster.png

*This takes <70 Minutes to run on an A100 GPU (generates 200 structures total + 200 foldseek files).* <br>

# How to Cite
Lee, M., Schafer, J.W., Prabakaran, J. et al. Large-scale predictions of alternative protein conformations by AlphaFold2-based sequence association. Nat Commun 16, 5622 (2025). https://doi.org/10.1038/s41467-025-60759-5
<br><br>

# License
Please see the LICENSE.md file.


