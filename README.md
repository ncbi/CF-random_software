# Data and code for CF-random
General installation and usage guidance of CF-random for predicting the alternative conformation and fold-switching proteins.<br>
To use the CF-random in Colab notebook, please use following [link](https://colab.research.google.com/drive/1LsSFe8FxJaLfNGUcE5HMgxxwGGlLfexk?usp=sharing).<br><br>
<a target="_blank" href="https://colab.research.google.com/drive/1LsSFe8FxJaLfNGUcE5HMgxxwGGlLfexk?usp=sharing">
 <img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open CF-random Colab"/>
</a>


# Installation
CF-random uses the [localcolabfold](https://github.com/YoshitakaMo/localcolabfold) and [Foldseek](https://github.com/steineggerlab/foldseek) under linux environment.<br>
For more details about localcolabfold, please visit [here.](https://github.com/YoshitakaMo/localcolabfold) <br>
We currently not support the Windows and MacOS environment.<br>

Installation process including localcolabfold, dependencies, and Foldseek is done with following commands.
```
wget https://raw.githubusercontent.com/YoshitakaMo/localcolabfold/main/install_colabbatch_linux.sh
bash install_colabbatch_linux.sh

** A bash script in install folder 
bash install_colabbatch_linux.sh
```
<br>


After the installation of localcolabfold, activate conda environment.<br>
Conda name would be "installed localcolabfold directory", e.g.) "current directory" + /localcolabfold/colabfold-conda <br>
```
conda activate " current directory "/localcolabfold/colabfold-conda
pip install textalloc tmtools adjustText thefuzz mdtraj
pip3 install -U scikit-learn
```
Once finalize the dependencies, install the Foldseek.
<br>
```
conda install -c conda-forge -c bioconda foldseek
foldseek databases PDB pdb tmp
```



# Usage
* CF-random provides the different prediction modes such as fold-switching default, alternative conformation, and blind mode.<br>
* To execute all modes of CF-random, a multiple sequence alignment (MSA) is required (Essential). To avoid the overwriting the output files, we would recommend to use different folder name having MSA. <br>
* PDB files for both fold1 (dominant conformation) and fold2 (alternative conformation) are required for TM-score measurement with reference files. Blind mode doesn't require the PDB files, but default fold-switching and alternative conformation modes require the PDB.<br>
* ### All required PDB files and MSA file should be in same directory with provided Python scripts.
* Please make sure that a PDB file should have a single chain, not multiple chains. If PDB file has multiple chains, CF-random would be stopped. <Pbr>

```
 -fname ####    |  folder name having a multiple sequence alignment (MSA)
 -pname ####    |  project name for running the blind mode (only for blind mode)
 -pdb1  ####    |  dmoninat reference model used to calculate TM-score with predicted models
 -pdb2  ####    |  alternative refernece model used to calculate TM-score with predicted models
 -option ###    |  AC: predicting alternative conformations of protein with references, FS: predicting the fold-switching protein with references, and blind: predicting the alternative conformations or fold-switching proteins without reference PDB files.
```
* In default mode (fold-switching and alternative conformation), CF-ramdon produces the results of TM-scores (csv and png files), plDDT, and information of selected random MSA. If CF-random predicts the both folds, generated prediction files are deposited under successed_prediction/pdb1_name and additional_sampling/pdb1_name . If not, it would not generate anything. <br>
* Before running the default mode of fold-switching, setting the "range_fs_pairs_all.txt" file is required. The name of reference PDB files, residue ranges of reference pdb files, and residue ranges of prediction files. ColabFold generates the residue index starting from 1, so please choose the residue range of fold-switching region correctly. CF-random reads the residue index in PDB file, make sure that selection of residue range is correct. <br>
 examples) pdb1, pdb2, XXX-XXX, XXX-XXX, XXX-XXX, XXX-XXX <br>
* In blind mode, predicted files are deposited under blind_prediction/pdb1_name . CF-random with blind mode produces the comparison result with Foldseek. <br><br>
* ### For running the foldseek in blind mode, Foldseek parameter files and running Python scripts should be in same directory. <br>

Before running the CF-random, please check out your conda environment.<br>
The name of conda environment would be the installed directory.
```
conda activate <your localfolabfold environment>
** examples of generated localcolabfold name would be "directory of installed the localcolabfold"/localcolabfold/colabfold-conda <br>
```
<br>

# Examples
We provide some examples how users can run the CF-random with different modes.<br>
First two modes such as fold-switching and alternative conformation are default modes of CF-random and the last one is a blind mode.
### For CF-random with fold-switching mode. <br>
For this example, RfaH would be predicted with two reference structures (i.e., 2oug_C.pdb and 6c6s_D.pdb) and a MSA file.
```
python main.py --fname 2oug_C-search/ --pdb1 2oug_C.pdb --pdb2 6c6s_D.pdb --option FS
```
* Used input files: <br>
PDB1: 2oug_C.pdb <br>
PDB2: 6c6s_D.pdb <br>
MSA: 2oug_C-search/0.a3m (MSA file should be in a folder) <br>
range_fs_pairs_all.txt (This file is required for reading the fold-switching region in refernece and predicted structures. Users should check the region before running this mode.) <br>

* Generated output files:
TM-score plot of whole structure: TMscore_fs-region_full-MSA_2oug_C.png <br>
TM-score plot of fold-switching region: TMscore_full-MSA_2oug_C.png <br>
TM-score plot of fold-switching region with label of prediction rank: TMscore_fs-region_full-MSA_2oug_C_label.png <br> 
plDDT scores of predictions with deep MSA plddt_full-MSA_2oug_C.csv <br>
plDDT scores of predictions with random MSAs plddt_random-MSA_2oug_C.csv <br>
plDDT scores of ensemble generation: plddt_additional-MSA_2oug_C.csv <br>


### For CF-random with alternative confroamtion mode. <br>
For this mode, Phosphoglycerate kinase 1 (PGK1) would be predicted with two reference structures (i.e., 2ybe.pdb and 2xe7.pdb) and a MSA file.
```
python main.py --fname 2ybe_A-search/ --pdb1 2ybe_A.pdb --pdb2 2xe7_A.pdb --option AC
```
* Used input files: <br>
PDB1: 2ybe_A.pdb <br>
PDB2: 2xe7_A.pdb <br>
MSA: 2ybe_A-search/0.a3m (MSA file should be in a folder) <br>

* Generated output files:

### For CF-random with blind mode covering both fold-switching and alternative conformation. <br>
```
<option 1>
python main.py --fname folder_containing_MSA/ --option blind
<option 2>
python main.py --fname folder_containing_MSA/ --pname Mad2 --option blind
```
* Used input files: <br>
MSA: 2vfx_L-search/0.a3m (MSA file should be in a folder) <br>

* Generated output files:

# How to Cite
To be updated
<br><br>

# License
To be updated


