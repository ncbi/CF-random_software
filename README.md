# Data and code for CF-random
General installation and usage guidance of CF-random for predicting the alternative conformation and fold-switching proteins.<br>
To run CF-random in a Colab notebook, please use following [link](https://colab.research.google.com/drive/1LsSFe8FxJaLfNGUcE5HMgxxwGGlLfexk?usp=sharing).<br><br>
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

** Or use a bash script in install folder 
bash install_colabbatch_linux.sh
```
<br>


After the installation of localcolabfold, add the localcolabfold path to your .bashrc file:.<br>
```
export PATH="/path/to/your/localcolabfold/colabfold-conda/bin:$PATH"
```
<br>

Then reactivate your .bashrc file <br>

Now create a conda new conda environment: 
```
conda create --name CF-random python=3.10
conda activate CF-random
pip install textalloc tmtools adjustText thefuzz mdtraj biopython
pip3 install -U scikit-learn
```
Once the dependencies are installed, install Foldseek.
<br>
```
conda install -c conda-forge -c bioconda foldseek
foldseek databases PDB pdb tmp
```
<br>

### We recommend running the foldseek databases command in a directory where the libraries can be stored. <br>


# Usage
* CF-random has different prediction modes such as fold-switching default, alternative conformation, and blind mode.<br>
* To execute all modes of CF-random, a multiple sequence alignment (MSA) is required. To avoid the overwriting the output files, we recommend using a different folder containing MSA. <br>
* PDB files for both fold1 (dominant conformation) and fold2 (alternative conformation) are required for TM-score measurement with reference files. Blind mode doesn't require PDB files, but default fold-switching and alternative conformation modes do.<br>
* ### All required PDB files and MSA file should be in same directory with provided Python scripts.
* Please make sure that a PDB file should have a single chain, not multiple chains. If PDB file has multiple chains, CF-random would be stopped. <Pbr>

```
 --fname ####    |  folder name having a multiple sequence alignment (MSA)
 --pname ####    |  project name for running blind mode (only for blind mode)
 --pdb1  ####    |  dominant reference model used to calculate TM-score with predicted models
 --pdb2  ####    |  alternative reference model used to calculate TM-score with predicted models
 --nMSA  ####    |  the number of additional samples for predicting the structure with MSAs, default = 0
 --nESN  ####    |  the number of additional samples for ensenble generation, default = 0
 --options ###    |  AC: predicting alternative conformations of protein with references, inAC: predicting the alternative conformation with the more number of samples, FS: predicting the fold-switching protein with references, and blind: predicting the alternative conformations or fold-switching proteins without reference PDB files.
```
* In default mode (fold-switching and alternative conformation), CF-ramdon produces the results of TM-scores (csv and png files), plDDT, and information of selected random MSA. If CF-random predicts the both folds, generated prediction files are deposited under successed_prediction/pdb1_name and additional_sampling/pdb1_name . If not, it would not generate anything. <br>
* Before running the default mode of fold-switching, setting the "range_fs_pairs_all.txt" file is required. The name of reference PDB files, residue ranges of reference pdb files, and residue ranges of prediction files. ColabFold generates the residue index starting from 1, so please choose the residue range of fold-switching region correctly. CF-random reads the residue index in PDB file, make sure that selection of residue range is correct. <br>
 examples) pdb1, pdb2, XXX-XXX, XXX-XXX, XXX-XXX, XXX-XXX <br>
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

*This takes <30 Minutes to run on an A100 GPU (generates 300 structures total).* <br>

### Generated output files: <br>
_Predicted files from deep and random MSAs are deposited in 'successed_prediction' directory, and ensembles were in 'additional_sampling' folder._ <br>
_If CF-random fails to find the selected random MSA, all generated files will be in 'failed_prediction' directory._ <br>
* TM-score plot of whole structure: TMscore_fs-region_full-MSA_2oug_C.png <br>
* TM-score plot of fold-switching region: TMscore_full-MSA_2oug_C.png <br>
* TM-score plot of fold-switching region with label of prediction rank: TMscore_fs-region_full-MSA_2oug_C_label.png <br> 
* TM-scores and plDDT scores of predictions with deep MSA: TMs_plDDT_full_all_2oug_C.csv <br>
* TM-scores and plDDT scores of predictions with random MSAs: TMs_plDDT_rand_all_2oug_C.csv <br>
* TM-scores and and plDDT scores of predictions with ensembles: TMs_plDDT_addi_all_2oug_C.csv
  - TM-scores of whole structure and fold-switching regions were saved in TMs_plDDT~ file with ensembles. <br>
* Selection of random MSA: selected_MSA-size_2oug_C.csv (When CF-random finds the MSA depth)
  - MSA depth information (e.g. # = max-seq:max-seq-extra) (0 = 1:2, 1 = 2:4, 2 = 4:8, 3 = 8:16, 4 = 16:32, 5 = 32:64, 6 = 64:128) <br>


## 2. For CF-random with alternative conformation mode. <br>
For this mode, Lactococcal OppA would be predicted with two reference structures (i.e., 3drh.pdb and 3drf.pdb) and an MSA file. <br>
```
python main.py --fname 5olw_A-search --pdb1 5olw_A.pdb --pdb2 5olx_A.pdb --option AC --nMSA 5 --nENS 5
```
### Used input files: <br>
* PDB1: 5olw_A.pdb <br>
* PDB2: 5olx_A.pdb <br>
* MSA: 5olw_A-search/0.a3m (MSA file should be in a folder) <br>

*This takes <70 Minutes to run on an A100 GPU (generates 175 structures total; protein is large: ~250 residues).* <br>

### Generated output files: <br>
_Predicted files from deep and random MSAs are deposited in 'successed_prediction' directory, and ensembles were in 'additional_sampling' folder._ <br>
_If CF-random fails to find the selected random MSA, all generated files will be in 'failed_prediction' directory._ <br>
* TM-score plot of whole structure: TMscore_full-MSA_5olw_A.png <br>
* TM-scores and plDDT scores of predictions with deep MSA: TMs_plDDT_full_all_5olw_A.csv <br>
* TM-scores and plDDT scores of predictions with random MSAs: TMs_plDDT_rand_all_5olw_A.csv <br>
* TM-scores and and plDDT scores of predictions with ensembles: TMs_plDDT_addi_all_5olw_A.csv
  - TM-scores of whole structure were saved in TMs_plDDT~ file with ensembles. <br>
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
* List of prediction files, foldseek searched pdb name, TM-score, and foldseek score: Mad2_test.csv
* Best hits of alternative conformations: Mad2_test_best_hits.txt
* Cluster analysis result as an image file: Mad2_test.png

*This takes <70 Minutes to run on an A100 GPU (generates 200 structures total + 200 foldseek files).* <br>

# How to Cite
To be updated
<br><br>

# License
To be updated


