# Data and code for CF-random
The general installation and usage guidance of CF-random for predicting the alternative conformation and fold-switching proteins.<br>
For using CF-random in Colab notebook, please use following [link](https://colab.research.google.com/drive/1LsSFe8FxJaLfNGUcE5HMgxxwGGlLfexk?usp=sharing).<br>


# Installation
Local CF-random uses the [localcolabfold](https://github.com/YoshitakaMo/localcolabfold) and [Foldseek](https://github.com/steineggerlab/foldseek) under linux environment.<br>
For more details about localcolabfold, please visit [here.](https://github.com/YoshitakaMo/localcolabfold) <br>
We currently not support the Windows and MacOS environment.<br>

Installation process including localcolabfold, dependencies, and Foldseek is done with following command.
```
wget https://raw.githubusercontent.com/YoshitakaMo/localcolabfold/main/install_colabbatch_linux.sh
bash install_colabbatch_linux.sh
```
<br>


After finishing the installation of localcolabfold, activate conda environment.<br>
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
CF-random provides the different prediction modes such as fold-switching default, alternative conformation, and blind mode.<br>
To execute all modes of CF-random, a multiple sequence alignment (MSA) is required (Essential).<br> PDB files for both fold1 (dominant conformation) and fold2 (alternative conformation) are required for TM-score measurement with reference files. Blind mode doesn't require the PDB files, but default fold-switching and alternative conformation modes require the PDB.<br>

```
 -fname ####    |  folder name having a multiple sequence alignment (MSA)
 -pdb1  ####    |  dmoninat reference model used to calculate TM-score with predicted models
 -pdb2  ####    |  alternative refernece model used to calculate TM-score with predicted models
 -option ###    |  AC: predicting alternative conformations of protein with references, FS: predicting the fold-switching protein with references, and blind: predicting the alternative conformations or fold-switching proteins without reference PDB files.
```
*In default mode (fold-switching and alternative conformation), CF-ramdon produces the results of TM-scores (csv and png files), plDDT, and information of selected random MSA. If CF-random predicts the both folds, generated prediction files are deposited under successed_prediction/pdb1_name and additional_sampling/pdb1_name . If not, it would not generate anything. <br>
*Before running the default mode of fold-switching, setting the "range_fs_pairs_all.txt" file is required. The name of reference PDB files, residue ranges of reference pdb files, and residue ranges of prediction files. ColabFold generates the residue index starting from 1, so please choose the residue range of fold-switching region correctly. CF-random reads the residue index in PDB file, make sure that selection of residue range is correct. <br>
 examples) pdb1, pdb2, XXX-XXX, XXX-XXX, XXX-XXX, XXX-XXX <br>
*In blind mode, predicted files are deposited under blind_prediction/pdb1_name . CF-random with blind mode produces the comparison result with Foldseek. <br><br>


Before running the CF-random, please check out your conda environment.<br>
The name of conda environment would be the installed directory.
```
conda activate <your localfolabfold environment>
** examples of generated localcolabfold name would be "directory of installed the localcolabfold"/localcolabfold/colabfold-conda
```
<br>

### For running the fold-switching default mode. <br>
```
python main.py --fname folder_containing_MSA/ --pdb1 fold1.pdb --pdb2 fold2.pdb --option FS
```
<br>

### For executing the alternative confroamtion mode. <br>
```
python main.py --fname folder_containing_MSA/ --pdb1 fold1.pdb --pdb2 fold2.pdb --option AC
```
<br>

### For running the CF-random with blind mode covering both fold-switching and alternative conformation. <br>
```
python main.py --fname folder_containing_MSA/ --option blind
```
<br>

# Examples
To be updated
<br><br>

# How to Cite
To be updated
<br><br>

# Lincese
To be updated


