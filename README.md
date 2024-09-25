# Data and code for CF-random
The general installation and usage guidance of CF-random for predicting the alternative conformation and fold-switching proteins.<br><br>

# Installation
CF-random uses the [localcolabfold](https://github.com/YoshitakaMo/localcolabfold) and [Foldseek](https://github.com/steineggerlab/foldseek).
It can be simply installed with following commands. <br>

1. Download and install localcolabfold <br>
```
'''
Download the bash script file in installation folder or use following command
'''
wget https://raw.githubusercontent.com/YoshitakaMo/localcolabfold/main/install_colabbatch_linux.sh
bash install_colabbatch_linux_custom.sh
```

2. Download Foldseek installation file and install <br>
```
wget https://mmseqs.com/foldseek/foldseek-linux-avx2.tar.gz
tar xvzf foldseek-linux-avx2.tar.gz
export PATH=$(pwd)/foldseek/bin/:$PATH
```

3. Add additional libries to conda environment<br>
Please check out your miniconda environment after finished the first step.<br>
The name of miniconda environment would be the installed directory.
```
conda activate <your localfolabfold environment>
conda env update --file CF-random.yml --prune
```
<br>

# Usage
CF-random provides the different prediction modes such as fold-switching default, alternative conformation, and blind mode.<br>
To execute all modes of CF-random, a multiple sequence alignment (MSA) is required.<br> PDB files for both fold1 (dominant conformation) and fold2 (alternative conformation) are required for running the default fold-switching and alternative conformation. Blind mode doesn't require the PDB files.<br>

For running the fold-switching default mode. <br>
```
python main.py --fname <folder_containing_MSA> --pdb1 fold1.pdb --pdb2 fold2.pdb --option FS
```

For executing the alternative confroamtion mode. <br>
```
python main.py --fname <folder_containing_MSA> --pdb1 fold1.pdb --pdb2 fold2.pdb --option AC
```

For running the CF-random with blind mode covering both fold-switching and alternative conformation. <br>
```
python main.py --fname <folder_containing_MSA> --option blind
```
<br>

# How to Cite
To be updated
<br>

# Lincese
To be updated


