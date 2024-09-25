# Data and code for CF-random
The general installation and usage guidance of CF-random for predicting the alternative conformation and fold-switching proteins.


# Installation
CF-random uses the [localcolabfold](https://github.com/YoshitakaMo/localcolabfold) and [Foldseek](https://github.com/steineggerlab/foldseek). It can be simply installed with following commands. <br><br>
1. Install miniconda <br>
``bash install_colabbatch_linux_custom.sh``

2. Download Foldseek installation file and install <br>
``wget https://mmseqs.com/foldseek/foldseek-linux-avx2.tar.gz``<br>
`tar xvzf foldseek-linux-avx2.tar.gz`
``export PATH=$(pwd)/foldseek/bin/:$PATH``


```wget https://mmseqs.com/foldseek/foldseek-linux-avx2.tar.gz
tar xvzf foldseek-linux-avx2.tar.gz
export PATH=$(pwd)/foldseek/bin/:$PATH```



# Usage
CF-random provides the different prediction modes such as fold-switching default, alternative conformation, and blind mode.<br>
To execute all modes of CF-random, a multiple sequence alignment (MSA) is required. PDB IDs for both fold1 (dominant conformation) and fold2 (alternative conformation) are required for running the default fold-switching and alternative conformation.<br>

For running the fold-switching default mode. <br>
``python main.py --fname folder_containing_MSA/ --pdb1 fold1.pdb --pdb2 fold2.pdb --option FS``<br>

For executing the alternative confroamtion mode. <br>
``python main.py --fname folder_containing_MSA/ --pdb1 fold1.pdb --pdb2 fold2.pdb --option AC``<br>

For running the CF-random with blind mode covering both fold-switching and alternative conformation. <br>
``python main.py --fname folder_containing_MSA/ --option blind``<br>


