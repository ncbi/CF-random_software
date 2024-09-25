# Data and code for CF-random
The general installation and usage guidance of CF-random for predicting the alternative conformation and fold-switching proteins.


# Installation
  1.
  2.
  3.


# Usage
CF-random provides the different prediction modes such as fold-switching default, alternative conformation, and blind mode.
To execute all modes of CF-random, a multiple sequence alignment (MSA) is required. PDB IDs for both fold1 (dominant conformation) and fold2 (alternative conformation) are required for running the default fold-switching and alternative conformation.

For running the fold-switching default mode
  > python main.py --fname folder_containing_MSA/ --pdb1 fold1.pdb --pdb2 fold2.pdb --option FS

For executing the alternative confroamtion mode
  > python main.py --fname folder_containing_MSA/ --pdb1 fold1.pdb --pdb2 fold2.pdb --option AC

For running the CF-random with blind mode covering both fold-switching and alternative conformation
  > python main.py --fname folder_containing_MSA/ --option blind


