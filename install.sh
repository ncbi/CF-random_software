#!/bin/bash
set -e

echo "=== CF-random Installation ==="

echo "[1/3] Creating conda environment..."
conda env create -f environment.yml
eval "$(conda shell.bash hook)"
conda activate cf-random

echo "[2/3] Downloading foldseek PDB database..."
mkdir -p pdb_database
foldseek databases PDB pdb_database/pdb pdb_database/tmp
rm -rf pdb_database/tmp

echo "[3/3] Verifying installation..."
python -c "from Bio.Data import SCOPData; print('  biopython ok')"
python -c "import numpy; print(f'  numpy {numpy.__version__}')"
python -c "import pandas; print(f'  pandas {pandas.__version__}')"
colabfold_batch --help > /dev/null && echo "  colabfold ok"
cf-random --help > /dev/null && echo "  cf-random ok"

echo ""
echo "=== Installation complete ==="
echo "Activate with: conda activate cf-random"