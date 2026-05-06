#!/bin/bash
set -e

echo "=== CF-random Installation ==="
echo "[1/2] Creating conda environment..."
conda env create -f environment.yml
eval "$(conda shell.bash hook)"
conda activate cf-random

# Install JAX with GPU or CPU depending on hardware
if command -v nvidia-smi &> /dev/null; then
    CUDA_VERSION=$(nvidia-smi | grep -oP "CUDA Version: \K[0-9]+" | head -1)
    echo "  GPU detected (CUDA $CUDA_VERSION), installing GPU-enabled JAX..."
    if [ "$CUDA_VERSION" -ge 12 ]; then
        pip install "jax[cuda12_pip]==0.4.25" \
            -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html
    else
        pip install "jax[cuda11_pip]==0.4.25" \
            -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html
    fi
else
    echo "  No GPU detected, installing CPU-only JAX..."
    pip install "jax==0.4.25" "jaxlib==0.4.25"
fi

echo "[2/2] Verifying installation..."
python -c "from Bio.Data import SCOPData; print('  biopython ok')"
python -c "import numpy; print(f'  numpy {numpy.__version__}')"
python -c "import pandas; print(f'  pandas {pandas.__version__}')"
python -c "import jax; print(f'  jax {jax.__version__} | devices: {jax.devices()}')"
colabfold_batch --help > /dev/null && echo "  colabfold ok"
cf-random --help > /dev/null && echo "  cf-random ok"

echo ""
echo "=== Installation complete ==="
echo "Activate with: conda activate cf-random"