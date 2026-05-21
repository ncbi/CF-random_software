#!/bin/bash
set -e

echo "=== CF-random Installation ==="

# Initialise conda for non-interactive shells
source "$(conda info --base)/etc/profile.d/conda.sh"

echo "[1/3] Setting up conda environment 'cf-random'..."
if conda env list | grep -q "^cf-random "; then
    echo "  Conda environment 'cf-random' already exists; skipping creation."
else
    echo "  Creating conda environment from environment.yml..."
    conda env create -f environment.yml
fi

conda activate cf-random

echo "[2/3] Installing pip dependencies..."
pip install --no-cache-dir --default-timeout=100 -e .

echo "[3/3] Installing JAX..."
JAX_VERSION=0.4.28
if command -v nvidia-smi &>/dev/null; then
    CUDA_VERSION=$(nvidia-smi | grep -oP "CUDA Version: \K[0-9]+" | head -1)
    echo "  GPU detected (CUDA $CUDA_VERSION), installing GPU-enabled JAX..."
    if [ "$CUDA_VERSION" -ge 12 ]; then
        pip install --retries 5 --timeout 120 "jax[cuda12_pip]==${JAX_VERSION}" \
            -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html
    else
        pip install --retries 5 --timeout 120 "jax[cuda11_pip]==${JAX_VERSION}" \
            -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html
    fi
else
    echo "  No GPU detected, installing CPU-only JAX..."
    pip install --retries 5 --timeout 120 "jax[cpu]==${JAX_VERSION}" \
        -f https://storage.googleapis.com/jax-releases/jax_releases.html
fi

echo ""
echo "=== Verifying installation ==="
python -c "import importlib.util; print('  biopython:', importlib.util.find_spec('Bio') is not None)"
python -c "import importlib.util; print('  numpy:', importlib.util.find_spec('numpy') is not None)"
python -c "import importlib.util; print('  jax:', importlib.util.find_spec('jax') is not None)"
for cmd in colabfold_batch cf-random foldseek; do
    if command -v "$cmd" &>/dev/null; then
        echo "  $cmd: ok"
    else
        echo "  $cmd: not found"
    fi
done

echo ""
echo "=== Installation complete ==="
echo "Activate with: conda activate cf-random"