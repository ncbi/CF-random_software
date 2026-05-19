#!/bin/bash
set -e
echo "=== CF-random Installation ==="
echo "[1/2] Ensuring conda environment 'cf-random' exists..."
eval "$(conda shell.bash hook)"
if conda env list | grep -q "cf-random"; then
    echo "  Conda environment 'cf-random' already exists; skipping creation."
else
    echo "  Creating conda environment from environment.yml..."
    conda env create -f environment.yml -n cf-random
fi
conda activate cf-random

# Install JAX with GPU or CPU depending on hardware
JAX_VERSION=0.4.28
if command -v nvidia-smi &> /dev/null; then
    CUDA_VERSION=$(nvidia-smi | grep -oP "CUDA Version: \K[0-9]+" | head -1)
    echo "  GPU detected (CUDA $CUDA_VERSION), installing GPU-enabled JAX..."
    if [ "$CUDA_VERSION" -ge 12 ]; then
        pip install "jax[cuda12_pip]==${JAX_VERSION}" \
            -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html
    else
        pip install "jax[cuda11_pip]==${JAX_VERSION}" \
            -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html
    fi
else
    echo "  No GPU detected, installing CPU-only JAX..."
    pip install "jax[cpu]==${JAX_VERSION}" \
        -f https://storage.googleapis.com/jax-releases/jax_releases.html
fi

echo "[2/2] Verifying installation (best-effort checks)..."
python -c "import importlib.util; print('  biopython:', importlib.util.find_spec('Bio') is not None)"
python -c "import importlib.util; print('  numpy:', importlib.util.find_spec('numpy') is not None)"
python -c "import importlib.util; print('  jax:', importlib.util.find_spec('jax') is not None)"
if command -v colabfold_batch &> /dev/null; then
    echo "  colabfold ok"
fi
if command -v cf-random &> /dev/null; then
    echo "  cf-random ok"
fi

echo ""
echo "=== Installation complete ==="
echo "Activate with: conda activate cf-random"