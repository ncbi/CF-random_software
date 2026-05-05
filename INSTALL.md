# Installation Guide for CF-random

## Requirements

- Python >= 3.10
- Conda (for managing environment)
- pip (for Python packages)

## System Requirements

CF-random is tested on Linux. For macOS and Windows, use the [Colab notebook](https://colab.research.google.com/drive/16pD2tUMkUx1gwDxZXcSr9WOosYp0ZU6j?authuser=0).

## Installation Steps (Recommended - Using Official ColabFold)

### 1. Create a conda environment

```bash
conda create --name cf-random python=3.10 -y
conda activate cf-random
```

### 2. Install CF-random with ColabFold

The easiest way to get started is with the full installation:

```bash
pip install -e ".[colabfold]"
```

This installs CF-random and all ColabFold dependencies including:
- colabfold (official package)
- alphafold
- dm-tree
- All other required Python packages

### 3. Install foldseek and optional visualization tools

```bash
# Install foldseek
conda install -c conda-forge -c bioconda foldseek

# Optional: Install PyMOL for visualization
conda install conda-forge::pymol-open-source

# Set up foldseek databases (run once)
foldseek databases PDB pdb tmp
```

**Note:** The foldseek databases command should be run in a directory where you want to store the database files.

## Troubleshooting Installation

### If you encounter NumPy compatibility issues:

The official colabfold package requires NumPy < 2.0. If you already have a conflicting NumPy version:

```bash
pip install 'numpy<2.0' --force-reinstall --no-cache-dir
```

### Alternative Installation Methods

#### Option A: Minimal Installation (core only)

If you want to install CF-random without ColabFold:

```bash
pip install -e .
```

#### Option B: Development Installation

For developers contributing to CF-random:

```bash
pip install -e ".[full]"
```

This includes ColabFold + dev tools (pytest, black, isort, flake8).

#### Option C: From PyPI (when available)

Once published:

```bash
pip install cf-random
```

## Verify Installation

Run these commands to verify everything is set up correctly:

```bash
# Verify CF-random package
pip show cf-random

# Test the CLI
cf-random --help

# Test Python import
python -c "import cf_random; print(f'CF-random version: {cf_random.__version__}')"

# Verify ColabFold is working
colabfold_batch --help

# Verify foldseek is installed
which foldseek
```

If all commands succeed, you're ready to use CF-random!

## Usage

After installation, you can use CF-random via the command line:

```bash
cf-random --help
```

Or import it in Python:
```python
import cf_random
```

### Examples

See the `examples/` directory for detailed examples:
- `1_FS_mode-RfaH/` - Fold-switching prediction example
- `2_Beta-phosphoglucomutase/` - Alternative conformation example
- `3_blind_mode-Mad2/` - Blind mode prediction example

## Troubleshooting

### Import errors
Make sure you've activated the conda environment:
```bash
conda activate cf-random
```

### Missing dependencies
If you get import errors, reinstall the package:
```bash
pip install -e .
```

### Foldseek not found
Ensure foldseek is installed and in your PATH:
```bash
which foldseek
```

## Documentation

For more information, see the main [README.md](README.md) or visit the [project repository](https://github.com/your-repo/cf-random).
