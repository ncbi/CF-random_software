# Installation Guide for CF-random

## Requirements

- Python >= 3.10
- Conda (for installing foldseek and pymol)
- localcolabfold (see system requirements below)

## System Requirements

CF-random currently supports Linux environments. For macOS and Windows, use the [Colab notebook](https://colab.research.google.com/drive/16pD2tUMkUx1gwDxZXcSr9WOosYp0ZU6j?authuser=0).

## Installation Steps

### 1. Set up localcolabfold (Linux)

```bash
wget https://raw.githubusercontent.com/YoshitakaMo/localcolabfold/main/install_colabbatch_linux.sh
bash install_colabbatch_linux.sh
```

Then add to your `.bashrc`:
```bash
export PATH="/path/to/your/localcolabfold/colabfold-conda/bin:$PATH"
```

Reactivate bash: `source ~/.bashrc`

### 2. Install foldseek and other dependencies

Create and activate a conda environment:
```bash
conda create --name cf-random python=3.10
conda activate cf-random
```

Install foldseek and pymol:
```bash
conda install -c conda-forge -c bioconda foldseek
conda install conda-forge::pymol-open-source
```

Set up foldseek databases:
```bash
foldseek databases PDB pdb tmp
```

### 3. Install CF-random as a package

#### Option A: From source (development mode)
```bash
git clone https://github.com/your-repo/cf-random.git
cd cf-random
pip install -e .
```

#### Option B: From PyPI (when available)
```bash
pip install cf-random
```

#### Option C: Standard installation
```bash
cd cf-random
pip install .
```

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
