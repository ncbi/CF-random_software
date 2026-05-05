# Troubleshooting CF-random

## Common Issues and Solutions

### 1. AttributeError: module 'numpy' has no attribute 'int' ⚠️ CRITICAL

**Error:**
```
AttributeError: module 'numpy' has no attribute 'int'.
`np.int` was a deprecated alias for the builtin `int`. 
```

**Cause:** AlphaFold/ColabFold uses deprecated NumPy syntax (`np.int`) that was removed in NumPy 2.0+. This is a **compatibility issue between AlphaFold and newer NumPy versions**.

**Solution:** Pin NumPy to a compatible version (1.26.x):

```bash
# Option 1: Install compatible numpy version
pip install 'numpy==1.26.4'

# Option 2: Force reinstall alphafold with compatible numpy
pip install --force-reinstall --no-cache-dir 'numpy==1.26.4' alphafold

# Option 3: If using conda (recommended)
conda install numpy=1.26.4
conda install alphafold
```

**Best Practice:** Create a fresh environment with tested versions:
```bash
conda create --name cf-random-fix python=3.10 numpy=1.26.4
conda activate cf-random-fix
pip install cf-random[colabfold]
```

---

### 2. ModuleNotFoundError: No module named 'tree'

**Error:**
```
ModuleNotFoundError: No module named 'tree'
```

**Cause:** This is a ColabFold dependency that's not automatically installed. The `tree` module is a Google internal library used by AlphaFold.

**Solution:** The official colabfold package should include this. If you get this error:
```bash
pip install dm-tree
```

---

### 3. FileNotFoundError: foldseek: command not found

**Error:**
```
FileNotFoundError: [Errno 2] No such file or directory: 'foldseek'
```

**Cause:** Foldseek is not installed or not in your PATH.

**Solution:** 
```bash
# Install foldseek via conda
conda install -c conda-forge -c bioconda foldseek

# Verify it's installed
which foldseek

# Set up databases (run once)
foldseek databases PDB pdb tmp
```

---

### 4. ImportError: No module named 'pymol'

**Error:**
```
ImportError: No module named 'pymol'
```

**Cause:** PyMOL is an optional dependency that must be installed via conda.

**Solution:**
```bash
conda install conda-forge::pymol-open-source
```

Or if you don't need PyMOL visualization:
```bash
pip install cf-random
# The package will gracefully handle missing pymol
```

---

### 5. ModuleNotFoundError during `cf-random` CLI

**Error:**
```
ModuleNotFoundError: No module named 'alphafold'
```

**Cause:** You're not using the localcolabfold conda environment which includes AlphaFold.

**Solution:**
```bash
# Activate the correct environment
conda activate cf-random

# Then run cf-random
cf-random --help
```

---

### 6. colabfold_batch: command not found

**Error:**
```
FileNotFoundError: [Errno 2] No such file or directory: 'colabfold_batch'
```

**Cause:** The ColabFold package is not installed or your conda environment is not activated.

**Solution:**
```bash
# Activate your CF-random conda environment
conda activate cf-random

# Verify colabfold is installed
which colabfold_batch
pip show colabfold

# If not installed, install it
pip install colabfold
```

---

### 7. colabfold_batch fails silently (no output directories created)

**Cause:** ColabFold is failing but the error is being swallowed. This often happens due to:
- NumPy version incompatibility (see Issue #1)
- Missing dependencies
- GPU/CUDA issues

**Solution:**
1. Check the NumPy version first: See Issue #1 above
2. Verify all dependencies:
   ```bash
   conda list | grep -E "colabfold|alphafold|foldseek|numpy"
   ```
3. Try running colabfold directly to see detailed errors:
   ```bash
   colabfold_batch --help
   ```

---

### 8. Permission errors when accessing MSA files

**Error:**
```
PermissionError: [Errno 13] Permission denied: '...'
```

**Cause:** MSA folder doesn't have read permissions or path is incorrect.

**Solution:**
```bash
# Check permissions
ls -la /path/to/msa/folder

# Fix permissions if needed
chmod -R 755 /path/to/msa/folder

# Ensure path is correct when calling cf-random
cf-random --fname /full/path/to/msa/folder --pdb1 structure.pdb --option AC
```

---

### 9. mv: rename failed: No such file or directory

**Error:**
```
mv: rename test_predicted_models_full_rand_8/ to blind_prediction/test/test_predicted_models_full_rand_8/: No such file or directory
```

**Cause:** ColabFold prediction failed, so the output directory was never created, and the move command failed.

**Solution:**
1. Look at the errors printed **above** this mv error - they will show why colabfold_batch failed
2. Fix the underlying ColabFold issue (usually NumPy compatibility - see Issue #1)
3. Retry the prediction:
   ```bash
   cf-random --pname test --fname test_data/ --option blind
   ```

---

## Environment Setup Verification

Run this script to verify all dependencies:

```bash
#!/bin/bash

echo "=== CF-random Environment Check ==="
echo ""

# Check conda environment
echo "1. Conda environment:"
conda info --envs | grep cf-random

# Check Python packages
echo ""
echo "2. CF-random package:"
pip show cf-random

# Check NumPy version (critical!)
echo ""
echo "3. NumPy version (must be < 2.0):"
python -c "import numpy; print(f'NumPy: {numpy.__version__}')"

# Check colabfold
echo ""
echo "4. ColabFold:"
which colabfold_batch

# Check foldseek
echo ""
echo "5. Foldseek:"
which foldseek

# Check optional packages
echo ""
echo "6. Optional packages:"
python -c "
try:
    import pymol
    print('✓ PyMOL: installed')
except ImportError:
    print('✗ PyMOL: not installed (optional)')

try:
    import alphafold
    print('✓ AlphaFold: installed')
except ImportError:
    print('✗ AlphaFold: not installed')

try:
    import colabfold
    print('✓ ColabFold: installed')
except ImportError:
    print('✗ ColabFold: not installed')
"

echo ""
echo "=== Check complete ==="
```

---

## Recommended Quick Fix

If you're experiencing issues, try this clean environment setup using the **official ColabFold package**:

```bash
# Create fresh environment with compatible NumPy
conda create -n cf-random-clean python=3.10 numpy=1.26.4 -y
conda activate cf-random-clean

# Install CF-random with ColabFold support
pip install -e ".[colabfold]"

# Install foldseek for blind mode
conda install -c conda-forge -c bioconda foldseek

# Test
cf-random --help
colabfold_batch --help
```

---

## Getting Help

If you encounter issues not listed here:

1. Check the examples in `examples/` directory
2. Review the main [README.md](README.md) for installation details
3. Check [INSTALL.md](INSTALL.md) for step-by-step setup
4. For ColabFold issues, see: https://github.com/sokrypton/ColabFold
5. For Foldseek issues, see: https://github.com/steineggerlab/foldseek
6. For AlphaFold NumPy issues, check: https://github.com/deepmind/alphafold/issues
