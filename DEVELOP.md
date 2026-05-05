# Development

This document describes how to set up a development environment for CF-random.

## Setup

1. Clone the repository:
```bash
git clone https://github.com/your-repo/cf-random.git
cd cf-random
```

2. Create a conda environment:
```bash
conda create --name cf-random-dev python=3.10 numpy=1.26.4
conda activate cf-random-dev
```

3. Install in development mode with all dependencies:
```bash
# Option 1: Development + ColabFold + full suite
pip install -e ".[full]"

# Option 2: Development only (minimal dependencies)
pip install -e ".[dev]"

# Option 3: ColabFold only
pip install -e ".[colabfold]"
```

## Running Tests

```bash
pytest tests/
```

## Code Style

We use Black and isort for code formatting:

```bash
# Format code
black cf_random/
isort cf_random/

# Check style
flake8 cf_random/
```

## Building the Package

To build the distribution packages:

```bash
pip install build
python -m build
```

This will create:
- `dist/cf-random-*.tar.gz` - Source distribution
- `dist/cf_random-*.whl` - Wheel distribution

## Publishing to PyPI

1. Build the package (see above)
2. Install twine: `pip install twine`
3. Upload to TestPyPI first: `twine upload -r testpypi dist/*`
4. Test installation: `pip install -i https://test.pypi.org/simple/ cf-random`
5. Upload to PyPI: `twine upload dist/*`

## Project Structure

```
cf-random/
├── cf_random/           # Main package
│   ├── __init__.py
│   ├── main.py         # Wrapper around code/main.py
│   └── cli.py          # CLI entry point
├── code/               # Original analysis code
├── Data/               # Sample datasets
├── examples/           # Usage examples
├── tests/              # Unit tests
├── pyproject.toml      # Project metadata
├── setup.py            # Setup script
├── MANIFEST.in         # Package manifest
├── INSTALL.md          # Installation guide
├── TROUBLESHOOTING.md  # Common issues and fixes
├── DEVELOP.md          # This file
└── README.md           # Project readme
```

## Modifying the Package

- Add new modules to the `cf_random/` directory
- Update `cf_random/__init__.py` to export new functionality
- Add tests to the `tests/` directory
- Update documentation in README.md

## Contributing

Please follow the code style guidelines and test your changes before submitting pull requests.

## Known Issues

### NumPy Compatibility
AlphaFold uses deprecated `np.int` syntax that was removed in NumPy 2.0. Always use NumPy < 2.0:
```bash
pip install 'numpy<2.0'
```

### ColabFold on macOS
ColabFold is primarily tested on Linux. For macOS/Windows, use the Colab notebook:
https://colab.research.google.com/drive/16pD2tUMkUx1gwDxZXcSr9WOosYp0ZU6j?authuser=0
