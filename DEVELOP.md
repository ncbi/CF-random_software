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
conda create --name cf-random-dev python=3.10
conda activate cf-random-dev
```

3. Install in development mode with dev dependencies:
```bash
pip install -e ".[dev]"
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
