#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Setup script for CF-random package
"""

from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="cf-random",
    version="0.1.0",
    author="Myeongsang (Samuel) Lee",
    description="Predicting alternative conformations and fold-switching proteins",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/your-repo/cf-random",
    packages=find_packages(),
    python_requires=">=3.10",
    install_requires=[
        "biopython>=1.79",
        "numpy",
        "matplotlib",
        "seaborn",
        "scikit-learn",
        "mdtraj",
        "MDAnalysis",
        "textalloc",
        "tmtools",
        "adjustText",
        "thefuzz",
    ],
    entry_points={
        "console_scripts": [
            "cf-random=cf_random.cli:main",
        ],
    },
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: Public Domain",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
