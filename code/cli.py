#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Command-line interface for CF-random

Wrapper around main.py functionality
"""

import sys
import os
from pathlib import Path

# Add the package directory to path for imports
PACKAGE_DIR = Path(__file__).parent.absolute()


def main():
    """Main entry point for cf-random CLI"""
    sys.path.insert(0, str(PACKAGE_DIR))
    from main import *

    # Call the main function from main.py
    # The main.py file uses argparse, so it will handle command-line arguments automatically
    if __name__ == "__main__":
        pass


if __name__ == "__main__":
    main()
