#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CF-random: Command-line interface

Main entry point for running CF-random from the command line.
"""

import sys
import warnings
from pathlib import Path

# Add the code directory to path so imports work
_CODE_DIR = Path(__file__).parent.parent / "code"
if str(_CODE_DIR) not in sys.path:
    sys.path.insert(0, str(_CODE_DIR))

warnings.filterwarnings('ignore')


def main():
    """Main entry point for the cf-random CLI
    
    This function is called by the entry point and delegates to the
    main.py script which handles all argument parsing and execution.
    """
    # Dynamically execute the main.py module
    import runpy
    main_script = _CODE_DIR / "main.py"
    runpy.run_path(str(main_script), run_name="__main__")


if __name__ == "__main__":
    main()

