#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CF-random: Command-line interface

Main entry point for running CF-random from the command line.
"""

import sys
from pathlib import Path

# Add the code directory to path so imports work
_CODE_DIR = Path(__file__).parent.parent / "code"
if str(_CODE_DIR) not in sys.path:
    sys.path.insert(0, str(_CODE_DIR))


def main():
    """Main entry point for the cf-random CLI"""
    # Import and run the main module which handles all CLI arguments
    import warnings
    warnings.filterwarnings('ignore')
    
    from pred_cal_tmscore_FS import *
    from pred_cal_tmscore_blind import *
    from pred_cal_tmscore_AC import *
    from cal_plddt_ACFS import *
    from PLOT_AC import *
    from PLOT_FS import *
    from search_w_foldseek_cluster import *
    
    # The imports above will register the argparse handlers
    # If you want to add additional CLI logic, do so here


if __name__ == "__main__":
    main()
