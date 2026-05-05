#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CF-random: Main module for predicting alternative conformations and fold-switching proteins

This module serves as the entry point for CF-random analysis.
"""

import sys
from pathlib import Path

# Add parent code directory to path for compatibility
_CODE_DIR = Path(__file__).parent.parent / "code"
if str(_CODE_DIR) not in sys.path:
    sys.path.insert(0, str(_CODE_DIR))

# Re-export main functionality from code directory
from main import *  # noqa: F401, F403
