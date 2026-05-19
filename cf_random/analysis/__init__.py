#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Analysis module for structural quality metrics.

Provides classes for computing TM-scores and pLDDT metrics focused on
fold-switching regions and structural comparisons.
"""

from .base import (
    BaseTMScore,
)
from .cal_plddt_ac_fs import (
    PlddtCal,
)
from .cal_tmscore_fs_flmsa import (
    TMScoreFS,
)
from .cal_tmscore_fs_multimer import (
    TMScoreFSMulti,
)
from .tmscore_all_var import (
    TMScoreCalAllVar,
)
from .tmscore_all_var_fs import (
    TMScoreCalAllVarFS,
)

__all__ = [
    "BaseTMScore",
    "TMScoreFS",
    "PlddtCal",
    "TMScoreCalAllVar",
    "TMScoreCalAllVarFS",
]
