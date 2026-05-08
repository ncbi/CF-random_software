#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Prediction orchestration module.

Provides classes for running ColabFold predictions with various MSA
configurations and computing quality metrics for predicted structures.
"""

from .base import (
    ColabFoldRunner,
    MSAMaxRunner,
    MSAVariableRunner,
)
from .pred_cal_tmscore_fs import (
    MSAMaxFoldSwitch,
    MSAVariableFoldSwitch,
    TMScore,
)
from .pred_cal_tmscore_multimer import (
    MSAMaxMultimer,
    MSAVariableMultimer,
    TMScoreMonomer,
)
from .prediction_all_var import (
    PredictionAll,
)

__all__ = [
    "ColabFoldRunner",
    "MSAMaxRunner",
    "MSAVariableRunner",
    "PredictionAll",
    "TMScore",
    "TMScoreMonomer",
    "MSAMaxFoldSwitch",
    "MSAVariableFoldSwitch",
    "MSAMaxMultimer",
    "MSAVariableMultimer",
]
