#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Duration.py: Helper functions for handling of paths and folders.
"""

import numpy as np

# =============================================================================
#                            CALCULATE DURATION
# =============================================================================

def get_duration(df):
    times = df.Time.values
    diff = np.diff(times)
    resets = np.where(diff < 0)
    duration = sum(times[resets]) + times[-1]
    return duration
