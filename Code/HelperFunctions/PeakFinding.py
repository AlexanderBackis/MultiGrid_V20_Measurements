#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PeakFinding.py: Helper functions for handling of paths and folders.
"""
import numpy as np
from scipy.signal import find_peaks

# =============================================================================
#                          FIND PEAKS IN HISTOGRAM
# =============================================================================

def get_peaks(hist, heights, number_bins):
    # Extract heights
    height_1 = heights[0]
    height_2 = heights[1]
    # Define vectors with heights, used to find peaks above height values
    bins_part_1 = number_bins - number_bins//3
    bins_part_2 = number_bins//3
    heights_part1 = np.ones(bins_part_1)*height_1
    heights_part2 = np.ones(bins_part_2)*height_2
    heights = np.append(heights_part1, heights_part2)
    # Get peaks
    peaks, *_ = find_peaks(hist, height=heights)
    return peaks
