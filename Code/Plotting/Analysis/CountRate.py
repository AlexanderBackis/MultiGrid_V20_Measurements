#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CountRate.py: Calculates the count rate
"""

import numpy as np

# =============================================================================
#                               COUNT RATE
# =============================================================================

def calculate_count_rate(ToF_values, measurement_time):
    """
    Calculates the count rate.

    Args:
        ToF_values (numpy array): Numpy array containing all of the ToF values
                                  corresponding to the clustered events
        measurement_time (float): Measurement time in seconds

    Returns:
        rate (float): Count rate value expressed in Hz
    """
    # Declare constants, period time is in [µs]
    PERIOD_TIME = (1/60)
    # Histogram data
    ToF_hist, bin_edges = np.histogram(ToF_values, bins=1000, range=[0, PERIOD_TIME])
    # Calculate range in ToF equivalent to FWHM of the peak
    peak_idxs = np.where(ToF_hist == max(ToF_hist))
    peak_idx = peak_idxs[len(peak_idxs)//2][0]
    start = find_nearest(hist[:peak_idx], hist[peak_idx]/2)
    stop = find_nearest(hist[peak_idx:], hist[peak_idx]/2)
    # Calculate counts in peak
    ToF_values_peak = ToF_values[(ToF_values >= start) & (ToF_values <= stop)]
    counts_peak = len(ToF_values_peak)
    # Calculate peak duration
    number_of_periods = measurement_time/PERIOD_TIME
    duration_peak_per_period = stop - start
    duration_peak = number_of_periods * duration_peak_per_period
    # Calculate rate
    rate = counts_peak/duration_peak
    return rate


# =============================================================================
#                               HELPER FUNCTIONS
# =============================================================================

def find_nearest(array, value):
    """
    Returns the index of the element in 'array' which is closest to 'value'.

    Args:
        array (numpy array): Numpy array with elements
        value (float): Value which we want to find the closest element to in
                       arrray

    Returns:
        idx (int): index of the element in 'array' which is closest to 'value'
    """
    idx = (np.abs(array - value)).argmin()
    return idx
