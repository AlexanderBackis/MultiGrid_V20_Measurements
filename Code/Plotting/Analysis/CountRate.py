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
    Calculates the count rate. Does this in five steps:

    1. Locates peak index (IMPORANT: The prompt peak must be filtered out before
       the procedure starts, if this is more intense than elastic peak)

    2. Finds edges of the elastic peak corresponding to the FWHM

    3. Calculate duration of peak for all periods

    4. Calculate number of counts within peak for all periods

    5. Calculate count rate according to counts/time

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
    bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
    # Calculate range in ToF equivalent to FWHM of the peak
    peak_idxs = np.where(ToF_hist == max(ToF_hist))
    peak_idx = peak_idxs[len(peak_idxs)//2][0]
    start_idx = find_nearest(ToF_hist[:peak_idx], ToF_hist[peak_idx]/2)
    stop_idx = find_nearest(ToF_hist[peak_idx:], ToF_hist[peak_idx]/2)
    start, stop = bin_centers[start_idx], bin_centers[peak_idx+stop_idx]
    # Calculate counts in peak
    ToF_values_peak = ToF_values[(ToF_values >= start) & (ToF_values <= stop)]
    counts_peak = len(ToF_values_peak)
    # Calculate peak duration
    number_of_periods = measurement_time/PERIOD_TIME
    duration_peak_per_period = stop - start
    duration_peak = number_of_periods * duration_peak_per_period
    # Calculate rate
    rate = counts_peak/duration_peak
    # Print statements for debugging purposes
    print()
    print('**** COUNT RATE ****')
    print('------------------')
    print('Start: %f' % start)
    print('Stop: %f' % stop)
    print('Counts peak: %f' % counts_peak)
    print('Number of periods: %f' % number_of_periods)
    print('Duration of peak per period: %f ' % duration_peak_per_period)
    print('Full duration of peak: %f' % duration_peak)
    print('------------------')
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
