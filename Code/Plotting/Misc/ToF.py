#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ToF.py: Histograms the ToF values in the clustered data.
"""

import matplotlib.pyplot as plt

# =============================================================================
#                               ToF - MG
# =============================================================================


def ToF_histogram(df, number_bins, label=None):
    """
    Histograms the ToF values in the clustered data.

    Args:
        df (DataFrame): Clustered events
        number_bins (int): The number of bins to histogram the data into

    Yields:
        Figure containing the ToF histogram, this can then be used in overlay
        or in sub plots.
    """
    # Declare parameters
    time_offset = (0.6e-3) * 1e6
    period_time = (1/14) * 1e6
    # Prepare figure
    plt.title('ToF')
    plt.xlabel('ToF [$\mu$s]')
    plt.ylabel('Counts')
    plt.yscale('log')
    plt.grid(True, which='major', linestyle='--', zorder=0)
    plt.grid(True, which='minor', linestyle='--', zorder=0)
    # Histogram data
    plt.hist((df.ToF * 62.5e-9 * 1e6 + time_offset) % period_time,
             bins=number_bins, zorder=4, histtype='step', label=label)
