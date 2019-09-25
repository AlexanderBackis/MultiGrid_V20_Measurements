#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ToF.py: Histograms the ToF values in the clustered data.
"""

import matplotlib.pyplot as plt

# =============================================================================
#                               ToF - MG
# =============================================================================


def ToF_histogram(df, number_bins):
    """
    Histograms the ToF values in the clustered data.

    Args:
        df (DataFrame): Clustered events
        number_bins (int): The number of bins to histogram the data into

    Returns:
        fig (Figure): Figure containing the ToF histogram
    """
    # Prepare figure
    fig = plt.figure()
    plt.title('ToF')
    plt.xlabel('ToF [$\mu$s]')
    plt.ylabel('Counts')
    plt.grid(True, which='major', linestyle='--', zorder=0)
    plt.grid(True, which='minor', linestyle='--', zorder=0)
    # Histogram data
    plt.hist(df.ToF * 62.5e-9 * 1e6, bins=number_bins, color='black',
             zorder=4, histtype='step')
    return fig
