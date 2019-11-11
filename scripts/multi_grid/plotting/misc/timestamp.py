#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Timestamp.py: Visualises the distrubtion of timestamps versus cluster index.
              Uses a color scale to show summed multiplicity of the clusters.
"""

import matplotlib.pyplot as plt
import numpy as np

# =============================================================================
#                                   Timestamp
# =============================================================================


def timestamp_plot(Time, number_bins, unit, label):
    """
    Scatter plot of cluster index vs timestamp, where every 100:th cluster is
    plotted. The color bar shows the summation of wire and grid event, which
    gives an indication of the 'healthiness' of clusters.

    Args:
        df (DataFrame): Clustered events

    Returns:
        fig (Figure): Figure containing the scatter plot of the timestamps
                      from all of the detectors.
    """


    # Prepare figure
    plt.title('Timestamp Histogram')
    plt.xlabel('Time (%s)' % unit)
    plt.ylabel('Rate (events/s)')
    plt.grid(True, which='major', zorder=0)
    plt.grid(True, which='minor', linestyle='--', zorder=0)
    # Plot
    #plt.hist(Time, histtype='step', bins=number_bins, label=label, zorder=5)
    hist, bin_edges = np.histogram(Time, bins=number_bins)
    bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
    delta_t = (60 ** 2) * (bin_centers[1] - bin_centers[0])
    plt.plot(bin_centers, hist/delta_t, label=label, zorder=5)
