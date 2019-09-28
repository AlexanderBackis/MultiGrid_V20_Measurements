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


def timestamp_plot(df, number_bins):
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
    fig = plt.figure()
    plt.title('Timestamp Histogram')
    plt.xlabel('Time [Hours]')
    plt.ylabel('Counts [events/bin]')
    plt.grid(True, which='major', zorder=0)
    plt.grid(True, which='minor', linestyle='--', zorder=0)
    # Plot
    plt.hist((df.Time*62.5e-9)/(60 ** 2), histtype='step', bins=number_bins,
             color='black', zorder=5)
    return fig
