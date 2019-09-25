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


def timestamp_plot(df):
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

    # Filter so we only plot every 100:th element, so it is not so slow
    df_temp = df #[df.index % 100 == 0]
    # Prepare figure
    fig = plt.figure()
    plt.title('Timestamp (every 100:th cluster)')
    plt.xlabel('Cluster [index]')
    plt.ylabel('Timestamp [TDC channels]')
    plt.grid(True, which='major', zorder=0)
    plt.grid(True, which='minor', linestyle='--', zorder=0)
    # Plot
    plt.scatter(df_temp.index, df_temp.Time, c=df_temp.wM+df_temp.gM,
                cmap='jet', zorder=5)
    cbar = plt.colorbar()
    cbar.set_label('Multiplicity summation [wM + gM]')
    plt.tight_layout()
    return fig
