#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Coincidences_Projections.py: Helper functions for handling of paths and folders.
"""

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import pandas as pd
import numpy as np

# =============================================================================
#             Coincidence Histogram Projections (Front, Top, Side)
# =============================================================================

def coincidences_projections_plot(df):
    """
    Histograms the hitposition, histogrammed over the three projections of the
    detector (Front, Top and Side).

    Args:
        df (DataFrame): Clustered events

    Returns:
        fig (Figure): Figure containing three 2D coincidences histograms, one
                      for each projection (Front, Top, Side)
        histograms (list): List containing three elements with a
                           2D numpy array 2D matrix each. Each matrix contains
                           the histogram from a specfic projection (front, top,
                           side).
    """
    # Ensure we only plot coincident events
    df = df[(df.wCh != -1) & (df.gCh != -1)]
    # Define figure and set figure properties
    fig = plt.figure()
    fig.suptitle('Coincident events', x=0.5, y=1.08)
    fig.set_figheight(4)
    fig.set_figwidth(14)
    # Calculate colorbar limits
    if df.shape[0] != 0:
        vmin = 1
        vmax = df.shape[0] // 200 + 5
    else:
        vmin = 1
        vmax = 1
    wires, grids, buses = df.wCh, df.gCh, df.Bus
    # Plot front view
    plt.subplot(1, 3, 1)
    h_front, *_ = plt.hist2d((wires + (80*buses)) // 20, grids, bins=[36, 40],
                             range=[[-0.5, 35.5], [79.5, 119.5]],
                             norm=LogNorm(),
                             cmap='jet', vmin=vmin, vmax=vmax)
    plt.title('Front view')
    plt.xlabel('Row')
    plt.ylabel('Grid')
    plt.colorbar()
    # Plot top view
    plt.subplot(1, 3, 2)
    h_top, *_ = plt.hist2d((wires + (80*buses)) // 20, wires % 20, bins=[36, 20],
                           range=[[-0.5, 35.5], [-0.5, 19.5]], norm=LogNorm(),
                           cmap='jet', vmin=vmin, vmax=vmax)
    plt.title('Top view')
    plt.xlabel('Row')
    plt.ylabel('Layer')
    plt.colorbar()
    # Plot side view
    plt.subplot(1, 3, 3)
    h_side, *_ = plt.hist2d(wires % 20, grids, bins=[20, 40],
                            range=[[-0.5, 19.5], [79.5, 119.5]], norm=LogNorm(),
                            cmap='jet', vmin=vmin, vmax=vmax)
    plt.title('Side view')
    plt.xlabel('Layer')
    plt.ylabel('Grid')
    plt.colorbar()
    # Collect all histograms and tighted layout
    plt.tight_layout()
    histograms = [h_front, h_top, h_side]
    return fig, histograms
