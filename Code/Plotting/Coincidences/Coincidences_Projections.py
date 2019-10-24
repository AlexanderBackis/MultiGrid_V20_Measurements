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

def coincidences_projections_plot(df, bus_start, bus_stop):
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
    fig.set_figheight(5)
    fig.set_figwidth(14)
    # Calculate colorbar limits
    if df.shape[0] != 0:
        vmin = 1
        vmax = df.shape[0] // 20 + 5
    else:
        vmin = 1
        vmax = 1
    # Plot
    wChs, gChs, Buses = df.wCh, df.gCh, df.Bus
    plt.subplot(1, 3, 1)
    h_front = plot_front(wChs, gChs, Buses, bus_start, bus_stop, 3e1, 5e6)
    plt.subplot(1, 3, 2)
    h_top = plot_top(wChs, gChs, Buses, bus_start, bus_stop, 2e3 ,2e6)
    plt.subplot(1, 3, 3)
    h_side = plot_side(wChs, gChs, Buses, 2e2, 6e5)
    # Collect all histograms and tighted layout
    plt.tight_layout()
    histograms = [h_front, h_top, h_side]
    return fig, histograms

# =============================================================================
#                              HELPER FUNCTIONS
# =============================================================================

def plot_front(wChs, gChs, Buses, bus_start, bus_stop, vmin=None, vmax=None):
    rows = ((bus_stop + 1) - bus_start) * 4
    h_front, *_ = plt.hist2d((wChs + (80*Buses)) // 20, gChs, bins=[rows, 40],
                             range=[[bus_start*4-0.5, (bus_stop+1)*4-0.5], [79.5, 119.5]],
                             norm=LogNorm(),
                             vmin=vmin, vmax=vmax,
                             cmap='jet')
    plt.title('Front view')
    plt.xlabel('Row')
    plt.ylabel('Grid')
    plt.colorbar()
    return h_front

def plot_top(wChs, gChs, Buses, bus_start, bus_stop, vmin=None, vmax=None):
    rows = ((bus_stop + 1) - bus_start) * 4
    h_top, *_ = plt.hist2d((wChs + (80*Buses)) // 20, wChs % 20, bins=[rows, 20],
                           range=[[bus_start*4-0.5, (bus_stop+1)*4-0.5], [-0.5, 19.5]],
                           norm=LogNorm(),
                           vmin=vmin, vmax=vmax,
                           cmap='jet')
    plt.title('Top view')
    plt.xlabel('Row')
    plt.ylabel('Layer')
    plt.colorbar()
    return h_top

def plot_side(wChs, gChs, Buses, vmin=None, vmax=None):
    h_side, *_ = plt.hist2d(wChs % 20, gChs, bins=[20, 40],
                            range=[[-0.5, 19.5], [79.5, 119.5]], norm=LogNorm(),
                            vmin=vmin, vmax=vmax,
                            cmap='jet')
    plt.title('Side view')
    plt.xlabel('Layer')
    plt.ylabel('Grid')
    plt.colorbar()
    return h_side
