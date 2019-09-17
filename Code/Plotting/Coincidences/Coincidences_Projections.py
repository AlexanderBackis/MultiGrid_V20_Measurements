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
    # Plot front view
    plt.subplot(1, 3, 1)
    h_front = plot_2D_Front(df, vmin, vmax)
    # Plot top view
    plt.subplot(1, 3, 2)
    h_top = plot_2D_Top(df, vmin, vmax)
    # Plot side view
    plt.subplot(1, 3, 3)
    h_side = plot_2D_Side(df, vmin, vmax)
    plt.tight_layout()
    # Collect all histograms
    histograms = [h_front, h_top, h_side]
    return fig, histograms


# =============================================================================
#                       Coincidence Histogram - Front
# =============================================================================


def plot_2D_Front(df, vmin, vmax):
    df_tot = pd.DataFrame()
    for i in range(0, 9):
        df_clu = df[df.Bus == i]
        df_clu['wCh'] += (80 * i) + (i // 3) * 80
        df_clu['gCh'] += (-80 + 1)
        df_tot = pd.concat([df_tot, df_clu])
    h, *_ = plt.hist2d(np.floor(df_tot['wCh'] / 20).astype(int) + 1,
                       df_tot.gCh,
                       bins=[12*3 + 8, 40],
                       range=[[0.5, 12*3 + 0.5 + 8],
                              [0.5, 40.5]],
                       norm=LogNorm(), cmap='jet', vmin=vmin, vmax=vmax)
    plt.title('Front view')
    locs_x = [1, 12, 17, 28, 33, 44]
    ticks_x = [1, 12, 13, 25, 26, 38]
    plt.xlabel('Layer')
    plt.ylabel('Grid')
    plt.xticks(locs_x, ticks_x)
    plt.colorbar()
    return h

# =============================================================================
#                       Coincidence Histogram - Top
# =============================================================================


def plot_2D_Top(df, vmin, vmax):
    df_tot = pd.DataFrame()
    for i in range(0, 9):
        df_clu = df[df.Bus == i]
        df_clu['wCh'] += (80 * i) + (i // 3) * 80
        df_tot = pd.concat([df_tot, df_clu])
    h, *_ = plt.hist2d(np.floor(df_tot['wCh'] / 20).astype(int) + 1,
                       df_tot['wCh'] % 20 + 1,
                       bins=[12*3 + 8, 20],
                       range=[[0.5, 12*3 + 0.5 + 8],
                              [0.5, 20.5]],
                       norm=LogNorm(), cmap='jet', vmin=vmin, vmax=vmax)
    plt.title('Top view')
    plt.xlabel('Layer')
    plt.ylabel('Wire')
    locs_x = [1, 12, 17, 28, 33, 44]
    ticks_x = [1, 12, 13, 25, 26, 38]
    plt.xticks(locs_x, ticks_x)
    plt.colorbar()
    return h


# =============================================================================
#                       Coincidence Histogram - Side
# =============================================================================


def plot_2D_Side(df, vmin, vmax):
    df_tot = pd.DataFrame()
    for i in range(0, 9):
        df_clu = df[df.Bus == i]
        df_clu['gCh'] += (-80 + 1)
        df_tot = pd.concat([df_tot, df_clu])
    h, *_ = plt.hist2d(df_tot['wCh'] % 20 + 1, df_tot['gCh'],
                       bins=[20, 40],
                       range=[[0.5, 20.5], [0.5, 40.5]],
                       norm=LogNorm(),
                       cmap='jet', vmin=vmin, vmax=vmax)
    plt.title('Side view')
    plt.xlabel('Wire')
    plt.ylabel('Grid')
    plt.colorbar()
    return h
