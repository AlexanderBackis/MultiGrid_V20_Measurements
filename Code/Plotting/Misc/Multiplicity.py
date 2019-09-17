#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Multiplicity.py: Histograms multiplicity of wires versus grids in the
                 clustered neutron events.
"""

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# =============================================================================
#                               Multiplicity
# =============================================================================

def multiplicity_plot(df):
    """
    Histograms multiplicity of wires versus grids in the clustered neutron
    events.

    Args:
        df (DataFrame): Clustered events

    Returns:
        fig (Figure): Figure containing nine 2D coincidences histograms, one
                      for each bus.

    """
    def plot_multiplicity_bus(df, bus, vmin, vmax):
        # Plot data
        plt.hist2d(df.wM, df.gM, bins=[80, 40], range=[[0, 80], [0, 40]],
                   norm=LogNorm(), vmin=vmin, vmax=vmax, cmap='jet')
        plt.xlabel('Wire Multiplicity')
        plt.ylabel('Grid Multiplicity')
        plt.colorbar()
        plt.tight_layout()
        plt.title('Bus %d\n(%d events)' % (bus, df.shape[0]))
    # Set limits
    if df.shape[0] != 0:
        vmin = 1
        vmax = df.shape[0] // 9 + 10
    else:
        vmin = 1
        vmax = 1
    # Prepare figure
    fig = plt.figure()
    fig.set_figheight(14)
    fig.set_figwidth(12)
    # Iterate through all buses
    for bus in range(0, 9):
        plt.subplot(3, 3, bus+1)
        df_bus = df[df.Bus == bus]
        plot_multiplicity_bus(df_bus, bus, vmin, vmax)
    plt.tight_layout()
    return fig
