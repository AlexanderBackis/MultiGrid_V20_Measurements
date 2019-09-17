#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PHS_Wires_Vs_Grids.py: Histograms ADC charge from wires vs grids, one for each
                       bus, showing the relationship between charge collected by
                       wires and charge collected by grids. In the ideal case
                       there should be linear relationship between these
                       two quantities.
"""

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# =============================================================================
#                           PHS (Wires vs Grids)
# =============================================================================


def PHS_wires_vs_grids_plot(ce):
    """
    Histograms ADC charge from wires vs grids, one for each bus, showing the
    relationship between charge collected by wires and charge collected by
    grids. In the ideal case there should be linear relationship between these
    two quantities.

    Args:
        ce (DataFrame): Clustered events

    Returns:
        fig (Figure): Figure containing nine 2D PHS plots, showing wire charge
                      versus grid charge, one plot for each bus.
    """
    def charge_scatter(fig, ce, sub_title, bus, vmin, vmax):
        plt.xlabel('Collected charge wires [ADC channels]')
        plt.ylabel('Collected charge grids [ADC channels]')
        plt.title(sub_title)
        bins = [200, 200]
        ADC_range = [[0, 5000], [0, 5000]]
        plt.hist2d(ce.wADC, ce.gADC, bins=bins, norm=LogNorm(), range=ADC_range,
                   vmin=vmin, vmax=vmax, cmap='jet')
        plt.colorbar()
        return fig
    
    # Plot data
    fig = plt.figure()
    fig.set_figheight(12)
    fig.set_figwidth(14)
    # Set color limits
    if ce.shape[0] != 0:
        vmin = 1
        vmax = ce.shape[0] // 4500 + 1000
    else:
        vmin = 1
        vmax = 1
    # Plot
    for bus in range(0, 9):
        events_bus = ce[ce.Bus == bus]
        sub_title = 'Bus %d\n(%d events)' % (bus, events_bus.shape[0])
        plt.subplot(3, 3, bus+1)
        fig = charge_scatter(fig, events_bus, sub_title, bus, vmin, vmax)
    plt.tight_layout()
    return fig
