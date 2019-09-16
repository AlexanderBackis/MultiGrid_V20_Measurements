#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Coincidences_2D.py: Helper functions for handling of paths and folders.
"""

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# =============================================================================
#                          Coincidence Histogram (2D)
# =============================================================================

def coincidences_2D_plot(ce, measurement_time):
    def plot_2D_bus(fig, sub_title, ce, vmin, vmax, duration):
        h, *_ = plt.hist2d(ce.wCh, ce.gCh, bins=[80, 40],
                           range=[[-0.5, 79.5], [79.5, 119.5]],
                           vmin=vmin, vmax=vmax, norm=LogNorm(), cmap='jet')
        plt.xlabel('Wire [Channel number]')
        plt.ylabel('Grid [Channel number]')
        plt.title(sub_title)
        plt.colorbar()
        return fig, h

    # Calculate color limits
    if ce.shape[0] != 0:
        duration = measurement_time
        vmin = 1
        vmax = ce.shape[0] // 4500 + 5
    else:
        duration = 1
        vmin = 1
        vmax = 1
    # Plot data
    fig = plt.figure()
    fig.set_figwidth(14)
    fig.set_figheight(12)
    histograms = []
    for bus in range(0, 9):
        ce_bus = ce[ce.Bus == bus]
        number_events = ce_bus.shape[0]
        events_per_s = round(number_events/duration, 4)
        sub_title = ('Bus %d\n(%d events, %f events/s)' % (bus, number_events, events_per_s))
        plt.subplot(3, 3, bus+1)
        fig, h = plot_2D_bus(fig, sub_title, ce_bus, vmin, vmax, duration)
        histograms.append(h)
    plt.tight_layout()
    return fig, histograms
