#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PHS_2D.py: Histograms the ADC-values from each channel individually and
           summarises it in a 2D histogram plot, where the color scale indicates
           number of counts. Each bus is presented in an individual plot.
"""

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# =============================================================================
#                                   PHS (2D)
# =============================================================================


def PHS_2D_plot(events):
    """
    Histograms the ADC-values from each channel individually and summarises it
    in a 2D histogram plot, where the color scale indicates number of counts.
    Each bus is presented in an individual plot.

    Args:
        events (DataFrame): Individual events

    Returns:
        fig (Figure): Figure containing 2D PHS plot
    """
    def PHS_2D_plot_bus(fig, events, sub_title, vmin, vmax):
        plt.xlabel('Channel')
        plt.ylabel('Charge [ADC channels]')
        plt.title(sub_title)
        bins = [120, 120]
        plt.hist2d(events.Ch, events.ADC, bins=bins, norm=LogNorm(),
                   range=[[-0.5, 119.5], [0, 4400]], vmin=vmin, vmax=vmax,
                   cmap='jet')
        plt.colorbar()
        return fig

    # Prepare figure
    fig = plt.figure()
    fig.set_figheight(5)
    fig.set_figwidth(14)
    # Calculate color limits
    vmin = 1
    vmax = events.shape[0] // 1000 + 100
    # Iterate through all buses
    for i, bus in enumerate(range(3, 6)):
        events_bus = events[events.Bus == bus]
        # Calculate number of grid and wire events in a specific bus
        wire_events = events_bus[events_bus.Ch < 80].shape[0]
        grid_events = events_bus[events_bus.Ch >= 80].shape[0]
        # Plot
        plt.subplot(1, 3, i+1)
        sub_title = 'Bus: %d, events: %d' % (bus, events_bus.shape[0])
        sub_title += '\nWire events: %d, Grid events: %d' % (wire_events, grid_events)
        fig = PHS_2D_plot_bus(fig, events_bus, sub_title, vmin, vmax)
    plt.tight_layout()
    return fig
