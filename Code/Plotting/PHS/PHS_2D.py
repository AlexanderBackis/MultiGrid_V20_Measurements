#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PHS_2D.py: Helper functions for handling of paths and folders.
"""
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# =============================================================================
#                                   PHS (2D)
# =============================================================================


def PHS_2D_plot(events):
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
    fig = plt.figure()
    fig.set_figheight(12)
    fig.set_figwidth(14)
    vmin = 1
    vmax = events.shape[0] // 1000 + 100
    for bus in range(0, 9):
        events_bus = events[events.Bus == bus]
        wire_events = events_bus[events_bus.Ch < 80].shape[0]
        grid_events = events_bus[events_bus.Ch >= 80].shape[0]
        plt.subplot(3, 3, bus+1)
        sub_title = 'Bus: %d, events: %d' % (bus, events_bus.shape[0])
        sub_title += '\nWire events: %d, Grid events: %d' % (wire_events, grid_events)
        fig = PHS_2D_plot_bus(fig, events_bus, sub_title, vmin, vmax)
    plt.tight_layout()
    return fig
