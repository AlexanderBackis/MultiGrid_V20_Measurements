#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PHS_2D.py: Helper functions for handling of paths and folders.
"""


# =============================================================================
#                                   PHS (2D)
# =============================================================================


def PHS_2D_plot(events, data_sets, module_order):
    def PHS_2D_plot_bus(fig, events, sub_title, vmin, vmax):
        xlabel = 'Channel'
        ylabel = 'Charge [ADC channels]'
        bins = [120, 120]
        fig = stylize(fig, xlabel, ylabel, title=sub_title, colorbar=True)
        plt.hist2d(events.Channel, events.ADC, bins=bins, norm=LogNorm(),
                   range=[[-0.5, 119.5], [0, 4400]], vmin=vmin, vmax=vmax,
                   cmap='jet')
        plt.colorbar()
        return fig

    fig = plt.figure()
    title = 'PHS (2D)\nData set(s): %s' % data_sets
    height = 12
    width = 14
    vmin = 1
    vmax = events.shape[0] // 1000 + 100
    for i, bus in enumerate(module_order):
        events_bus = events[events.Bus == bus]
        wire_events = events_bus[events_bus.Channel < 80].shape[0]
        grid_events = events_bus[events_bus.Channel >= 80].shape[0]
        plt.subplot(3, 3, i+1)
        sub_title = 'Bus: %d, events: %d' % (bus, events_bus.shape[0])
        sub_title += ('\nWire events: %d, Grid events: %d'
                      % (wire_events, grid_events)
                      )
        fig = PHS_2D_plot_bus(fig, events_bus, sub_title, vmin, vmax)
    fig = set_figure_properties(fig, title, height, width)
    return fig
