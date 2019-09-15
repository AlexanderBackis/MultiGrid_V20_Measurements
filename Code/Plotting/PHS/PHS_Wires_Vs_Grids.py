#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PHS_Wires_Vs_Grids.py: Helper functions for handling of paths and folders.
"""

# =============================================================================
#                           PHS (Wires vs Grids)
# =============================================================================


def PHS_wires_vs_grids_plot(ce, data_sets, module_order, window):
    def charge_scatter(fig, ce, sub_title, bus, vmin, vmax):
        xlabel = 'Collected charge wires [ADC channels]'
        ylabel = 'Collected charge grids [ADC channels]'
        bins = [200, 200]
        ADC_range = [[0, 5000], [0, 5000]]
        fig = stylize(fig, xlabel, ylabel, title=sub_title, colorbar=True)
        plt.hist2d(ce.wADC, ce.gADC, bins=bins,
                   norm=LogNorm(), range=ADC_range,
                   vmin=vmin, vmax=vmax, cmap='jet')
        plt.colorbar()
        return fig
    # Intial filter
    ce = filter_ce_clusters(window, ce)
    # Plot data
    fig = plt.figure()
    title = 'PHS (Wires vs Grids)\nData set(s): %s' % data_sets
    height = 12
    width = 14
    if ce.shape[0] != 0:
        vmin = 1
        vmax = ce.shape[0] // 4500 + 1000
    else:
        vmin = 1
        vmax = 1
    for i, bus in enumerate(module_order):
        events_bus = ce[ce.Bus == bus]
        sub_title = 'Bus %d\n(%d events)' % (bus, events_bus.shape[0])
        plt.subplot(3, 3, i+1)
        fig = charge_scatter(fig, events_bus, sub_title, bus, vmin, vmax)
    fig = set_figure_properties(fig, title, height, width)
    return fig
