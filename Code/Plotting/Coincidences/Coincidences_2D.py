#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PathsAndFolders.py: Helper functions for handling of paths and folders.
"""


# =============================================================================
# Coincidence Histogram (2D)
# =============================================================================

def Coincidences_2D_plot(ce, data_sets, module_order, window):
    def plot_2D_bus(fig, sub_title, ce, vmin, vmax, duration):
        h, *_ = plt.hist2d(ce.wCh, ce.gCh, bins=[80, 40],
                           range=[[-0.5, 79.5], [79.5, 119.5]],
                           vmin=vmin, vmax=vmax, norm=LogNorm(), cmap='jet')
        xlabel = 'Wire [Channel number]'
        ylabel = 'Grid [Channel number]'
        fig = stylize(fig, xlabel, ylabel, title=sub_title, colorbar=True)
        plt.colorbar()
        return fig, h

    # Filter clusters
    ce = filter_ce_clusters(window, ce)
    # Declare parameters (added with condition if empty array)
    if ce.shape[0] != 0:
        duration = window.measurement_time
        vmin = 1
        vmax = ce.shape[0] // 4500 + 5
    else:
        duration = 1
        vmin = 1
        vmax = 1
    title = 'Coincident events (2D)\nData set(s): %s' % data_sets
    height = 12
    width = 14
    # Ensure only coincident events are plotted
    ce = ce[(ce.wCh != -1) & (ce.gCh != -1)]
    # Retrieve text path
    dir_name = os.path.dirname(__file__)
    folder_path = os.path.join(dir_name,
                               '../../text/%s/Coincidences_2D/' % window.data_sets)
    mkdir_p(folder_path)
    # Plot data and save matrices in text
    fig = plt.figure()
    for i, bus in enumerate(module_order):
        ce_bus = ce[ce.Bus == bus]
        number_events = ce_bus.shape[0]
        events_per_s = round(number_events/duration, 4)
        sub_title = ('Bus %d\n(%d events, %f events/s)' % (bus, number_events,
                                                           events_per_s)
                     )
        plt.subplot(3, 3, i+1)
        fig, h = plot_2D_bus(fig, sub_title, ce_bus, vmin, vmax, duration)
        # Save matrix in text
        path = folder_path + 'Bus_' + str(bus) + '.txt'
        if window.createText.isChecked():
            np.savetxt(path, h, fmt="%d", delimiter=",")
    fig = set_figure_properties(fig, title, height, width)
    return fig
