#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Coincidences_Projections.py: Helper functions for handling of paths and folders.
"""

# =============================================================================
# Coincidence Histogram (Front, Top, Side)
# =============================================================================

def Coincidences_Front_Top_Side_plot(df, data_sets, module_order,
                                     number_of_detectors, window):
    # Ensure we only plot coincident events
    df = df[(df.wCh != -1) & (df.gCh != -1)]
    df = filter_ce_clusters(window, df)
    # Define figure and set figure properties
    fig = plt.figure()
    title = ('Coincident events (Front, Top, Side)' +
             '\nData set(s): %s' % data_sets
             )
    height = 4
    width = 14
    fig = set_figure_properties(fig, title, height, width)
    if df.shape[0] != 0:
        vmin = 1
        vmax = df.shape[0] // 200 + 5
    else:
        vmin = 1
        vmax = 1
    # Retrieve text path
    dir_name = os.path.dirname(__file__)
    folder_path = os.path.join(dir_name,
                               '../../text/%s/Projections_2D/' % window.data_sets)
    mkdir_p(folder_path)
    # Plot front view
    plt.subplot(1, 3, 1)
    __, h_front = plot_2D_Front(module_order, df, fig, number_of_detectors, vmin, vmax)
    path_front = folder_path + 'Front.txt'
    if window.createText.isChecked():
        np.savetxt(path_front, h_front, fmt="%d", delimiter=",")
    # Plot top view
    plt.subplot(1, 3, 2)
    __, h_top = plot_2D_Top(module_order, df, fig, number_of_detectors, vmin, vmax)
    path_top = folder_path + 'Top.txt'
    if window.createText.isChecked():
        np.savetxt(path_top, h_top, fmt="%d", delimiter=",")
    # Plot side view
    plt.subplot(1, 3, 3)
    __, h_side = plot_2D_Side(module_order, df, fig, number_of_detectors, vmin, vmax)
    path_side = folder_path + 'Side.txt'
    if window.createText.isChecked():
        np.savetxt(path_side, h_side, fmt="%d", delimiter=",")
    return fig


# =============================================================================
# Coincidence Histogram - Front
# =============================================================================


def plot_2D_Front(bus_vec, df, fig, number_of_detectors, vmin, vmax):
    df_tot = pd.DataFrame()
    for i, bus in enumerate(bus_vec):
        df_clu = df[df.Bus == bus]
        df_clu['wCh'] += (80 * i) + (i // 3) * 80
        df_clu['gCh'] += (-80 + 1)
        df_tot = pd.concat([df_tot, df_clu])
    h, *_ = plt.hist2d(np.floor(df_tot['wCh'] / 20).astype(int) + 1,
                       df_tot.gCh,
                       bins=[12*number_of_detectors + 8, 40],
                       range=[[0.5, 12*number_of_detectors + 0.5 + 8],
                              [0.5, 40.5]
                              ],
                       norm=LogNorm(), cmap='jet', vmin=vmin, vmax=vmax
                       )
    title = 'Front view'
    locs_x = [1, 12, 17, 28, 33, 44]
    ticks_x = [1, 12, 13, 25, 26, 38]
    xlabel = 'Layer'
    ylabel = 'Grid'
    fig = stylize(fig, xlabel, ylabel, title=title, colorbar=True,
                  locs_x=locs_x, ticks_x=ticks_x)
    plt.colorbar()
    return fig, h

# =============================================================================
# Coincidence Histogram - Top
# =============================================================================


def plot_2D_Top(bus_vec, df, fig, number_of_detectors, vmin, vmax):
    df_tot = pd.DataFrame()
    for i, bus in enumerate(bus_vec):
        df_clu = df[df.Bus == bus]
        df_clu['wCh'] += (80 * i) + (i // 3) * 80
        df_tot = pd.concat([df_tot, df_clu])
    h, *_ = plt.hist2d(np.floor(df_tot['wCh'] / 20).astype(int) + 1,
                       df_tot['wCh'] % 20 + 1,
                       bins=[12*number_of_detectors + 8, 20],
                       range=[[0.5, 12*number_of_detectors + 0.5 + 8],
                              [0.5, 20.5]
                              ],
                       norm=LogNorm(), cmap='jet', vmin=vmin, vmax=vmax
                       )
    title = 'Top view'
    locs_x = [1, 12, 17, 28, 33, 44]
    ticks_x = [1, 12, 13, 25, 26, 38]
    xlabel = 'Layer'
    ylabel = 'Wire'
    fig = stylize(fig, xlabel, ylabel, title=title, colorbar=True,
                  locs_x=locs_x, ticks_x=ticks_x)
    plt.colorbar()
    return fig, h


# =============================================================================
# Coincidence Histogram - Side
# =============================================================================


def plot_2D_Side(bus_vec, df, fig, number_of_detectors, vmin, vmax):

    df_tot = pd.DataFrame()
    for i, bus in enumerate(bus_vec):
        df_clu = df[df.Bus == bus]
        df_clu['gCh'] += (-80 + 1)
        df_tot = pd.concat([df_tot, df_clu])
    h, *_ = plt.hist2d(df_tot['wCh'] % 20 + 1, df_tot['gCh'],
                       bins=[20, 40],
                       range=[[0.5, 20.5], [0.5, 40.5]],
                       norm=LogNorm(),
                       cmap='jet', vmin=vmin, vmax=vmax
                       )

    title = 'Side view'
    xlabel = 'Wire'
    ylabel = 'Grid'
    fig = stylize(fig, xlabel, ylabel, title=title, colorbar=True)
    plt.colorbar()
    return fig, h
