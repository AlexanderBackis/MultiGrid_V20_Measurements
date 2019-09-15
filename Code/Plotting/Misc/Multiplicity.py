#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Multiplicity.py: Helper functions for handling of paths and folders.
"""

# =============================================================================
#                               Multiplicity
# =============================================================================

def Multiplicity_plot(df, data_sets, module_order, number_of_detectors,
                      window):
    # Initial filter
    df = filter_ce_clusters(window, df)
    # Declare parameters
    name = 'Multiplicity\nData set(s): ' + str(data_sets)
    fig = plt.figure()
    figwidth = 14
    figheight = 12
    if df.shape[0] != 0:
        vmin = 1
        vmax = df.shape[0] // 9 + 10
    else:
        vmin = 1
        vmax = 1
    # Prepare figure
    fig.suptitle(name, x=0.5, y=1.08)
    fig.set_figheight(figheight)
    fig.set_figwidth(figwidth)
    # Iterate through all buses
    for loc, bus in enumerate(module_order):
        plt.subplot(3, 3, loc+1)
        df_bus = df[df.Bus == bus]
        plot_multiplicity_bus(df_bus, bus, fig, vmin, vmax)

    plt.tight_layout()
    fig.show()


def plot_multiplicity_bus(df, bus, fig, vmin, vmax):
    # Declare parameters
    m_range = [0, 80, 0, 40]
    # Plot data
    hist, xbins, ybins, im = plt.hist2d(df.wM, df.gM,
                                        bins=[m_range[1]-m_range[0]+1,
                                              m_range[3]-m_range[2]+1],
                                        range=[[m_range[0], m_range[1]+1],
                                               [m_range[2], m_range[3]+1]],
                                        norm=LogNorm(),
                                        vmin=vmin,
                                        vmax=vmax,
                                        cmap='jet')
    # Iterate through all squares and write percentages
    #tot = df.shape[0]
    #font_size = 12
    #for i in range(len(ybins)-1):
    #    for j in range(len(xbins)-1):
    #        if hist[j, i] > 0:
    #            text = plt.text(xbins[j]+0.5, ybins[i]+0.5,
    #                            '%.1f%%' % (100*(hist[j, i]/tot)),
    #                            color="w", ha="center", va="center",
    #                            fontweight="bold", fontsize=font_size)
    #            text.set_path_effects([path_effects.Stroke(linewidth=1,
    #                                                       foreground='black'),
    #                                   path_effects.Normal()])
    # Set ticks on axis
    ticks_x = np.arange(m_range[0], m_range[1]+1, 1)
    locs_x = np.arange(m_range[0] + 0.5, m_range[1]+1.5, 1)
    ticks_y = np.arange(m_range[2], m_range[3]+1, 1)
    locs_y = np.arange(m_range[2] + 0.5, m_range[3]+1.5, 1)
    #plt.xticks(locs_x, ticks_x)
    #plt.yticks(locs_y, ticks_y)
    plt.xlabel("Wire Multiplicity")
    plt.ylabel("Grid Multiplicity")
    plt.colorbar()
    plt.tight_layout()
    name = 'Bus ' + str(bus) + '\n(' + str(df.shape[0]) + ' events)'
    plt.title(name)
