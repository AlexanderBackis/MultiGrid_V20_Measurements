#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Coincidences_2D.py: Histograms the hitposition, expressed in wire- and grid
                    channel coincidences in all of the detector voxels.
"""

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from HelperFunctions.AreaAndSolidAngle import get_multi_grid_area_and_solid_angle

# =============================================================================
#                          Coincidence Histogram (2D)
# =============================================================================

def coincidences_2D_plot(ce, measurement_time, bus_start, bus_stop):
    """
    Histograms the hitposition, expressed in wire- and grid channel coincidences
    in all of the detector voxels.

    Args:
        ce (DataFrame): Clustered events
        measurement_time (float): Measurement time expressed in seconds

    Returns:
        fig (Figure): Figure containing nine 2D coincidences histograms, one
                      for each bus.
        histograms (numpy array): 2D numpy array containing a 2D matrix.
                                  Columns are grids, rows are wires, and the
                                  values in the cells are the number of counts
                                  in that specific voxel.
    """

    def plot_2D_bus(fig, sub_title, ce, vmin, vmax, duration):
        h, *_ = plt.hist2d(ce.wCh, ce.gCh, bins=[80, 40],
                           range=[[-0.5, 79.5], [79.5, 119.5]],
                           vmin=vmin, vmax=vmax,
                           norm=LogNorm(), cmap='jet')
        plt.xlabel('Wire [Channel number]')
        plt.ylabel('Grid [Channel number]')
        plt.title(sub_title)
        plt.colorbar()
        return fig, h

    # Perform initial filter
    ce = ce[(ce.wCh != -1) & (ce.gCh != -1)]
    # Calculate color limits
    if ce.shape[0] != 0:
        duration = measurement_time
        vmin = 1
        vmax = ce.shape[0] // 450 + 5
    else:
        duration = 1
        vmin = 1
        vmax = 1
    # Plot data
    fig = plt.figure()
    number_detectors = ((bus_stop + 1) - bus_start)//3
    fig.set_figheight(5*number_detectors)
    fig.set_figwidth(17)
    histograms = []
    for i, bus in enumerate(range(bus_start, bus_stop+1)):
        ce_bus = ce[ce.Bus == bus]
        # Calculate number of events and rate in a specific bus
        number_events = ce_bus.shape[0]
        events_per_s = number_events/duration
        events_per_m_s = events_per_s
        sub_title = ('Bus %d\n(%d events, %.3f events/s)' % (bus, number_events, events_per_s))
        plt.subplot(number_detectors, 3, i+1)
        fig, h = plot_2D_bus(fig, sub_title, ce_bus, vmin, vmax, duration)
        histograms.append(h)
    plt.tight_layout()
    return fig, histograms
