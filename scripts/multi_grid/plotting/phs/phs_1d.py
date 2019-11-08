#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" PHS_1D.py: Module which plots the PHS spectra for wires and grids separatly.
               Plots PHS from indivdual events and summation from clusters
               overlaid.
"""

import matplotlib.pyplot as plt
import numpy as np


# ============================================================================
#                                   PHS (1D)
# ============================================================================


def PHS_1D_plot(clusters, number_bins, label='', norm=1, ylabel='', intervals=None):
    """
    Histograms the ADC-values from wires and grids individually, and overlays
    the results from indiviual events and clustered events. For the clustered
    events, each data point is a summation of the ADC values contained in that
    cluster, i.e. the summed ADC is being histogrammed.

    Args:
        events (DataFrame): Individual events
        clusters (DataFrame): Clustered events
        number_bins (int): Number of bins to use in the histogram

    Returns:
        fig (Figure): Figure containing 1D PHS plot
    """
    # Declare parameters
    titles = ['Wires', 'Grids']
    limits = [[0, 79], [80, 119]]
    ADC_types = ['wADC', 'gADC']
    weights = norm*np.ones(clusters.shape[0])
    bin_centers_vec = []
    hists = []
    for i, (title, limit, ADC_type) in enumerate(zip(titles, limits, ADC_types)):
        plt.subplot(1, 2, i+1)
        # Set figure properties
        plt.title('PHS (1D) - %s' % title)
        plt.xlabel('Collected charge [ADC channels]')
        plt.ylabel('Counts %s' % ylabel)
        plt.yscale('log')
        plt.grid(True, which='major', linestyle='--', zorder=0)
        plt.grid(True, which='minor', linestyle='--', zorder=0)
        if intervals is not None:
            range = intervals[i]
        else:
            range = None
        # Plot
        hist, bins, *_ = plt.hist(clusters[ADC_type], bins=number_bins,
                                  histtype='step', label='Clusters %s' % label,
                                  zorder=5, weights=weights, range=range)
        bin_centers = 0.5 * (bins[1:] + bins[:-1])
        bin_centers_vec.append(bin_centers)
        hists.append(hist)
        plt.legend()
    plt.tight_layout()
    return bin_centers_vec, hists
