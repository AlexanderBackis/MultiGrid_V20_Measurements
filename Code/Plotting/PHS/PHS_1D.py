#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PHS_1D.py: Module which plots the PHS spectra for wires and grids separatly.
           Plots PHS from indivdual events and summation from clusters overlaid.
"""

import matplotlib.pyplot as plt
import os


# ============================================================================
#                                   PHS (1D)
# ============================================================================


def PHS_1D_plot(events, clusters, number_bins):
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
    fig = plt.figure()
    titles = ['Wires', 'Grids']
    limits = [[0, 79], [80, 119]]
    ADC_types = ['wADC', 'gADC']
    print('Number bins: %d' % number_bins)
    for i, (title, limit, ADC_type) in enumerate(zip(titles, limits, ADC_types)):
        plt.subplot(1, 2, i+1)
        # Set figure properties
        plt.title('PHS (1D) - %s' % title)
        plt.xlabel('Collected charge [ADC channels]')
        plt.ylabel('Counts')
        plt.yscale('log')
        plt.grid(True, which='major', linestyle='--', zorder=0)
        plt.grid(True, which='minor', linestyle='--', zorder=0)
        # Separate wires and grids
        #indices_to_use = ((events['Ch'] >= limit[0]) & (events['Ch'] <= limit[1]))
        # Calculate histogram
        #plt.hist(events[indices_to_use].ADC, bins=number_bins, range=[0, 4400],
        #         histtype='step', color='blue', label='Events', zorder=5)
        plt.hist(clusters[ADC_type], bins=number_bins,
                 histtype='step', color='red', label='Clusters', zorder=5)
        print(clusters[ADC_type])
        plt.legend()
    plt.tight_layout()
    return fig
