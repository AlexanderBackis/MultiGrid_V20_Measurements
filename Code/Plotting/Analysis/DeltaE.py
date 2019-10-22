#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
DeltaE.py: Function which histograms energy transfer data
"""

from HelperFunctions.EnergyTransfer import calculate_energy
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks

# =============================================================================
#                             ENERGY TRANSFER
# =============================================================================


def energy_plot(df, origin_voxel, number_bins, start=1, stop=10,
                plot_energy=False, label=None, useMaxNorm=False):
    """
    Histograms the energy transfer values from a measurement

    Args:
        df (DataFrame): Clustered events
        Ei (float): Incident energy in meV
        number_bins (int): Number of bins to histogram energy transfer data

    Returns:
        fig (Figure): Figure containing nine 2D coincidences histograms, one
                      for each bus.
        dE_hist (numpy array): Numpy array containing the histogram data
        bin_centers (numpy array): Numpy array containing the bin centers
    """
    def meV_to_A(energy):
        return np.sqrt(81.81/energy)

    def A_to_meV(wavelength):
        return (81.81/(wavelength ** 2))
    # Calculate DeltaE
    energy = calculate_energy(df, origin_voxel)
    # Select normalization
    if useMaxNorm is False:
        norm = 1 * np.ones(len(energy))
    else:
        if plot_energy:
            hist_temp, _ = np.histogram(energy, bins=number_bins, range=[A_to_meV(stop), A_to_meV(start)])
            norm = (1/max(hist_temp))*np.ones(len(energy))
        else:
            hist_temp, _ = np.histogram(meV_to_A(energy), bins=number_bins, range=[start, stop])
            norm = (1/max(hist_temp))*np.ones(len(energy))
    # Plot data
    if plot_energy:
        plt.xlabel('Energy [meV]')
        plt.title('Energy Distribution')
        plt.xscale('log')
        hist, bin_edges, *_ = plt.hist(energy, bins=number_bins,
                                       range=[A_to_meV(stop), A_to_meV(start)],
                                       zorder=5, histtype='step',
                                       label=label,
                                       weights=norm)
    else:
        plt.xlabel('Wavelength [Å]')
        plt.title('Wavelength Distribution')
        hist, bin_edges, *_ = plt.hist(meV_to_A(energy), bins=number_bins,
                                       range=[start, stop], zorder=5,
                                       histtype='step',
                                       label=label,
                                       weights=norm)
    bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
    plt.grid(True, which='major', linestyle='--', zorder=0)
    plt.grid(True, which='minor', linestyle='--', zorder=0)
    plt.ylabel('Counts')
    plt.yscale('log')
    return hist, bin_centers
