#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
DeltaE.py: Function which histograms energy transfer data
"""

from HelperFunctions.EnergyTransfer import calculate_energy
import matplotlib.pyplot as plt
import numpy as np

# =============================================================================
#                             ENERGY TRANSFER
# =============================================================================


def energy_plot(df, detector_type, origin_voxel, number_bins):
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
    energy = calculate_energy(df, detector_type, origin_voxel)
    # Histogram DeltaE
    print(energy)
    # Plot data
    fig = plt.figure()
    plt.hist(meV_to_A(energy), bins=number_bins, range=[0, 10], zorder=5,
             histtype='step', color='black')
    plt.grid(True, which='major', linestyle='--', zorder=0)
    plt.grid(True, which='minor', linestyle='--', zorder=0)
    plt.xlabel('Wavelength [Å]')
    plt.ylabel('Counts')
    plt.title('Wavelength Distribution')
    return fig
