#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ImportHe3.py: Imports He3 data taken using the MCA4 Multichannel Analyzer
"""

import os
import struct
import shutil
import zipfile
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from HeliumTubes.EnergyHe3 import calculate_He3_energy
from scipy.signal import find_peaks

# =============================================================================
#                               PHS - HELIUM-3
# =============================================================================

def He3_PHS_plot(df, number_bins):
    plt.hist(df['ADC'], histtype='step', color='blue', zorder=5,
             bins=number_bins)
    plt.grid(True, which='major', linestyle='--', zorder=0)
    plt.grid(True, which='minor', linestyle='--', zorder=0)
    plt.xlabel('Charge [ADC Channels]')
    plt.ylabel('Counts')
    plt.title('PHS')

# =============================================================================
#                         CHANNEL HISTOGRAM - HELIUM-3
# =============================================================================

def He3_Ch_plot(df):
    plt.hist(df['Ch'], histtype='step', color='red', zorder=5, bins=20)
    plt.grid(True, which='major', linestyle='--', zorder=0)
    plt.grid(True, which='minor', linestyle='--', zorder=0)
    plt.xlabel('Channel')
    plt.ylabel('Counts')
    plt.title('Channel')

# =============================================================================
#                         TOF HISTOGRAM - HELIUM-3
# =============================================================================

def He3_ToF_plot(df, number_bins, label=None):
    # Declare parameters
    time_offset = (0.6e-3) * 1e6
    period_time = (1/14) * 1e6
    plt.hist((df.ToF * (8e-9) * 1e6 + time_offset) % period_time,
             histtype='step', zorder=5, bins=number_bins, label=label)
    plt.xlabel('ToF [µs]')
    plt.ylabel('Counts')
    plt.title('ToF')
    plt.grid(True, which='major', linestyle='--', zorder=0)
    plt.grid(True, which='minor', linestyle='--', zorder=0)


# =============================================================================
#                             ENERGY - HELIUM-3
# =============================================================================


def energy_plot_He3(df, number_bins, plot_energy=False, label=None, useMaxNorm=False):
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
    energy = calculate_He3_energy(df)
    # Define lambda range
    start = 1
    end = 10
    # Select normalization
    if useMaxNorm is False:
        norm = 1 * np.ones(len(energy))
    else:
        if plot_energy:
            hist_temp, _ = np.histogram(energy, bins=number_bins, range=[A_to_meV(end), A_to_meV(start)])
            norm = (1/max(hist_temp))*np.ones(len(energy))
        else:
            hist_temp, _ = np.histogram(meV_to_A(energy), bins=number_bins, range=[start, end])
            norm = (1/max(hist_temp))*np.ones(len(energy))
    # Plot data
    if plot_energy:
        plt.xlabel('Energy [meV]')
        plt.title('Energy Distribution')
        plt.xscale('log')
        hist, bin_edges, *_ = plt.hist(energy, bins=number_bins,
                                       range=[A_to_meV(end), A_to_meV(start)],
                                       zorder=5, histtype='step',
                                       label=label,
                                       weights=norm,
                                       linestyle='-')
        bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
        # Extract heights
        heights = np.zeros(len(bin_centers))
        heights[(bin_centers >= 0) & (bin_centers <= 2.5)] = 2000
        heights[(bin_centers >= 2.5) & (bin_centers <= 50.0)] = 20000
        heights[(bin_centers >= 50.0) & (bin_centers <= 70.0)] = 2000
        heights[(bin_centers >= 70.0) & (bin_centers <= 100.0)] = 286
        # Get peaks
        peaks, *_ = find_peaks(hist, height=heights)
        plt.plot(bin_centers[peaks], hist[peaks], marker='x', linestyle='', color='red')
        plt.plot(bin_centers, heights, color='black')
    else:
        plt.xlabel('Wavelength [Å]')
        plt.title('Wavelength Distribution')
        hist, bin_edges, *_ = plt.hist(meV_to_A(energy), bins=number_bins,
                                       range=[start, end], zorder=5,
                                       histtype='step',
                                       label=label,
                                       weights=norm,
                                       linestyle='-')
    bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
    plt.grid(True, which='major', linestyle='--', zorder=0)
    plt.grid(True, which='minor', linestyle='--', zorder=0)
    plt.ylabel('Counts')
    plt.yscale('log')
    return hist, bin_centers
