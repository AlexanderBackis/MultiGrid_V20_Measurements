#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Animation3D.py: Helper functions for handling of paths and folders.
"""

import plotly as py
import numpy as np
import matplotlib.pyplot as plt
import os
import imageio
import shutil

from Plotting.Misc.Timestamp import timestamp_plot
from HelperFunctions.PathsAndFolders import mkdir_p

from HelperFunctions.EnergyTransfer import calculate_energy, get_distances
from HelperFunctions.Fitting import fit_data, get_fit_parameters_guesses, get_hist


# =============================================================================
#                      LAYERS INVESTIGATION - ToF
# =============================================================================

def investigate_layers_ToF(df):
    # Define parameters for ToF-histogram
    time_offset = (0.6e-3) * 1e6
    period_time = (1/14) * 1e6
    number_bins = 200
    # Define beam hit position
    GRID = 88
    ROW = 6
    # Iterate through the first ten layers and compare
    fig = plt.figure()
    for layer in range(0, 10):
        # Filter data so that we are only using data from a single voxel
        df_red = df[((df.wCh % 20) == layer) &
                    (((df.Bus * 4) + df.wCh//20) == ROW) &
                    (df.gCh == GRID)]
        # Plot ToF-histogram
        plt.hist((df_red.ToF * 62.5e-9 * 1e6 + time_offset) % period_time,
                 bins=number_bins,
                 range=[40500, 40900],
                 zorder=4,
                 histtype='step',
                 label='Layer: %d' % layer)
        plt.grid(True, which='major', linestyle='--', zorder=0)
        plt.grid(True, which='minor', linestyle='--', zorder=0)
        plt.xlabel('ToF [µs]')
        plt.ylabel('Counts')
        plt.ylim(0, 300)
        plt.title('ToF from different layers (gCh = 88, row = 6)')
    plt.legend()
    fig.show()


# =============================================================================
#                      LAYERS INVESTIGATION - FWHM
# =============================================================================

def investigate_layers_FWHM(df, origin_voxel):
    def linear(x, k, m):
        return k*x + m
    # Declare parameters
    number_bins = 100
    GRID = 88
    ROW = 6
    voxel_to_distance_dict = get_distances(origin_voxel)
    # Define region of the peak we want to study
    peak = 13.32 # meV
    left, right = peak - 0.07, peak + 0.07
    # Iterate through all layers, a single voxel from each layer, saving the FWHMs
    fig = plt.figure()
    plt.subplot(1, 2, 1)
    FWHMs = []
    distances = []
    for layer in range(0, 20):
        # Filter data so that we are only using data from a single voxel
        df_red = df[((df.wCh % 20) == layer) &
                    (((df.Bus * 4) + df.wCh//20) == ROW) &
                    (df.gCh == GRID)]
        print(df_red)
        # Caluclate energies
        energies = calculate_energy(df_red, origin_voxel)
        hist, bins = get_hist(energies, number_bins, left, right)
        a_guess, x0_guess, sigma_guess = get_fit_parameters_guesses(hist, bins)
        a, x0, sigma, *_ = fit_data(hist, bins, a_guess, x0_guess, sigma_guess)
        # Plot
        plt.plot(bins, hist, marker='.', linestyle='-', zorder=5, label='Layer: %d' % layer)
        # Extract important parameters
        FWHM = 2 * np.sqrt(2*np.log(2)) * sigma
        FWHMs.append(FWHM)
        distance = voxel_to_distance_dict[df_red.Bus, df_red.gCh, df_red.wCh][0]
        distances.append(distance)
    # Stylize plot
    plt.grid(True, which='major', linestyle='--', zorder=0)
    plt.grid(True, which='minor', linestyle='--', zorder=0)
    plt.xlabel('Energy [meV]')
    plt.ylabel('Counts')
    plt.title('Energy Histogram, peak at 13.32 meV')
    plt.legend()
    # Plot second subplot, containing FWHM for the different layers, as well as
    # linear fit on MG data.
    fit_parameters = np.polyfit(distances, FWHMs, 1)
    k, m = fit_parameters[0], fit_parameters[1]
    plt.subplot(1, 2, 2)
    plt.title('FWHM vs distance')
    plt.ylabel('FWHM [meV]')
    plt.xlabel('Distance to Source Chopper [m]')
    plt.grid(True, which='major', linestyle='--', zorder=0)
    plt.grid(True, which='minor', linestyle='--', zorder=0)
    for distance, FWHM in zip(distances, FWHMs):
        plt.plot(distance, FWHM, marker='o', zorder=5, label=None)
    plt.plot(28.239, 0.015, marker='x', color='red', zorder=5, label='He-3')
    xx = np.linspace(28, 30, 100)
    plt.plot(xx, linear(xx, k, m), color='black', linestyle='--', label='Multi-Grid Fit')
    plt.legend()
    fig.show()
