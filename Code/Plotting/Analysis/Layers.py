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


# =============================================================================
#                      LAYERS INVESTIGATION - FWHM
# =============================================================================

def investigate_layers_delta_ToF(df_MG, df_He3, origin_voxel):
    def linear(x, k, m):
        return k*x + m

    def Gaussian(x, a, x0, sigma):
        return a*np.exp(-(x-x0)**2/(2*sigma**2))
    # Declare parameters
    time_offset = (0.6e-3) * 1e6
    period_time = (1/14) * 1e6
    start = 16.5e3  # us
    stop = 16.9e3  # us
    number_bins = 500
    GRID = 88
    ROW = 6
    distance_He3 = 28.239
    voxel_to_distance_dict = get_distances(origin_voxel, 0)
    # Get ToF He-3
    ToF_He3 = (df_He3.ToF * (8e-9) * 1e6 + time_offset) % period_time
    # Declare storage vectors
    FWHMs = []
    distances = []
    # Iterate through all layers
    fig = plt.figure()
    plt.subplot(1, 2, 1)
    for layer in range(0, 20):
        print(layer)
        # Filter data so that we are only using data from a single voxel
        MG_red = df_MG[((df_MG.wCh % 20) == layer) &
                    (((df_MG.Bus * 4) + df_MG.wCh//20) == ROW) &
                    (df_MG.gCh == GRID)]
        ToF_MG = (MG_red.ToF * (62.5e-9) * 1e6 + time_offset) % period_time
        hist_MG, bins_MG = np.histogram(ToF_MG, range=[start, stop], bins=number_bins)
        bin_centers_MG = 0.5 * (bins_MG[1:] + bins_MG[:-1])
        plt.errorbar(bin_centers_MG, hist_MG, np.sqrt(hist_MG), fmt='.-',
                     capsize=5, zorder=5, label='MG, layer: %d' % layer)
        # Get background level
        background_MG_events = hist_MG[(bin_centers_MG >= 16.5e3) & (bin_centers_MG <= 16.6e3)]
        background_MG_per_bin = sum(background_MG_events)/len(background_MG_events)
        plt.axhline(y=background_MG_per_bin, color='green', linewidth=2,
                    label=None, linestyle='--')
        # Get FWHM and plot it
        FWHM, idx1, idx2, max_idx = calculate_FWHM(bin_centers_MG, hist_MG)
        plt.plot(bin_centers_MG[max_idx], hist_MG[max_idx], 'o', color='red')
        plt.plot(bin_centers_MG[[idx1, idx2]],
                 [hist_MG[max_idx]/2, hist_MG[max_idx]/2], color='red')
        # Use fitting procedure to get parameters instead
        middle = bin_centers_MG[max_idx]
        hist_fit, bin_edges_fit = np.histogram(ToF_MG, range=[middle-10, middle+10], bins=50)
        bins_fit =  0.5 * (bin_edges_fit[1:] + bin_edges_fit[:-1])
        a_guess, x0_guess, sigma_guess = get_fit_parameters_guesses(hist_fit, bins_fit)
        a, x0, sigma, *_ = fit_data(hist_fit, bins_fit, a_guess, x0_guess, sigma_guess)
        FWHM_fit = 2 * np.sqrt(2*np.log(2)) * sigma
        xx = np.linspace(middle-10, middle+10, 1000)
        norm = hist_MG[max_idx]/a
        plt.plot(xx, Gaussian(xx, a, x0, sigma)*norm, color='black',
                 linestyle='-', label=None, zorder=50)
        # Store important values
        distance = voxel_to_distance_dict[MG_red.Bus, MG_red.gCh, MG_red.wCh][0]
        distances.append(distance)
        FWHMs.append(FWHM_fit)
        print('MG velocity [cm/µs], layer %d: %f' % (layer, (distance*100)/bin_centers_MG[max_idx]))
    # Plot He-3
    hist_He3, bins_He3 = np.histogram(ToF_He3, range=[start, stop], bins=number_bins)
    bin_centers_He3 = 0.5 * (bins_He3[1:] + bins_He3[:-1])
    # Calculate FWHM
    FWHM_He3, idx1, idx2, max_idx = calculate_FWHM(bin_centers_He3, hist_He3)
    # Fit data for He-3 instead
    middle = bin_centers_He3[max_idx]
    hist_fit, bin_edges_fit = np.histogram(ToF_He3, range=[middle-10, middle+10], bins=50)
    bins_fit =  0.5 * (bin_edges_fit[1:] + bin_edges_fit[:-1])
    a_guess, x0_guess, sigma_guess = get_fit_parameters_guesses(hist_fit, bins_fit)
    a, x0, sigma, *_ = fit_data(hist_fit, bins_fit, a_guess, x0_guess, sigma_guess)
    FWHM_he3_fit = 2 * np.sqrt(2*np.log(2)) * sigma
    # Plot Gaussian
    xx = np.linspace(middle-10, middle+10, 1000)
    norm = hist_He3[max_idx]/a
    plt.plot(xx, Gaussian(xx, a, x0, sigma)*norm, color='black', linestyle='-', label=None, zorder=50)
    # Plot data
    plt.errorbar(bin_centers_He3, hist_He3, np.sqrt(hist_He3), fmt='.-',
                 capsize=5, zorder=5, label='He-3', color='red')
    #print('He-3 velocity [cm/µs]: %f' % (distance_He3/100)/bin_centers_He3[max_idx])
    background_He3_events = hist_He3[(bin_centers_He3 >= 16.8e3) & (bin_centers_He3 <= 16.9e3)]
    background_He3_per_bin = sum(background_He3_events)/len(background_He3_events)
    plt.axhline(y=background_He3_per_bin, color='red', linewidth=2,
                label=None, linestyle='--')
    print('Distance in cm: %f' % (distance_He3*100))
    print('Flight time in us: %f' % bin_centers_He3[max_idx])
    print((distance_He3*100)/bin_centers_He3[max_idx])
    plt.legend(loc=2)
    plt.ylim(1, 70e3)
    plt.xlabel('ToF [µs]')
    plt.ylabel('Counts')
    plt.yscale('log')
    plt.title('Comparison ∆ToF')
    plt.grid(True, which='major', linestyle='--', zorder=0)
    plt.grid(True, which='minor', linestyle='--', zorder=0)
    # Plot FWHM
    plt.subplot(1, 2, 2)
    plt.xlabel('Distance [m]')
    plt.ylabel('FWHM [µs]')
    plt.title('FWHM ∆ToF in peaks as a function of distance')
    plt.plot(distances, FWHMs, label='MG', marker='o', linestyle='', color='blue', zorder=5)
    plt.plot(distance_He3, FWHM_he3_fit, label='He-3', marker='x', color='red', linestyle='', zorder=5)
    # Fit MG data
    fit_parameters = np.polyfit(distances, FWHMs, 1)
    k, m = fit_parameters[0], fit_parameters[1]
    xx = np.linspace(28, 30, 100)
    plt.plot(xx, linear(xx, k, m), color='black', linestyle='--', label='Multi-Grid Fit')
    plt.grid(True, which='major', linestyle='--', zorder=0)
    plt.grid(True, which='minor', linestyle='--', zorder=0)
    plt.legend(loc=2)
    fig.show()


# =============================================================================
#                               HELPER FUNCTIONS
# =============================================================================

def find_nearest(array, value):
    """
    Returns the index of the element in 'array' which is closest to 'value'.

    Args:
        array (numpy array): Numpy array with elements
        value (float): Value which we want to find the closest element to in
                       arrray

    Returns:
        idx (int): index of the element in 'array' which is closest to 'value'
    """
    idx = (np.abs(array - value)).argmin()
    return idx


def calculate_FWHM(bins, hist):
    # Extract relavant parameters
    maximum = max(hist)
    maximum_idx = find_nearest(hist, maximum)
    half_maximum = maximum/2
    half_maximum_idx_1 = find_nearest(hist[:maximum_idx], half_maximum)
    half_maximum_idx_2 = find_nearest(hist[maximum_idx:], half_maximum) + maximum_idx
    FWHM = bins[half_maximum_idx_2] - bins[half_maximum_idx_1]
    return FWHM, half_maximum_idx_1, half_maximum_idx_2, maximum_idx
