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
from scipy.optimize import curve_fit

from multi_grid.plotting.misc.timestamp import timestamp_plot
from multi_grid.helper_functions.misc import mkdir_p
from multi_grid.helper_functions.fitting import get_hist, get_fit_parameters_guesses, fit_data
from multi_grid.helper_functions.energy_calculation import calculate_energy, get_distances

from helium_tube.energy_he3 import calculate_He3_energy



# =============================================================================
#                      LAYERS INVESTIGATION - TOF
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

def investigate_layers_FWHM(df, df_He3, origin_voxel):
    def linear(x, k, m):
        return k*x + m

    def Gaussian(x, a, x0, sigma):
        return a*np.exp(-(x-x0)**2/(2*sigma**2))

    # Declare parameters
    number_bins = 50
    GRID = 88
    ROW = 6
    voxel_to_distance_dict = get_distances(origin_voxel, 0)
    # Define region of the peak we want to study
    peak = 15.03 # meV
    left, right = peak - 0.07, peak + 0.07
    # Iterate through all layers, a single voxel from each layer, saving the FWHMs
    fig = plt.figure()
    plt.subplot(1, 3, 1)
    FWHMs = []
    distances = []
    errors = []
    for layer in range(0, 20):
        # Filter data so that we are only using data from a single voxel
        df_red = df[((df.wCh % 20) == layer) &
                    (((df.Bus * 4) + df.wCh//20) == ROW) &
                    (df.gCh == GRID)]
        print(layer)
        # Caluclate energies
        energies = calculate_energy(df_red, origin_voxel)
        hist, bins = get_hist(energies, number_bins, left, right)
        a_guess, x0_guess, sigma_guess = get_fit_parameters_guesses(hist, bins)
        fit_left, fit_right = x0_guess - 2 * sigma_guess, x0_guess + 2 * sigma_guess
        hist_fit, bins_fit = get_hist(energies[(energies >= fit_left) & (energies <= fit_right)],
                                      10, fit_left, fit_right)
        a, x0, sigma, __, __, pcov = fit_data(hist_fit, bins_fit, a_guess, x0_guess, sigma_guess)
        a_err, x0_err, sigma_err = np.sqrt(np.diag(pcov))
        # Plot
        plt.errorbar(bins, hist, np.sqrt(hist), marker='.', linestyle='-', fmt='.-', capsize=5,
                     zorder=5, label='Layer: %d' % layer)
        xx = np.linspace(fit_left, fit_right, 1000)
        norm = (max(hist)/max(Gaussian(xx, a, x0, sigma)))
        plt.plot(xx, Gaussian(xx, a, x0, sigma)*norm, color='black', linestyle='-', label=None, zorder=50)
        # Extract important parameters
        FWHM = 2 * np.sqrt(2*np.log(2)) * sigma
        FWHM_err = 2 * np.sqrt(2*np.log(2)) * sigma_err
        distance = voxel_to_distance_dict[df_red.Bus, df_red.gCh, df_red.wCh][0]
        FWHMs.append(FWHM)
        errors.append(FWHM_err)
        distances.append(distance)
    # Stylize plot
    plt.grid(True, which='major', linestyle='--', zorder=0)
    plt.grid(True, which='minor', linestyle='--', zorder=0)
    plt.xlabel('Energy [meV]')
    plt.ylabel('Counts')
    plt.title('Energy Histogram MG, peak at %.2f meV' % peak)
    plt.legend()
    # Plot He-3
    plt.subplot(1, 3, 2)
    energies_He3 = calculate_He3_energy(df_He3)
    hist, bins = get_hist(energies_He3, number_bins, left, right)
    a_guess, x0_guess, sigma_guess = get_fit_parameters_guesses(hist, bins)
    fit_left, fit_right = x0_guess - 2 * sigma_guess, x0_guess + 2 * sigma_guess
    hist_fit, bins_fit = get_hist(energies_He3[(energies_He3 >= fit_left) & (energies_He3 <= fit_right)],
                                  10, fit_left, fit_right)
    a, x0, sigma, __, __, pcov = fit_data(hist_fit, bins_fit, a_guess, x0_guess, sigma_guess)
    a_err, x0_err, sigma_err  = np.sqrt(np.diag(pcov))
    # Plot
    plt.errorbar(bins, hist, np.sqrt(hist), marker='.', fmt='.-', capsize=5, linestyle='-', zorder=5, label='He-3', color='red')
    xx = np.linspace(fit_left, fit_right, 1000)
    norm = (max(hist)/max(Gaussian(xx, a, x0, sigma)))
    plt.plot(xx, Gaussian(xx, a, x0, sigma)*norm, color='black', linestyle='-', label=None, zorder=50)
    # Extract important parameters
    FWHM_He3 = 2 * np.sqrt(2*np.log(2)) * sigma
    FWHM_He3_err = 2 * np.sqrt(2*np.log(2)) * sigma_err
    # Stylize plot
    plt.grid(True, which='major', linestyle='--', zorder=0)
    plt.grid(True, which='minor', linestyle='--', zorder=0)
    plt.xlabel('Energy [meV]')
    plt.ylabel('Counts')
    plt.title('Energy Histogram He-3, peak at %.2f meV' % peak)
    plt.legend()
    # Plot third subplot, containing FWHM for the different layers, as well as
    # linear fit on MG data.
    paras, pcov = np.polyfit(distances, FWHMs, 1, w=1/np.array(errors), cov=True)
    k, m = paras[0], paras[1]
    k_err, m_err = np.sqrt(np.diag(pcov))
    upper_k, upper_m = k - k_err, m + m_err
    lower_k, lower_m = k + k_err, m - m_err
    plt.subplot(1, 3, 3)
    plt.title('FWHM vs distance')
    plt.ylabel('FWHM [meV]')
    plt.xlabel('Distance to Source Chopper [m]')
    plt.grid(True, which='major', linestyle='--', zorder=0)
    plt.grid(True, which='minor', linestyle='--', zorder=0)
    #plt.plot(distances, FWHMs, marker='o', zorder=5, label='MG', color='blue', linestyle='')
    plt.errorbar(distances, FWHMs, errors, fmt='.-', capsize=5, zorder=5, label='MG', color='blue', linestyle='')
    plt.errorbar(28.239, FWHM_He3, FWHM_He3_err, fmt='.-', color='red', capsize=5, zorder=5, label='He-3')
    xx = np.linspace(26, 30, 100)
    plt.plot(xx, linear(xx, upper_k, upper_m), color='black', linestyle='-.', label='Multi-Grid Fit (Upper Bound)')
    plt.plot(xx, linear(xx, k, m), color='black', linestyle='-', label='Multi-Grid Fit')
    plt.plot(xx, linear(xx, lower_k, lower_m), color='black', linestyle='-.', label='Multi-Grid Fit (Lower Bound)')
    plt.xlim(28, 28.65)
    plt.ylim(0.015, 0.035)
    plt.legend()
    fig.show()


# =============================================================================
#                      LAYERS INVESTIGATION - TOF
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
    stop = 17e3  # us
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
    errors = []
    # Prepare output paths
    dirname = os.path.dirname(__file__)
    output_folder = os.path.join(dirname, '../../../../output/Layers/')
    animation_path = os.path.join(dirname, '../../../../output/layers_animation.gif')
    mkdir_p(output_folder)
    # Iterate through all layers
    for layer in range(0, 20):
        print(layer)
        fig = plt.figure()
        # Filter data so that we are only using data from a single voxel
        MG_red = df_MG[((df_MG.wCh % 20) == layer) &
                    (((df_MG.Bus * 4) + df_MG.wCh//20) == ROW) &
                    (df_MG.gCh == GRID)]
        print(MG_red)
        ToF_MG = (MG_red.ToF * (62.5e-9) * 1e6 + time_offset) % period_time
        hist_MG, bins_MG = np.histogram(ToF_MG, range=[start, stop], bins=number_bins)
        bin_centers_MG = 0.5 * (bins_MG[1:] + bins_MG[:-1])
        plt.errorbar(bin_centers_MG, hist_MG, np.sqrt(hist_MG), fmt='.-',
                     capsize=5, zorder=20, label='MG', color='blue')
        # Get background level
        background_MG_events = hist_MG[(bin_centers_MG >= 16.5e3) & (bin_centers_MG <= 16.6e3)]
        background_MG_per_bin = sum(background_MG_events)/len(background_MG_events)
        #plt.axhline(y=background_MG_per_bin, color='green', linewidth=2, label=None, linestyle='--')
        # Get FWHM and plot it
        FWHM, idx1, idx2, max_idx = calculate_FWHM(bin_centers_MG, hist_MG)
        plt.plot(bin_centers_MG[max_idx], hist_MG[max_idx], 'o', color='red')
        plt.plot(bin_centers_MG[[idx1, idx2]],
                 [hist_MG[max_idx]/2, hist_MG[max_idx]/2], color='red')
        # Use fitting procedure to get parameters instead
        middle = bin_centers_MG[max_idx]
        #hist_fit, bin_edges_fit = np.histogram(ToF_MG, range=[middle-15, middle+15], bins=50)
        #bins_fit =  0.5 * (bin_edges_fit[1:] + bin_edges_fit[:-1])
        a_guess, x0_guess, sigma_guess = get_fit_parameters_guesses(hist_MG, bin_centers_MG)
        fit_left, fit_right = x0_guess - sigma_guess*2, x0_guess + sigma_guess*2
        hist_fit, bins_fit = get_hist(ToF_MG[(ToF_MG >= fit_left) & (ToF_MG <= fit_right)],
                                      10, fit_left, fit_right)
        a, x0, sigma, __, __, pcov = fit_data(hist_fit, bins_fit, a_guess, x0_guess, sigma_guess)
        a_err, x0_err, sigma_err = np.sqrt(np.diag(pcov))
        FWHM_fit = 2 * np.sqrt(2*np.log(2)) * sigma
        FWHM_err = 2 * np.sqrt(2*np.log(2)) * sigma_err
        xx = np.linspace(fit_left, fit_right, 1000)
        norm = hist_MG[max_idx]/max(Gaussian(xx, a, x0, sigma))
        plt.plot(xx, Gaussian(xx, a, x0, sigma)*norm, color='black',
                 linestyle='-', label=None, zorder=50)
        # Store important values
        distance = voxel_to_distance_dict[MG_red.Bus, MG_red.gCh, MG_red.wCh][0]
        print('Layer: %d, distance: %')
        distances.append(distance)
        FWHMs.append(FWHM_fit)
        errors.append(FWHM_err)
        print('MG velocity [cm/µs], layer %d: %f' % (layer, (distance*100)/bin_centers_MG[max_idx]))
        # Plot He-3
        hist_He3, bins_He3 = np.histogram(ToF_He3, range=[start, stop], bins=number_bins)
        bin_centers_He3 = 0.5 * (bins_He3[1:] + bins_He3[:-1])
        # Calculate FWHM
        FWHM_He3, idx1, idx2, max_idx = calculate_FWHM(bin_centers_He3, hist_He3)
        # Fit data for He-3 instead
        middle = bin_centers_He3[max_idx]
        #hist_fit, bin_edges_fit = np.histogram(ToF_He3, range=[middle-15, middle+15], bins=50)
        #bins_fit =  0.5 * (bin_edges_fit[1:] + bin_edges_fit[:-1])
        a_guess, x0_guess, sigma_guess = get_fit_parameters_guesses(hist_He3, bin_centers_He3)
        fit_left, fit_right = x0_guess - sigma_guess*2, x0_guess + sigma_guess*2
        hist_fit, bins_fit = get_hist(ToF_He3[(ToF_He3 >= fit_left) & (ToF_He3 <= fit_right)],
                                      10, fit_left, fit_right)
        a, x0, sigma, __, __, pcov = fit_data(hist_fit, bins_fit, a_guess, x0_guess, sigma_guess)
        a_err, x0_err, sigma_err = np.sqrt(np.diag(pcov))
        FWHM_he3_fit = 2 * np.sqrt(2*np.log(2)) * sigma
        FWHM_he3_fit_err = 2 * np.sqrt(2*np.log(2)) * sigma_err
        # Plot Gaussian
        xx = np.linspace(fit_left, fit_right, 1000)
        norm = hist_He3[max_idx]/max(Gaussian(xx, a, x0, sigma))
        plt.plot(xx, Gaussian(xx, a, x0, sigma)*norm, color='black', linestyle='-', label=None, zorder=50)
        # Plot data
        plt.errorbar(bin_centers_He3, hist_He3, np.sqrt(hist_He3), fmt='.-',
                     capsize=5, zorder=5, label='He-3', color='red')
        background_He3_events = hist_He3[(bin_centers_He3 >= 16.8e3) & (bin_centers_He3 <= 16.9e3)]
        background_He3_per_bin = sum(background_He3_events)/len(background_He3_events)
        #plt.axhline(y=b, color='red', linewidth=2, label=None, linestyle='--')
        print('Distance in cm: %f' % (distance_He3*100))
        print('Flight time in us: %f' % bin_centers_He3[max_idx])
        print((distance_He3*100)/bin_centers_He3[max_idx])
        plt.legend(loc=2)
        plt.ylim(1, 1e5)
        plt.xlabel('ToF [µs]')
        plt.ylabel('Counts')
        plt.yscale('log')
        plt.title('Comparison ∆ToF, layer: %d' % layer)
        plt.grid(True, which='major', linestyle='--', zorder=0)
        plt.grid(True, which='minor', linestyle='--', zorder=0)
        # Save plot
        file_name = '%d.png' % layer
        output_path = output_folder + file_name
        fig.savefig(output_path, bbox_inches='tight')
        plt.close()
    # Animate
    images = []
    files = os.listdir(output_folder)
    files = [file[:-4] for file in files if file[-9:] != '.DS_Store' and file != '.gitignore']
    for file in sorted(files, key=int):
        images.append(imageio.imread(output_folder + file + '.png'))
    imageio.mimsave(animation_path, images)
    shutil.rmtree(output_folder, ignore_errors=True)
    # Plot FWHM
    fig = plt.figure()
    plt.xlabel('Distance [m]')
    plt.ylabel('FWHM [µs]')
    plt.title('FWHM ∆ToF in peaks as a function of distance')
    plt.errorbar(distances, FWHMs, errors, label='MG', fmt='.-', capsize=5, linestyle='', color='blue', zorder=5)
    plt.errorbar(distance_He3, FWHM_he3_fit, FWHM_he3_fit_err, label='He-3',
                 fmt='.-', capsize=5, color='red', linestyle='', zorder=5)
    # Fit MG data
    paras, pcov = np.polyfit(distances, FWHMs, 1, w=1/np.array(errors), cov=True)
    k, m = paras[0], paras[1]
    k_err, m_err = np.sqrt(np.diag(pcov))
    upper_k, upper_m = k - k_err, m + m_err
    lower_k, lower_m = k + k_err, m - m_err
    xx = np.linspace(28, 30, 100)
    plt.plot(xx, linear(xx, upper_k, upper_m), color='black', linestyle='-.', label='Multi-Grid Fit (Upper Bound)')
    plt.plot(xx, linear(xx, k, m), color='black', linestyle='-', label='Multi-Grid Fit')
    plt.plot(xx, linear(xx, lower_k, lower_m), color='black', linestyle='-.', label='Multi-Grid Fit (Lower Bound)')
    plt.grid(True, which='major', linestyle='--', zorder=0)
    plt.grid(True, which='minor', linestyle='--', zorder=0)
    plt.legend(loc=2)
    print('ToF spread: %f +/- %f us/m' % (k, k_err))
    fig.show()


# =============================================================================
#                      LAYERS INVESTIGATION - COUNTS
# =============================================================================

def investigate_layers_counts(df, duration):
    # Get count as a function of layer
    layers = np.arange(0, 20, 1)
    counts = [df[(df.wCh % 20) == layer].shape[0] for layer in layers]
    # Plot data
    fig = plt.figure()
    plt.grid(True, which='major', linestyle='--', zorder=0)
    plt.grid(True, which='minor', linestyle='--', zorder=0)
    plt.xlabel('Layer')
    plt.ylabel('Counts (Normalized by duration)')
    plt.title('Counts vs layer')
    plt.errorbar(layers, counts/duration, np.sqrt(counts)*(1/duration),
                 fmt='.-', capsize=5, zorder=5, color='black')
    fig.show()


# =============================================================================
#                      LAYERS INVESTIGATION - PHS
# =============================================================================

def investigate_layers_phs(df, duration):
    # Get count as a function of layer
    layers = np.arange(0, 20, 1)
    counts = [df[(df.wCh % 20) == layer].shape[0] for layer in layers]
    # Plot data
    fig = plt.figure()
    plt.grid(True, which='major', linestyle='--', zorder=0)
    plt.grid(True, which='minor', linestyle='--', zorder=0)
    plt.xlabel('Layer')
    plt.ylabel('Counts (Normalized by duration)')
    plt.title('Counts vs layer')
    plt.errorbar(layers, counts/duration, np.sqrt(counts)*(1/duration),
                 fmt='.-', capsize=5, zorder=5, color='black')
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
