#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Lineshape.py: Analyses the lineshape using our Figure-of-Merit
"""

import sys
import os
from Plotting.Analysis.DeltaE import energy_plot
from HelperFunctions.EnergyTransfer import calculate_energy
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.signal import find_peaks, peak_widths, peak_prominences
from scipy.optimize import curve_fit

from HelperFunctions.PathsAndFolders import mkdir_p

from HeliumTubes.EnergyHe3 import calculate_He3_energy
from HeliumTubes.PlottingHe3 import energy_plot_He3

from HelperFunctions.Filtering import filter_clusters
from HeliumTubes.FilteringHe3 import filter_He3

# =============================================================================
#                         LINESHAPE INVESTIGATION
# =============================================================================

def analyze_Lineshape(ce_MG, ce_He3, origin_voxel):
    """

    Non-Coated One Voxel: 800, 100
    Coated, Full Volume: 20000, 10000

    """
    def Gaussian(x, a, x0, sigma):
        return a*np.exp(-(x-x0)**2/(2*sigma**2))

    # Declare parameters
    number_bins = 5000
    plot_energy = True
    label_MG, label_He3 = 'Multi-Grid', 'He-3'
    # Plot Multi-Grid
    energies = calculate_energy(ce_MG, origin_voxel)
    energy_hist, bin_centers = energy_plot(ce_MG, origin_voxel,
                                           number_bins, 1, 10, plot_energy,
                                           label_MG)
    # Plot He-3
    He3_energies = calculate_He3_energy(ce_He3)
    hist_He3, bin_centers_He3 = energy_plot_He3(ce_He3, number_bins,
                                                plot_energy, label_He3)
    plt.legend()
    # Find peaks
    bins_part_1 = number_bins - number_bins//3
    bins_part_2 = number_bins//3
    heights_part1 = np.ones(bins_part_1)*20000
    heights_part2 = np.ones(bins_part_2)*10000
    heights = np.append(heights_part1, heights_part2)
    plt.plot(bin_centers[:bins_part_1], heights_part1, color='purple')
    plt.plot(bin_centers[bins_part_1:], heights_part2, color='purple')
    peaks, *_ = find_peaks(energy_hist, height=heights)
    widths, *_ = peak_widths(energy_hist, peaks)
    plt.plot(bin_centers[peaks], energy_hist[peaks], color='red', zorder=5,
             linestyle='', marker='o')
    figs = []
    titles = []
    print('Number of peaks: %d' % len(peaks))
    for width, peak in zip(widths, peaks):
        left, right = bin_centers[peak]-width/20, bin_centers[peak]+width/20
        #left_fit, right_fit = bin_centers[peak]-width/60, bin_centers[peak]+width/60
        # Perform new histogram to find correct peak limits
        peak_bins_temp = 100
        peak_energies_temp = energies[(energies >= left) & (energies <= right)]
        peak_hist_temp, peak_edges_temp = np.histogram(peak_energies_temp, bins=peak_bins_temp, range=[left, right])
        peak_bin_centers_temp = 0.5 * (peak_edges_temp[1:] + peak_edges_temp[:-1])
        # Prepare guesses
        maximum_temp = max(peak_hist_temp)
        maximum_idx = find_nearest(peak_hist_temp, maximum_temp)
        half_maximum = maximum_temp/2
        half_maximum_idx_1 = find_nearest(peak_hist_temp[:maximum_idx],
                                          half_maximum)
        half_maximum_idx_2 = (find_nearest(peak_hist_temp[maximum_idx:],
                                           half_maximum) + maximum_idx)
        FWHM = peak_bin_centers_temp[half_maximum_idx_2] - peak_bin_centers_temp[half_maximum_idx_1]
        a_guess = maximum_temp
        x0_guess = peak_bin_centers_temp[maximum_idx]
        sigma_guess = FWHM/(2*np.sqrt(2*np.log(2)))
        # Prepare new histograms within +/- 7 of our guessed sigma
        peak_bins = 100
        left_fit = x0_guess - (7 * sigma_guess)
        right_fit = x0_guess + (7 * sigma_guess)
        peak_energies = energies[(energies >= left_fit) & (energies <= right_fit)]
        peak_hist, peak_edges = np.histogram(peak_energies, bins=peak_bins, range=[left_fit, right_fit])
        peak_bin_centers = 0.5 * (peak_edges[1:] + peak_edges[:-1])
        maximum = max(peak_hist)
        # Define points to to fit procedure on
        fig_new = plt.figure()
        # Try fit
        try:
            # Fit
            popt, __ = curve_fit(Gaussian,
                                 peak_bin_centers,
                                 peak_hist,
                                 p0=[a_guess, x0_guess, sigma_guess])
            a, x0, sigma = popt[0], popt[1], abs(popt[2])
            # Plot sigma values
            plt.axvline(x=x0 - 5*sigma, color='orange', linewidth=2, label='-5σ')
            plt.axvline(x=x0 - 3*sigma, color='purple', linewidth=2, label='-3σ')
            plt.axvline(x=x0 - sigma, color='green', linewidth=2, label='-σ')
            plt.axvline(x=x0 + sigma, color='green', linewidth=2, label='σ')
            # Plot Gaussian
            xx = np.linspace(left_fit, right_fit, 1000)
            #plt.plot(xx, Gaussian(xx, a, x0, sigma)/a, color='green', label='Gaussian fit')
        except:
            print("Unexpected error:", sys.exc_info())
        # Define MG normalization
        MG_norm = 1/maximum
        # Print fit edges
        #plt.axvline(x=left_fit, color='green', linewidth=0.5)
        #plt.axvline(x=right_fit, color='green', linewidth=0.5)
        # Plot parameters from scipy to double-check
        #plt.axvline(x=left, color='orange', linewidth=0.5)
        #plt.axvline(x=right, color='orange', linewidth=0.5)
        plt.plot(peak_bin_centers, peak_hist*MG_norm, color='blue', marker='o',
                 linestyle='-', label=label_MG, zorder=5)
        #plt.plot(peak_bin_centers[fit_start:fit_stop],
        #         peak_hist[fit_start:fit_stop]*MG_norm,
        #         color='red', marker='x', linestyle='')
        plt.title('Peak at: %.2f meV (%.2f Å)' % (bin_centers[peak], meV_to_A(bin_centers[peak])))
        plt.xlabel('Energy [meV]')
        plt.ylabel('Counts (Normalized to maximum height)')
        #plt.yscale('log')
        # Plot He-3 data
        He3_peak_energies = He3_energies[(He3_energies >= left_fit) & (He3_energies <= right_fit)]
        He3_peak_hist, He3_peak_edges = np.histogram(He3_peak_energies,
                                                     bins=peak_bins, range=[left_fit, right_fit])
        He3_peak_bin_centers = 0.5 * (peak_edges[1:] + peak_edges[:-1])
        He3_norm = 1/max(He3_peak_hist)
        plt.plot(He3_peak_bin_centers, He3_peak_hist*He3_norm, color='red',
                 label=label_He3, marker='o', linestyle='-', zorder=5)
        plt.grid(True, which='major', linestyle='--', zorder=0)
        plt.grid(True, which='minor', linestyle='--', zorder=0)
        plt.legend()
        # Append figure
        figs.append(fig_new)
        titles.append('peak_at_%.2f_meV.pdf' % bin_centers[peak])
        plt.close()


    # Save all figures
    dirname = os.path.dirname(__file__)
    output_folder = os.path.join(dirname, '../../../Output/Lineshape/')
    for fig_temp, title in zip(figs, titles):
        output_path = output_folder + title
        fig_temp.savefig(output_path, bbox_inches='tight')



def analyze_all_lineshapes(origin_voxel, MG_filter_parameters, He3_filter_parameters):
    # Define parameters
    colors = {'MG_Coated': 'blue', 'MG_Non_Coated': 'green', 'He3': 'red'}
    # Prepare data
    MG_coated_data, MG_non_coated_data, He3_data = prepare_data(origin_voxel, MG_filter_parameters, He3_filter_parameters)
    # Plot all individual peaks
    plot_all_peaks(MG_coated_data, 'MG_Coated', colors['MG_Coated'])
    plot_all_peaks(MG_non_coated_data, 'MG_Non_Coated', colors['MG_Non_Coated'])
    plot_all_peaks(He3_data, 'He3', colors['He3'])
    # Plot MG compared to He-3
    plot_all_peaks_from_two_data_sets(MG_coated_data, 'MG_Coated', colors['MG_Coated'],
                                      He3_data, 'He3', colors['He3'])
    plot_all_peaks_from_two_data_sets(MG_non_coated_data, 'MG_Non_Coated', colors['MG_Non_Coated'],
                                      He3_data, 'He3', colors['He3'])
    # Plot Coated Radial blades compared to non-coated radial blades
    plot_all_peaks_from_two_data_sets(MG_non_coated_data, 'MG_Non_Coated', colors['MG_Non_Coated'],
                                      MG_coated_data, 'MG_Coated', colors['MG_Coated'])


def plot_all_peaks_from_two_data_sets(data_1, label_1, color_1, data_2, label_2, color_2):
    # Prepare output paths
    dirname = os.path.dirname(__file__)
    output_folder = os.path.join(dirname, '../../../Output/Comparison_%s_and_%s/' % (label_1, label_2))
    mkdir_p(output_folder)
    # Extract parameters
    energies, peaks, widths, hist, bins = data_1[0], data_1[1], data_1[2], data_1[3], data_1[4]
    energies_2 = data_2[0]
    number_bins = 100
    # Iterate through all peaks
    for width, peak in zip(widths, peaks):
        # Extract fit guesses and peak borders
        left, right = bins[peak]-width/20, bins[peak]+width/20
        hist_peak, bins_peak = get_hist(energies, number_bins, left, right)
        a_guess, x0_guess, sigma_guess = get_fit_parameters_guesses(hist_peak, bins_peak)
        # Prepare peak within +/- 7 of our estimated sigma
        left_fit, right_fit = (x0_guess - (7 * sigma_guess)), (x0_guess + (7 * sigma_guess))
        hist_fit, bins_fit = get_hist(energies, number_bins, left_fit, right_fit)
        try:
            # Fit data
            a, x0, sigma = fit_data(hist_fit, bins_fit, a_guess, x0_guess, sigma_guess)
            fig = plt.figure()
            plot_sigma_borders(x0, sigma)
            # Plot from main data
            norm_1 = 1/max(hist_fit)
            plt.plot(bins_fit, hist_fit*norm_1, marker='o', linestyle='-', label=label_1,
                     zorder=5, color=color_1)
            # Plot from second data
            hist_2, bins_2 = get_hist(energies_2, number_bins, left_fit, right_fit)
            norm_2 = 1/max(hist_2)
            plt.plot(bins_2, hist_2*norm_2, marker='o', linestyle='-', label=label_2,
                     zorder=5, color=color_2)
            # Stylise plot
            plt.grid(True, which='major', linestyle='--', zorder=0)
            plt.grid(True, which='minor', linestyle='--', zorder=0)
            plt.title('Peak at: %.2f meV (%.2f Å)' % (bins[peak], meV_to_A(bins[peak])))
            plt.xlabel('Energy [meV]')
            plt.ylabel('Counts (Normalized to maximum)')
            plt.xlim(x0 - (7 * sigma), x0 + (7 * sigma))
            plt.yscale('log')
            plt.legend(loc=1)
            # Save plot
            file_name = '%s_Peak_at_%.2f_meV_(%.2f_Å).pdf' % (label_1, bins[peak], meV_to_A(bins[peak]))
            output_path = output_folder + file_name
            fig.savefig(output_path, bbox_inches='tight')
            plt.close()
        except:
            print("Unexpected error:", sys.exc_info())



def plot_all_peaks(data, label, color):
    # Prepare output paths
    dirname = os.path.dirname(__file__)
    output_folder = os.path.join(dirname, '../../../Output/%s/' % label)
    mkdir_p(output_folder)
    # Extract parameters
    energies, peaks, widths, hist, bins = data[0], data[1], data[2], data[3], data[4]
    number_bins = 100
    # Iterate through all peaks
    for width, peak in zip(widths, peaks):
        # Extract fit guesses and peak borders
        left, right = bins[peak]-width/20, bins[peak]+width/20
        hist_peak, bins_peak = get_hist(energies, number_bins, left, right)
        a_guess, x0_guess, sigma_guess = get_fit_parameters_guesses(hist_peak, bins_peak)
        # Prepare peak within +/- 7 of our estimated sigma
        left_fit, right_fit = (x0_guess - (7 * sigma_guess)), (x0_guess + (7 * sigma_guess))
        hist_fit, bins_fit = get_hist(energies, number_bins, left_fit, right_fit)
        try:
            fig = plt.figure()
            # Fit data
            a, x0, sigma = fit_data(hist_fit, bins_fit, a_guess, x0_guess, sigma_guess)
            plot_sigma_borders(x0, sigma)
            plt.plot(bins_fit, hist_fit, marker='o', linestyle='-', label=label,
                     zorder=5, color=color)
            # Stylise plot
            plt.grid(True, which='major', linestyle='--', zorder=0)
            plt.grid(True, which='minor', linestyle='--', zorder=0)
            plt.title('Peak at: %.2f meV (%.2f Å)' % (bins[peak], meV_to_A(bins[peak])))
            plt.xlabel('Energy [meV]')
            plt.ylabel('Counts')
            plt.xlim(x0 - (7 * sigma), x0 + (7 * sigma))
            plt.legend(loc=1)
            # Save plot
            file_name = '%s_Peak_at_%.2f_meV_(%.2f_Å).pdf' % (label, bins[peak], meV_to_A(bins[peak]))
            output_path = output_folder + file_name
            fig.savefig(output_path, bbox_inches='tight')
            plt.close()
        except:
            print("Unexpected error:", sys.exc_info())


def fit_data(hist, bins, a_guess, x0_guess, sigma_guess):
    def Gaussian(x, a, x0, sigma):
        return a*np.exp(-(x-x0)**2/(2*sigma**2))

    popt, __ = curve_fit(Gaussian, bins, hist, p0=[a_guess, x0_guess, sigma_guess])
    a, x0, sigma = popt[0], popt[1], abs(popt[2])
    # Plot Gaussian
    xx = np.linspace(bins[0], bins[-1], 1000)
    plt.plot(xx, Gaussian(xx, a, x0, sigma), label='Gaussian fit', color='black')
    return a, x0, sigma


def plot_sigma_borders(x0, sigma):
    plt.axvline(x=x0 - 5*sigma, color='orange', linewidth=2, label='-5σ')
    plt.axvline(x=x0 - 3*sigma, color='purple', linewidth=2, label='-3σ')
    plt.axvline(x=x0 - sigma, color='green', linewidth=2, label='-σ')
    plt.axvline(x=x0 + sigma, color='green', linewidth=2, label='σ')





def get_fit_parameters_guesses(hist, bins):
    # Extract relavant parameters
    maximum = max(hist)
    maximum_idx = find_nearest(hist, maximum)
    half_maximum = maximum/2
    half_maximum_idx_1 = find_nearest(hist[:maximum_idx], half_maximum)
    half_maximum_idx_2 = find_nearest(hist[maximum_idx:], half_maximum) + maximum_idx
    FWHM = bins[half_maximum_idx_2] - bins[half_maximum_idx_1]
    # Calculate guesses
    a_guess = maximum
    x0_guess = bins[maximum_idx]
    sigma_guess = FWHM/(2*np.sqrt(2*np.log(2)))
    return a_guess, x0_guess, sigma_guess




def prepare_data(origin_voxel, MG_filter_parameters, He3_filter_parameters):
    # Declare parameters
    dirname = os.path.dirname(__file__)
    number_bins = 5000
    start = 0.8  # [meV]
    end = 80  # [meV]
    heights_MG_coated = [20000, 10000]
    heights_MG_non_coated = [12000, 1000]
    heights_He3 = [20000, 1000]
    # Declare file names
    MG_COATED = 'mvmelst_165_191002_111641_Det2_overnight3.h5'
    MG_NON_COATED = 'mvmelst_135_190930_141618_Det1_overnight2_30x80_14x60.h5'
    HE_3 = '2019_09_HZB_He3InBeam54304s_overnight.h5'
    # Declare paths to data
    MG_COATED_RADIAL_PATH = os.path.join(dirname, '../../../Data/Lineshape/%s' % MG_COATED)
    MG_NON_COATED_RADIAL_PATH = os.path.join(dirname, '../../../Data/Lineshape/%s' % MG_NON_COATED)
    HE_3_PATH = os.path.join(dirname, '../../../Data/Lineshape/%s' % HE_3)
    # Import data
    print('Importing data...')
    df_MG_coated = pd.read_hdf(MG_COATED_RADIAL_PATH, 'ce')
    df_MG_non_coated = pd.read_hdf(MG_NON_COATED_RADIAL_PATH, 'ce')
    df_He3 = pd.read_hdf(HE_3_PATH, 'df')
    # Filter data
    print('Filtering data...')
    df_MG_coated_filtered = filter_clusters(df_MG_coated, MG_filter_parameters)
    df_MG_non_coated_filtered = filter_clusters(df_MG_non_coated, MG_filter_parameters)
    df_He3_filtered = filter_He3(df_He3, He3_filter_parameters)
    # Extract energies
    print('Extracting energies...')
    energies_MG_coated = calculate_energy(df_MG_coated_filtered, origin_voxel)
    energies_MG_non_coated = calculate_energy(df_MG_non_coated_filtered, origin_voxel)
    energies_He3 = calculate_He3_energy(df_He3_filtered)
    # Histogram data
    print('Histograming data...')
    hist_MG_coated, bins_MG_coated = get_hist(energies_MG_coated, number_bins, start, end)
    hist_MG_non_coated, bins_MG_non_coated = get_hist(energies_MG_non_coated, number_bins, start, end)
    hist_He3, bins_He3 = get_hist(energies_He3, number_bins, start, end)
    # Get peaks
    print('Extracting peaks...')
    peaks_MG_coated = get_peaks(hist_MG_coated, heights_MG_coated, number_bins)
    peaks_MG_non_coated = get_peaks(hist_MG_non_coated, heights_MG_non_coated, number_bins)
    peaks_He3 = get_peaks(hist_He3, heights_He3, number_bins)
    # Get widths
    print('Calculating widths...')
    widths_MG_coated, *_ = peak_widths(hist_MG_coated, peaks_MG_coated)
    widths_MG_non_coated, *_ = peak_widths(hist_MG_non_coated, peaks_MG_non_coated)
    widths_He3, *_ = peak_widths(hist_He3, peaks_He3)
    # Store data from each detector in a separate list
    MG_coated_data = [energies_MG_coated,
                      peaks_MG_coated,
                      widths_MG_coated,
                      hist_MG_coated,
                      bins_MG_coated]
    MG_non_coated_data = [energies_MG_non_coated,
                          peaks_MG_non_coated,
                          widths_MG_non_coated,
                          hist_MG_non_coated,
                          bins_MG_non_coated]
    He3_data = [energies_He3,
                peaks_He3,
                widths_He3,
                hist_He3,
                bins_He3]
    return MG_coated_data, MG_non_coated_data, He3_data


# =============================================================================
#                                HELPER FUNCTIONS
# =============================================================================

def get_peaks(hist, heights, number_bins):
    # Extract heights
    height_1 = heights[0]
    height_2 = heights[1]
    # Define vectors with heights, used to find peaks above height values
    bins_part_1 = number_bins - number_bins//3
    bins_part_2 = number_bins//3
    heights_part1 = np.ones(bins_part_1)*height_1
    heights_part2 = np.ones(bins_part_2)*height_2
    heights = np.append(heights_part1, heights_part2)
    # Histogram data
    start = A_to_meV(10)
    stop = A_to_meV(1)
    # Get peaks
    peaks, *_ = find_peaks(hist, height=heights)
    return peaks


def get_hist(energies, number_bins, start, stop):
    hist, bin_edges = np.histogram(energies, bins=number_bins, range=[start, stop])
    bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
    return hist, bin_centers


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


def meV_to_A(energy):
    return np.sqrt(81.81/energy)


def A_to_meV(wavelength):
    return (81.81/(wavelength ** 2))
