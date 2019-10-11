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


def analyze_all_lineshapes(origin_voxel, MG_filter_parameters, He3_filter_parameters):
    # Define parameters
    colors = {'MG_Coated': 'blue', 'MG_Non_Coated': 'green', 'He3': 'red'}
    # Prepare data
    full_data = prepare_data(origin_voxel, MG_filter_parameters, He3_filter_parameters)
    MG_coated_data, MG_non_coated_data, He3_data = full_data[0], full_data[1], full_data[4]
    MG_coated_background, MG_non_coated_background, He3_background = full_data[2], full_data[3], full_data[5]
    # Plot all individual peaks
    #plot_all_peaks(MG_coated_data, 'MG_Coated', colors['MG_Coated'])
    #plot_all_peaks(MG_non_coated_data, 'MG_Non_Coated', colors['MG_Non_Coated'])
    #plot_all_peaks(He3_data, 'He3', colors['He3'])
    # Plot MG compared to He-3
    plot_all_peaks_from_two_data_sets(MG_coated_data, 'MG_Coated', colors['MG_Coated'], He3_data, 'He3', colors['He3'])
    plot_all_peaks_from_two_data_sets(MG_non_coated_data, 'MG_Non_Coated', colors['MG_Non_Coated'], He3_data, 'He3', colors['He3'])
    # Plot Coated Radial blades compared to non-coated radial blades
    plot_all_peaks_from_two_data_sets(MG_non_coated_data, 'MG_Non_Coated', colors['MG_Non_Coated'], MG_coated_data, 'MG_Coated', colors['MG_Coated'])
    # Plot data with background and extract important values
    FWHM_Coated, FoM_Coated, err_Coated, energies_Coated = plot_all_peaks_with_background(MG_coated_data, 'MG_Coated', colors['MG_Coated'], MG_coated_background)
    FWHM_NonCoated, FoM_NonCoated, err_NonCoated, energies_NonCoated = plot_all_peaks_with_background(MG_non_coated_data, 'MG_Non_Coated', colors['MG_Non_Coated'], MG_non_coated_background)
    FWHM_He3, FoM_He3, err_He3, energies_He3 = plot_all_peaks_with_background(He3_data, 'He3', colors['He3'], He3_background)
    # Plot important values
    FWHMs = [FWHM_Coated, FWHM_NonCoated, FWHM_He3]
    energies = [energies_Coated, energies_NonCoated, energies_He3]
    labels = ['MG_Coated', 'MG_Non_Coated', 'He3']
    FoMs = [FoM_Coated, FoM_NonCoated, FoM_He3]
    errors = [err_Coated, err_NonCoated, err_He3]
    # Plot FWHM
    fig = plt.figure()
    for energy, FWHM, label in zip(energies, FWHMs, labels):
        plot_FWHM(energy, FWHM, label)
    plt.legend()
    fig.show()
    # Plot FoM
    fig = plt.figure()
    for energy, FoM, error, label in zip(energies, FoMs, errors, labels):
        plot_FoM(energy, FoM, error, label)
    plt.legend()
    fig.show()

def plot_FoM(energies, FoMs, errors, label):
    plt.title('Figure-of-Merit')
    plt.ylabel('FoM [shoulder/peak]')
    plt.xlabel('Peak energy [meV]')
    plt.grid(True, which='major', linestyle='--', zorder=0)
    plt.grid(True, which='minor', linestyle='--', zorder=0)
    plt.xscale('log')
    plt.errorbar(energies, FoMs, errors, fmt='.-', capsize=5, zorder=5, label=label)


def plot_FWHM(energies, FWHMs, label):
    plt.title('FWHM')
    plt.ylabel('FWHM [meV]')
    plt.xlabel('Peak energy [meV]')
    plt.grid(True, which='major', linestyle='--', zorder=0)
    plt.grid(True, which='minor', linestyle='--', zorder=0)
    plt.xscale('log')
    plt.yscale('log')
    plt.plot(energies, FWHMs, label=label, marker='o', linestyle='-', zorder=5)


def plot_all_peaks_with_background(data, label_data, color_data, background):
    # Prepare output paths
    dirname = os.path.dirname(__file__)
    output_folder = os.path.join(dirname, '../../../Output/%s_with_background/' % label_data)
    mkdir_p(output_folder)
    # Extract parameters
    duration, energies, hist, bins, peaks, widths = data[0], data[1], data[2], data[3], data[4], data[5]
    duration_background, energies_background = background[0], background[1]
    number_bins = 100
    # Define lists to store important data
    FWHM_list = []
    FoM_list = []
    uncertainites = []
    peak_energies = []
    # Iterate through all peaks
    for width, peak in zip(widths, peaks):
        # Extract fit guesses and peak borders
        left, right = bins[peak]-width/20, bins[peak]+width/20
        hist_peak, bins_peak = get_hist(energies, number_bins, left, right)
        a_guess, x0_guess, sigma_guess = get_fit_parameters_guesses(hist_peak, bins_peak)
        # Prepare peak within +/- 5 of our estimated sigma, we'll use this to fit
        left_fit, right_fit = (x0_guess - (5 * sigma_guess)), (x0_guess + (5 * sigma_guess))
        hist_fit, bins_fit = get_hist(energies, number_bins, left_fit, right_fit)
        try:
            # Fit data
            a, x0, sigma = fit_data(hist_fit, bins_fit, a_guess, x0_guess, sigma_guess)
            fig = plt.figure()
            plot_sigma_borders(x0, sigma)
            # Prepare data within +/- 25 of our estimated sigma, we'll use this to plot
            left_plot, right_plot = (x0 - (25 * sigma)), (x0+ (25 * sigma))
            hist_plot, bins_plot = get_hist(energies, number_bins, left_plot, right_plot)
            # Plot from beam
            norm_beam = 1/duration
            plt.plot(bins_plot, hist_plot*norm_beam, marker='o', linestyle='-', label=label_data,
                     zorder=5, color=color_data)
            # Plot from background
            hist_background, bins_background = get_hist(energies_background, number_bins, left_plot, right_plot)
            norm_background = 1/duration_background
            plt.plot(bins_background, hist_background*norm_background, marker='o',
                     linestyle='-', label='Background', zorder=5, color='black')
            # Stylise plot
            plt.grid(True, which='major', linestyle='--', zorder=0)
            plt.grid(True, which='minor', linestyle='--', zorder=0)
            plt.title('Peak at: %.2f meV (%.2f Å)' % (bins[peak], meV_to_A(bins[peak])))
            plt.xlabel('Energy [meV]')
            plt.ylabel('Counts (Normalized to duration)')
            plt.xlim(x0 - (25 * sigma), x0 + (25 * sigma))
            plt.yscale('log')
            plt.legend(loc=1)
            # Save plot
            file_name = '%s_Peak_at_%.2f_meV_(%.2f_Å).pdf' % (label_data, bins[peak], meV_to_A(bins[peak]))
            output_path = output_folder + file_name
            fig.savefig(output_path, bbox_inches='tight')
            plt.close()
            # Store important values
            FWHM = 2 * np.sqrt(2*np.log(2)) * sigma
            FoM, uncertainity = get_FoM(energies, energies_background, x0, sigma, duration, duration_background)
            FWHM_list.append(FWHM)
            FoM_list.append(FoM)
            uncertainites.append(uncertainity)
            peak_energies.append(bins[peak])
        except:
            print("Unexpected error:", sys.exc_info())
    return FWHM_list, FoM_list, uncertainites, peak_energies


def plot_all_peaks_from_two_data_sets(data_1, label_1, color_1, data_2, label_2, color_2):
    # Prepare output paths
    dirname = os.path.dirname(__file__)
    output_folder = os.path.join(dirname, '../../../Output/Comparison_%s_and_%s/' % (label_1, label_2))
    mkdir_p(output_folder)
    # Extract parameters
    energies, hist, bins, peaks, widths = data_1[1], data_1[2], data_1[3], data_1[4], data_1[5]
    energies_2 = data_2[1]
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
            # Prepare data within +/- 25 of our estimated sigma, we'll use this to plot
            left_plot, right_plot = (x0 - (25 * sigma)), (x0+ (25 * sigma))
            hist_plot, bins_plot = get_hist(energies, number_bins, left_plot, right_plot)
            # Plot from main data
            norm_1 = 1/max(hist_plot)
            plt.plot(bins_plot, hist_plot*norm_1, marker='o', linestyle='-', label=label_1,
                     zorder=5, color=color_1)
            # Plot from second data
            hist_2, bins_2 = get_hist(energies_2, number_bins, left_plot, right_plot)
            norm_2 = 1/max(hist_2)
            plt.plot(bins_2, hist_2*norm_2, marker='o', linestyle='-', label=label_2,
                     zorder=5, color=color_2)
            # Stylise plot
            plt.grid(True, which='major', linestyle='--', zorder=0)
            plt.grid(True, which='minor', linestyle='--', zorder=0)
            plt.title('Peak at: %.2f meV (%.2f Å)' % (bins[peak], meV_to_A(bins[peak])))
            plt.xlabel('Energy [meV]')
            plt.ylabel('Counts (Normalized to maximum)')
            plt.xlim(x0 - (25 * sigma), x0 + (25 * sigma))
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
    energies, hist, bins, peaks, widths = data[1], data[2], data[3], data[4], data[5]
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
    """
    Data is returned in following order:

    1. Multi-Grid Coated Radial Blades (beam)
    2. Multi-Grid Non-Coated Radial Blades (beam)
    3. Multi-Grid Coated Radial Blades (background)
    4. Multi-Grid Non-Coated Radial Blades (background)
    5. He-3 (beam)
    6. He-3 (background)

    Within each data-list, data is returned in following order

    1. Duration
    2. Energies
    3. Histogram
    4. Bin centers
   (5. Peaks)
   (6. Widths)

    """
    # Declare parameters
    dirname = os.path.dirname(__file__)
    number_bins = 5000
    start = 0.8  # [meV]
    end = 80  # [meV]
    heights_MG_coated = [20000, 10000]
    heights_MG_non_coated = [12000, 1000]
    heights_He3 = [20000, 1000]
    heights_vec_MG = [heights_MG_coated, heights_MG_non_coated]
    # Declare file names
    MG_COATED = 'mvmelst_165_191002_111641_Det2_overnight3.h5'
    MG_COATED_BACKGROUND = 'mvmelst_169_191003_075039_Det2_He3InBeam_overnight4.h5'
    MG_NON_COATED = 'mvmelst_135_190930_141618_Det1_overnight2_30x80_14x60.h5'
    MG_NON_COATED_BACKGROUND = 'mvmelst_141_191001_120405_He3InBeam_overnight3.h5'
    HE_3 = '2019_09_HZB_He3InBeam54304s_overnight.h5'
    HE_3_BACKGROUND = '2019_09_HZB_out_of_beam_overnight_58094s.h5'
    MG_file_names = [MG_COATED, MG_NON_COATED,
                     MG_COATED_BACKGROUND, MG_NON_COATED_BACKGROUND]
    He3_file_names = [HE_3, HE_3_BACKGROUND]
    # Declare list to store all data
    full_data = []
    # Store Multi-Grid data
    print('Multi-Grid...')
    for i, file_name in enumerate(MG_file_names):
        path = os.path.join(dirname, '../../../Data/Lineshape/%s' % file_name)
        df = pd.read_hdf(path, 'ce')
        df_red = filter_clusters(df, MG_filter_parameters)
        duration = get_duration(df)
        energies = calculate_energy(df_red, origin_voxel)
        hist, bins = get_hist(energies, number_bins, start, end)
        data = [duration, energies, hist, bins]
        if i < 2:
            # If it is a beam measurement, extract peaks
            peaks = get_peaks(hist, heights_vec_MG[i], number_bins)
            widths, *_ = peak_widths(hist, peaks)
            data.extend([peaks, widths])
        full_data.append(data)
    # Store He-3 data
    He3_durations = [54304, 58094]
    print('He-3...')
    for i, (file_name, duration) in enumerate(zip(He3_file_names, He3_durations)):
        path = os.path.join(dirname, '../../../Data/Lineshape/%s' % file_name)
        df = pd.read_hdf(path, 'df')
        df_red = filter_He3(df, He3_filter_parameters)
        energies = calculate_He3_energy(df_red)
        hist, bins = get_hist(energies, number_bins, start, end)
        data = [duration, energies, hist, bins]
        if i < 1:
            # If it is a beam measurement, extract peaks
            peaks = get_peaks(hist, heights_He3, number_bins)
            widths, *_ = peak_widths(hist, peaks)
            data.extend([peaks, widths])
        full_data.append(data)
    return full_data


# =============================================================================
#                                HELPER FUNCTIONS
# =============================================================================

def get_FoM(beam, background, x0, sigma, duration_beam, duration_background):
    # Extract number of counts from regions of interest
    peak_beam_counts = beam[(beam >= (x0 - 3*sigma)) & (beam <= (x0 + 3*sigma))]
    peak_background_counts = background[(background >= (x0 - 3*sigma)) & (background <= (x0 + 3*sigma))]
    shoulder_beam_counts = beam[(beam >= (x0 - 5*sigma)) & (beam <= (x0 - 3*sigma))]
    shoulder_background_counts = background[(background >= (x0 - 5*sigma)) & (background <= (x0 - 3*sigma))]
    # Rename for easier calculation of uncertainties
    a = len(peak_beam_counts)
    b = len(peak_background_counts)
    c = len(shoulder_beam_counts)
    d = len(shoulder_background_counts)
    # Subtract background
    norm = (duration_beam/duration_background)
    e = a - b * norm
    f = c - d * norm
    # Calculate FoM
    g = f/e
    # Calculate uncertainites
    da = np.sqrt(a)
    db = np.sqrt(b)
    dc = np.sqrt(c)
    dd = np.sqrt(d)
    de = np.sqrt(da ** 2 + (db*norm) ** 2)
    df = np.sqrt(dc ** 2 + (dd*norm) ** 2)
    dg = np.sqrt((de/e) ** 2 + (df/f) ** 2)
    uncertainty = dg * g
    FoM = g
    return FoM, uncertainty


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

def get_duration(df):
    times = df.Time.values
    duration_in_seconds = (times[-1] - times[0]) * 62.5e-9
    return duration_in_seconds
