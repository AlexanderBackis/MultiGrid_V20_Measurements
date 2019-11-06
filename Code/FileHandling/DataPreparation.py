#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
DataPreparation.py:
"""
import os
import pandas as pd
import numpy as np
from scipy.signal import peak_widths
import matplotlib.pyplot as plt

from HelperFunctions.Filtering import filter_clusters
from HelperFunctions.EnergyCalculation import calculate_energy
from HelperFunctions.Fitting import get_hist, get_fit_parameters_guesses, fit_data
from HelperFunctions.Misc import get_duration, mkdir_p, meV_to_A
from HelperFunctions.PeakFinding import get_peaks

from Plotting.Analysis.Lineshape import get_FoM, calculate_distance_borders
from Plotting.Analysis.Efficiency import get_peak_area

from HeliumTubes.FilteringHe3 import filter_He3
from HeliumTubes.EnergyHe3 import calculate_He3_energy


# =============================================================================
#                               PREPARE DATA
# =============================================================================

def prepare_data(origin_voxel, MG_filter_parameters, He3_filter_parameters):
    """
    Data is returned in following order:

    1. Multi-Grid Coated Radial Blades, beam
    2. Multi-Grid Non-Coated Radial Blades, beam
    3. Multi-Grid Coated Radial Blades, background
    4. Multi-Grid Non-Coated Radial Blades, background
    5. He-3, beam
    6. He-3, background

    Within each data-list, data is returned in following order

    1. Duration
    2. Energies
    3. Histogram
    4. Bin centers
   (5. Peaks)
   (6. Widths)

    """
    # Declare parameters, such as distance offset and He3 duration
    dirname = os.path.dirname(__file__)
    number_bins = 5000
    start = 0.8  # [meV]
    end = 80  # [meV]
    MG_distance_offsets = [1.5e-3, 0, 1.5e-3, 0]
    He3_distance_offset = 3e-3
    He3_durations = [54304, 58094]
    # Declare heights used as threshold in peak finding algorithm
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
    MG_file_names = [MG_COATED, MG_NON_COATED, MG_COATED_BACKGROUND, MG_NON_COATED_BACKGROUND]
    He3_file_names = [HE_3, HE_3_BACKGROUND]
    # Declare list to store all data
    full_data = []
    # Store Multi-Grid data
    print('Multi-Grid...')
    for i, file_name in enumerate(MG_file_names):
        path = os.path.join(dirname, '../../Data/Lineshape/%s' % file_name)
        df = pd.read_hdf(path, 'ce')
        df_red = filter_clusters(df, MG_filter_parameters)
        duration = get_duration(df)
        energies = calculate_energy(df_red, origin_voxel, MG_distance_offsets[i])
        hist, bins = get_hist(energies, number_bins, start, end)
        data = [duration, energies, hist, bins]
        if i < 2:
            # If it is a beam measurement, extract peaks
            peaks = get_peaks(hist, heights_vec_MG[i], number_bins)
            widths, *_ = peak_widths(hist, peaks)
            data.extend([peaks, widths])
        full_data.append(data)
    # Store He-3 data
    print('He-3...')
    for i, (file_name, duration) in enumerate(zip(He3_file_names, He3_durations)):
        path = os.path.join(dirname, '../../Data/Lineshape/%s' % file_name)
        df = pd.read_hdf(path, 'df')
        df_red = filter_He3(df, He3_filter_parameters)
        energies = calculate_He3_energy(df_red, He3_distance_offset)
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
#                          EXTRACT KEY PARAMETERS
# =============================================================================

def plot_all_peaks(data, label, color, chopper_to_detector_distance):
    # Prepare output paths
    dirname = os.path.dirname(__file__)
    output_folder = os.path.join(dirname, '../../Output/%s/' % label)
    mkdir_p(output_folder)
    # Extract parameters
    energies, hist, bins, peaks, widths = data[1], data[2], data[3], data[4], data[5]
    number_bins = 100
    # Declarea vectors to store data
    peak_energies = []
    FoMs = []
    FoM_uncertainites = []
    peak_areas = []
    peak_area_uncertainties = []
    # Declare peak locations and widths
    peak_values = bins[peaks]
    width_values = widths
    # Add last two peaks
    if label == 'He3':
        peak_values = np.append(peak_values, [77, 104])
        width_values = np.append(width_values, [widths[-1]/2, widths[-1]/2])
    elif label == 'MG_Non_Coated':
        peak_values = np.append(peak_values, [77.5, 105])
        width_values = np.append(width_values, [widths[-1]/2, widths[-1]/2])
    # Iterate through all peaks
    for width, peak in zip(width_values, peak_values):
        # Extract fit guesses and peak borders
        left, right = peak-width/20, peak+width/20
        hist_peak, bins_peak = get_hist(energies, number_bins, left, right)
        a_guess, x0_guess, sigma_guess = get_fit_parameters_guesses(hist_peak, bins_peak)
        # Prepare peak within +/- 7 of our estimated sigma
        if peak < 70:
            left_fit, right_fit = (x0_guess - (7 * sigma_guess)), (x0_guess + (7 * sigma_guess))
        else:
            left_fit, right_fit = (x0_guess - (4 * sigma_guess)), (x0_guess + (4 * sigma_guess))
        hist_fit, bins_fit = get_hist(energies, number_bins, left_fit, right_fit)
        fig = plt.figure()
        # Fit data
        a, x0, sigma, x_fit, y_fit, *_ = fit_data(hist_fit, bins_fit, a_guess, x0_guess, sigma_guess)
        # Plot data
        if peak < 70:
            left_plot, right_plot = (x0 - (31 * sigma)), (x0 + (7 * sigma))
            plt.xlim(x0 - (31 * sigma), x0 + (7 * sigma))
        else:
            left_plot, right_plot = (x0 - (10 * sigma)), (x0 + (10 * sigma))
            plt.xlim(x0 - (10 * sigma), x0 + (10 * sigma))
        hist_plot, bins_plot = get_hist(energies, number_bins, left_plot, right_plot)
        plt.errorbar(bins_plot, hist_plot, np.sqrt(hist_plot), fmt='.-', capsize=5, zorder=5, label=label, color=color)
        plt.plot(x_fit, y_fit*(max(hist_plot)/max(y_fit)), label='Gaussian fit', color='black')
        # Plot where the shoulder should be
        reduced_energies, distances, linestyles = calculate_distance_borders(bins_plot, hist_plot, chopper_to_detector_distance)
        for distance, linestyle in zip(distances[1:3], linestyles[1:3]):
            E_new = reduced_energies[distance]
            plt.axvline(x=E_new, linewidth=2, zorder=10, color='black', linestyle=linestyle,
                        label='Extra distance: %d cm' % (distance*100))
        # Extract FoM
        start, end = reduced_energies[0.20], reduced_energies[0.10]
        bin_width = bins_plot[1] - bins_plot[0]
        FoM, FoM_uncertainity = get_FoM(energies, x0, sigma, start, end, bin_width)
        # Extract area
        peak_area, peak_area_uncertainity = get_peak_area(energies, x0, sigma, bin_width)
        # Save all important values
        peak_energies.append(peak)
        FoMs.append(FoM)
        FoM_uncertainites.append(FoM_uncertainity)
        peak_areas.append(peak_area)
        peak_area_uncertainties.append(peak_area_uncertainity)
        # Stylise plot
        plt.grid(True, which='major', linestyle='--', zorder=0)
        plt.grid(True, which='minor', linestyle='--', zorder=0)
        plt.title('Peak at: %.2f meV (%.2f Å)' % (peak, meV_to_A(peak)))
        plt.xlabel('Energy [meV]')
        plt.ylabel('Counts')
        plt.yscale('log')
        plt.legend(loc=2)
        # Save plot
        file_name = '%s_Peak_at_%.2f_meV_(%.2f_Å).pdf' % (label, peak, meV_to_A(peak))
        output_path = output_folder + file_name
        fig.savefig(output_path, bbox_inches='tight')
        plt.close()
    return peak_energies, FoMs, FoM_uncertainites, peak_areas, peak_area_uncertainties



# =============================================================================
#                   PREPARE COMPARISON BETWEEN ALL DETECTORS
# =============================================================================

def plot_all_peaks_from_three_data_sets(data_1, label_1, color_1,
                                        data_2, label_2, color_2,
                                        data_3, label_3, color_3):
    # Prepare output paths
    dirname = os.path.dirname(__file__)
    output_folder = os.path.join(dirname, '../../Output/Comparison_%s_and_%s/' % (label_1, label_2))
    mkdir_p(output_folder)
    # Extract parameters
    energies, hist, bins, peaks, widths = data_1[1], data_1[2], data_1[3], data_1[4], data_1[5]
    energies_2 = data_2[1]
    energies_3 = data_3[1]
    number_bins = 100
    # Iterate through all peaks
    FoM_list = []
    energies_list = []
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
            a, x0, sigma, *_ = fit_data(hist_fit, bins_fit, a_guess, x0_guess, sigma_guess)
            fig = plt.figure()
            #plot_sigma_borders(x0, sigma)
            # Prepare data within +/- 25 of our estimated sigma, we'll use this to plot
            left_plot, right_plot = (x0 - (50 * sigma)), (x0 + (10 * sigma))
            hist_plot, bins_plot = get_hist(energies, number_bins, left_plot, right_plot)
            # Plot from main data
            norm_1 = 1/max(hist_plot)
            plt.errorbar(bins_plot, hist_plot*norm_1, np.sqrt(hist_plot)*norm_1, fmt='.-',
                         capsize=5, zorder=5, label=label_1, color=color_1)
            # Plot from second data
            hist_2, bins_2 = get_hist(energies_2, number_bins, left_plot, right_plot)
            norm_2 = 1/max(hist_2)
            plt.errorbar(bins_2, hist_2*norm_2, np.sqrt(hist_2)*norm_2, fmt='.-',
                         capsize=5, zorder=5, label=label_2, color=color_2)
            # Plot from third data
            hist_3, bins_3 = get_hist(energies_3, number_bins, left_plot, right_plot)
            norm_3 = 1/max(hist_3)
            plt.errorbar(bins_3, hist_3*norm_3, np.sqrt(hist_3)*norm_3, fmt='.-',
                         capsize=5, zorder=5, label=label_3, color=color_3)
            # Plot where the shoulder should be
            E_reduced, distances, linestyles = calculate_distance_borders(bins_plot, hist_plot)
            for distance, linestyle in zip(distances, linestyles):
                E_new = E_reduced[distance]
                plt.axvline(x=E_new, linewidth=2, zorder=10, color='black', linestyle=linestyle,
                            label='Extra distance: %d cm' % (distance*100))
            # Stylise plot
            plt.grid(True, which='major', linestyle='--', zorder=0)
            plt.grid(True, which='minor', linestyle='--', zorder=0)
            plt.title('Peak at: %.2f meV (%.2f Å)' % (bins[peak], meV_to_A(bins[peak])))
            plt.xlabel('Energy [meV]')
            plt.ylabel('Counts (Normalized to maximum)')
            plt.xlim(x0 - (50 * sigma), x0 + (10 * sigma))
            plt.yscale('log')
            plt.legend(loc=2)
            # Save plot
            file_name = '%s_Peak_at_%.2f_meV_(%.2f_Å).pdf' % (label_1, bins[peak], meV_to_A(bins[peak]))
            output_path = output_folder + file_name
            fig.savefig(output_path, bbox_inches='tight')
            plt.close()
        except:
            print("Unexpected error:", sys.exc_info())
