#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Lineshape.py: Analyses the lineshape using our Figure-of-Merit
"""

import sys
import os
from Plotting.Analysis.DeltaE import energy_plot
from HelperFunctions.EnergyTransfer import calculate_energy, get_distances
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.signal import find_peaks, peak_widths, peak_prominences
from scipy.optimize import curve_fit

from HelperFunctions.PathsAndFolders import mkdir_p
from HelperFunctions.Filtering import filter_clusters
from HelperFunctions.Fitting import fit_data, get_fit_parameters_guesses, get_hist

from HeliumTubes.EnergyHe3 import calculate_He3_energy
from HeliumTubes.PlottingHe3 import energy_plot_He3
from HeliumTubes.FilteringHe3 import filter_He3

# =============================================================================
#                         LINESHAPE INVESTIGATION
# =============================================================================


def analyze_all_lineshapes(origin_voxel, MG_filter_parameters, He3_filter_parameters):
    # Define parameters
    colors = {'MG_Coated': 'blue', 'MG_Non_Coated': 'green', 'He3': 'red'}
    monitor_norm_coated = 1/11411036
    monitor_norm_non_coated = 1/9020907
    monitor_norm_He3 = 1/10723199
    # Prepare data
    full_data = prepare_data(origin_voxel, MG_filter_parameters, He3_filter_parameters)
    MG_coated_data, MG_non_coated_data, He3_data = full_data[0], full_data[1], full_data[4]
    MG_coated_background, MG_non_coated_background, He3_background = full_data[2], full_data[3], full_data[5]

    # Plot all individual peaks
    #Coated_values = plot_all_peaks(MG_coated_data, 'MG_Coated', colors['MG_Coated'], 28.413)
    NonCoated_values =  plot_all_peaks(MG_non_coated_data, 'MG_Non_Coated', colors['MG_Non_Coated'], 28.413+1.5e-3)
    He3_values= plot_all_peaks(He3_data, 'He3', colors['He3'], 28.239+3e-3)
    # Extract values
    #energies_Coated, FoM_Coated, FoM_err_Coated, peak_areas_Coated, peak_err_Coated = Coated_values
    energies_NonCoated, FoM_NonCoated, FoM_err_NonCoated, peak_areas_NonCoated, peak_err_NonCoated = NonCoated_values
    energies_He3, FoM_He3, FoM_err_He3, peak_areas_He3, peak_err_He3 = He3_values
    # Plot Coated Radial blades compared to non-coated radial blades
    #plot_all_peaks_from_three_data_sets(MG_non_coated_data, 'MG_Non_Coated', colors['MG_Non_Coated'],
    #                                    MG_coated_data, 'MG_Coated', colors['MG_Coated'],
    #                                    He3_data, 'He3', colors['He3'])
    # Plot difference between NonCoated and Coated
    #plot_difference_from_two_data_sets(MG_non_coated_data, MG_coated_data, 'MG_Non_Coated', 'MG_Coated')


    # Plot important values
    #FWHMs = [FWHM_Coated, FWHM_NonCoated, FWHM_He3]
    #energies = [energies_Coated, energies_NonCoated, energies_He3]
    #labels = ['MG_Coated', 'MG_Non_Coated', 'He3']
    #FoMs = [FoM_Coated, FoM_NonCoated, FoM_He3]
    #errors = [FoM_err_Coated, FoM_err_NonCoated, FoM_err_He3]

    # Plot FWHM
    #fig = plt.figure()
    #for energy, FWHM, label in zip(energies, FWHMs, labels):
    #    plot_FWHM(energy, FWHM, label, colors[label])

    # Plot SEQUOIA data
    #dirname = os.path.dirname(__file__)
    #overview_folder = os.path.join(dirname, '../../../Tables/')
    #MG_SEQ_FWHM = np.loadtxt(overview_folder + 'MG_SEQ_FWHM.txt', delimiter=",")
    #MG_SEQ_Ei = np.loadtxt(overview_folder + 'MG_SEQ_Ei.txt', delimiter=",")
    #plt.plot(MG_SEQ_Ei, MG_SEQ_FWHM, marker='o', linestyle='-', zorder=5,
    #         label='Multi-Grid: SEQUOIA', color='black')
    #plt.legend()
    #fig.show()

    # Plot FoM
    #fig = plt.figure()
    #for energy, FoM, error, label in zip(energies, FoMs, errors, labels):
    #    plot_FoM(energy, FoM, error, label, colors[label])
    #plt.legend()
    #fig.show()

    # Plot efficiency
    fig = plt.figure()
    fig.set_figheight(5)
    fig.set_figwidth(15)
    plot_efficiency(np.array(energies_He3), np.array(energies_NonCoated),
                    np.array(peak_areas_He3), np.array(peak_areas_NonCoated),
                    np.array(peak_err_He3), np.array(peak_err_NonCoated),
                    monitor_norm_He3, monitor_norm_non_coated)
    fig.show()

def plot_FoM(energies, FoMs, errors, label, color):
    plt.title('Figure-of-Merit')
    plt.ylabel('FoM [shoulder/peak]')
    plt.xlabel('Peak energy [meV]')
    plt.grid(True, which='major', linestyle='--', zorder=0)
    plt.grid(True, which='minor', linestyle='--', zorder=0)
    plt.xscale('log')
    plt.yscale('log')
    plt.errorbar(energies, FoMs, errors, fmt='.-', capsize=5, zorder=5,
                 label=label, color=color)


def plot_FWHM(energies, FWHMs, label, color):
    # Plot V20 data
    plt.title('FWHM')
    plt.ylabel('FWHM [meV]')
    plt.xlabel('Peak energy [meV]')
    plt.grid(True, which='major', linestyle='--', zorder=0)
    plt.grid(True, which='minor', linestyle='--', zorder=0)
    plt.xscale('log')
    plt.yscale('log')
    plt.plot(energies, FWHMs, label=label, marker='o', linestyle='-',
             zorder=5, color=color)

def plot_efficiency(He3_energies, MG_energies,
                    He3_areas, MG_areas,
                    He3_err, MG_err,
                    monitor_norm_He3, monitor_norm_MG):
    # Load calculated efficiencies, as a function of lambda
    dirname = os.path.dirname(__file__)
    He3_efficiency_path = os.path.join(dirname, '../../../Tables/He3_efficiency.txt')
    MG_efficiency_path = os.path.join(dirname, '../../../Tables/MG_efficiency.txt')
    He3_efficiency = np.loadtxt(He3_efficiency_path, delimiter=",", unpack=True)
    MG_efficiency_calc = np.loadtxt(MG_efficiency_path, delimiter=",", unpack=True)[[0, 2]]
    # Remove elements in MG data which are not recorded in He-3
    MG_energies = np.delete(MG_energies, [0, 2])
    MG_areas = np.delete(MG_areas, [0, 2])
    MG_err = np.delete(MG_err, [0, 2])
    # Iterate through energies to find matching efficiency from calculation to
    # our measured data points
    He3_efficiency_datapoints = []
    for energy in He3_energies:
        # Save He3 efficiencies for data points
        idx = find_nearest(A_to_meV(He3_efficiency[0]), energy)
        He3_efficiency_datapoints.append(He3_efficiency[1][idx])
    He3_efficiency_datapoints = np.array(He3_efficiency_datapoints)
    # Rescale our curve to fit calibration
    idx = find_nearest(He3_efficiency[0], 2.5)
    calculated_efficiency_at_2_5_A = He3_efficiency[1][idx]
    He3_calculation_norm_upper = 0.964/calculated_efficiency_at_2_5_A
    He3_calculation_norm_average = 0.957/calculated_efficiency_at_2_5_A
    He3_calculation_norm_lower = 0.950/calculated_efficiency_at_2_5_A
    # Calculate average, as well as upper and lower bound for uncertainity estimation
    He3_efficiency_datapoints_upper = He3_efficiency_datapoints * He3_calculation_norm_upper
    He3_efficiency_datapoints_average = He3_efficiency_datapoints * He3_calculation_norm_average
    He3_efficiency_datapoints_lower = He3_efficiency_datapoints * He3_calculation_norm_lower
    # Calculated measured efficiency
    MG_efficiency = (MG_areas*monitor_norm_MG)/(He3_areas*(1/He3_efficiency_datapoints_average)*monitor_norm_He3)
    MG_efficiency_stat_unc = np.sqrt((MG_err/MG_areas) ** 2 + (He3_err/He3_areas) ** 2) * MG_efficiency
    # Calculate uncertainities
    MG_efficiency_upper = (MG_areas*monitor_norm_MG)/(He3_areas*(1/He3_efficiency_datapoints_upper)*monitor_norm_He3)
    MG_efficiency_lower = (MG_areas*monitor_norm_MG)/(He3_areas*(1/He3_efficiency_datapoints_lower)*monitor_norm_He3)
    upper_errors = MG_efficiency_upper - MG_efficiency + MG_efficiency_stat_unc
    lower_errors = MG_efficiency - MG_efficiency_lower + MG_efficiency_stat_unc
    full_errors = np.array([lower_errors, upper_errors])
    # Plot areas
    plt.subplot(1, 3, 1)
    plt.errorbar(He3_energies,
                 He3_areas*monitor_norm_He3,
                 He3_err*monitor_norm_He3,
                 fmt='.-', capsize=5,  color='red', label='He-3', zorder=5)
    plt.errorbar(MG_energies,
                 MG_areas*monitor_norm_MG,
                 MG_err*monitor_norm_MG,
                 fmt='.-', capsize=5,  color='blue', label='Multi-Grid', zorder=5)
    plt.xlabel('Energy (meV)')
    plt.ylabel('Peak area (Counts normalized by beam monitor counts)')
    plt.xlim(2, 120)
    plt.grid(True, which='major', linestyle='--', zorder=0)
    plt.grid(True, which='minor', linestyle='--', zorder=0)
    plt.title('Comparison MG and He-3')
    plt.legend()
    plt.xscale('log')
    plt.subplot(1, 3, 2)
    plt.xlabel('Energy (meV)')
    plt.ylabel('Efficiency')
    plt.xlim(2, 120)
    plt.plot(A_to_meV(MG_efficiency_calc[0]), MG_efficiency_calc[1], color='black',
             label='MG (90° incident angle)', zorder=5)
    plt.errorbar(MG_energies, MG_efficiency, full_errors, fmt='.-',
                capsize=5, color='blue', label='Measured MG efficiency', zorder=5)
    #plt.plot(He3_energies, He3_efficiency_datapoints, color='red',
    #         marker='o', linestyle='', label='He-3, Calculated', zorder=5)
    plt.grid(True, which='major', linestyle='--', zorder=0)
    plt.grid(True, which='minor', linestyle='--', zorder=0)
    plt.title('Efficiency measurement')
    plt.xscale('log')
    plt.legend()
    plt.subplot(1, 3, 3)
    plt.xlabel('Wavelength (Å)')
    plt.ylabel('Efficiency')
    MG_efficiency = (MG_areas*monitor_norm_MG)/(He3_areas*(1/He3_efficiency_datapoints)*monitor_norm_He3)
    MG_efficiency_unc = np.sqrt((MG_err/MG_areas) ** 2 + (He3_err/He3_areas) ** 2) * MG_efficiency
    plt.errorbar(meV_to_A(MG_energies), MG_efficiency, full_errors, fmt='.-',
                 capsize=5, color='blue', label='Measured MG efficiency', zorder=5)
    plt.grid(True, which='major', linestyle='--', zorder=0)
    plt.grid(True, which='minor', linestyle='--', zorder=0)
    plt.plot(MG_efficiency_calc[0], MG_efficiency_calc[1], color='black',
             label='MG (90° incident angle)', zorder=5)
    #plt.plot(meV_to_A(He3_energies), He3_efficiency_datapoints, color='red',
    #         marker='o', linestyle='', label='He-3, Calculated', zorder=5)
    plt.title('Efficiency measurement')
    plt.legend()
    plt.tight_layout()



def plot_difference_from_two_data_sets(data_1, data_2, label_1, label_2):
    # Prepare output paths
    dirname = os.path.dirname(__file__)
    output_folder = os.path.join(dirname, '../../../Output/Difference_%s_and_%s/' % (label_1, label_2))
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
            a, x0, sigma, *_ = fit_data(hist_fit, bins_fit, a_guess, x0_guess, sigma_guess)
            fig = plt.figure()
            plot_sigma_borders(x0, sigma)
            # Prepare data within +/- 25 of our estimated sigma, we'll use this to plot
            left_plot, right_plot = (x0 - (31 * sigma)), (x0 + (31 * sigma))
            hist_plot, bins_plot = get_hist(energies, number_bins, left_plot, right_plot)
            # Plot from main data
            norm_1 = 1/max(hist_plot)
            # Plot from second data
            hist_2, bins_2 = get_hist(energies_2, number_bins, left_plot, right_plot)
            norm_2 = 1/max(hist_2)
            difference = hist_plot*norm_1 - hist_2*norm_2
            error_1 = np.sqrt(hist_plot)*norm_1
            error_2 = np.sqrt(hist_fit)*norm_2
            error_full = np.sqrt(error_1 ** 2 + error_2 ** 2)
            plt.errorbar(bins_plot, difference, error_full, fmt='.-',
                         capsize=5, zorder=5, label='NonCoated - Coated', color='black')
            # Stylise plot
            plt.grid(True, which='major', linestyle='--', zorder=0)
            plt.grid(True, which='minor', linestyle='--', zorder=0)
            plt.title('NonCoated Radial - Coated Radial, Peak at: %.2f meV (%.2f Å)' % (bins[peak], meV_to_A(bins[peak])))
            plt.xlabel('Energy [meV]')
            plt.ylabel('Counts (Normalized to maximum)')
            plt.xlim(x0 - (31 * sigma), x0 + (31 * sigma))
            #plt.yscale('log')
            plt.legend(loc=1)
            # Save plot
            file_name = '%s_Peak_at_%.2f_meV_(%.2f_Å).pdf' % (label_1, bins[peak], meV_to_A(bins[peak]))
            output_path = output_folder + file_name
            fig.savefig(output_path, bbox_inches='tight')
            plt.close()
        except:
            print("Unexpected error:", sys.exc_info())


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
            a, x0, sigma, *_ = fit_data(hist_fit, bins_fit, a_guess, x0_guess, sigma_guess)
            fig = plt.figure()
            plot_sigma_borders(x0, sigma)
            # Prepare data within +/- 25 of our estimated sigma, we'll use this to plot
            left_plot, right_plot = (x0 - (31 * sigma)), (x0 + (31 * sigma))
            hist_plot, bins_plot = get_hist(energies, number_bins, left_plot, right_plot)
            # Plot from beam
            norm_beam = 1/duration
            plt.errorbar(bins_plot, hist_plot, np.sqrt(hist_plot), fmt='.-',
                         capsize=5, zorder=5, label=label_data, color=color_data)
            # Plot from background
            hist_background, bins_background = get_hist(energies_background, number_bins, left_plot, right_plot)
            norm_background = 1/duration_background
            #plt.errorbar(bins_background, hist_background*norm_background,
            #             np.sqrt(hist_background)*norm_background, fmt='.-',
            #             capsize=5, zorder=5, label='Background', color='black')
            plt.axvline(x=x0 - 30*sigma, color='black', linewidth=2, label='-30σ')
            plt.axvline(x=x0 - 25*sigma, color='black', linewidth=2, label='-25σ')
            # Stylise plot
            plt.grid(True, which='major', linestyle='--', zorder=0)
            plt.grid(True, which='minor', linestyle='--', zorder=0)
            plt.title('Peak at: %.2f meV (%.2f Å)' % (bins[peak], meV_to_A(bins[peak])))
            plt.xlabel('Energy [meV]')
            plt.ylabel('Counts')
            plt.xlim(x0 - (31 * sigma), x0 + (31 * sigma))
            plt.yscale('log')
            plt.legend(loc=1)
            # Store important values
            bin_width = bins_plot[1] - bins_plot[0]
            FWHM = 2 * np.sqrt(2*np.log(2)) * sigma
            FoM, uncertainity = get_FoM(energies, x0, sigma, bin_width)
            FWHM_list.append(FWHM)
            FoM_list.append(FoM)
            uncertainites.append(uncertainity)
            peak_energies.append(bins[peak])
            # Save plot
            file_name = '%s_Peak_at_%.2f_meV_(%.2f_Å).pdf' % (label_data, bins[peak], meV_to_A(bins[peak]))
            output_path = output_folder + file_name
            fig.savefig(output_path, bbox_inches='tight')
            plt.close()
        except:
            print("Unexpected error:", sys.exc_info())
    return FWHM_list, FoM_list, uncertainites, peak_energies


def plot_all_peaks_from_three_data_sets(data_1, label_1, color_1,
                                        data_2, label_2, color_2,
                                        data_3, label_3, color_3):
    # Prepare output paths
    dirname = os.path.dirname(__file__)
    output_folder = os.path.join(dirname, '../../../Output/Comparison_%s_and_%s/' % (label_1, label_2))
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



def plot_all_peaks(data, label, color, chopper_to_detector_distance):
    # Prepare output paths
    dirname = os.path.dirname(__file__)
    output_folder = os.path.join(dirname, '../../../Output/%s/' % label)
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
        #for distance, linestyle in zip(distances[1:3], linestyles[1:3]):
        #    E_new = reduced_energies[distance]
        #    plt.axvline(x=E_new, linewidth=2, zorder=10, color='black', linestyle=linestyle,
        #                label='Extra distance: %d cm' % (distance*100))
        # Extract FoM
        start, end = reduced_energies[0.20], reduced_energies[0.10]
        bin_width = bins_plot[1] - bins_plot[0]
        FoM, FoM_uncertainity = get_distance_FoM(energies, x0, sigma, start, end, bin_width)
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
    get_heights = [get_heights_Coated, get_heights_NonCoated]
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
    MG_distance_offsets = [1.5e-3, 0, 1.5e-3, 0]
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
        energies = calculate_energy(df_red, origin_voxel, MG_distance_offsets[i])
        hist, bins = get_hist(energies, number_bins, start, end)
        data = [duration, energies, hist, bins]
        if i < 2:
            # If it is a beam measurement, extract peaks
            peaks = get_peaks(hist, heights_vec_MG[i], number_bins)
            #heights = get_heights[i](bins)
            #peaks, *_ = find_peaks(hist, height=heights)
            widths, *_ = peak_widths(hist, peaks)
            data.extend([peaks, widths])
        full_data.append(data)
    # Store He-3 data
    He3_durations = [54304, 58094]
    He3_distance_offset = 3e-3
    print('He-3...')
    for i, (file_name, duration) in enumerate(zip(He3_file_names, He3_durations)):
        path = os.path.join(dirname, '../../../Data/Lineshape/%s' % file_name)
        df = pd.read_hdf(path, 'df')
        df_red = filter_He3(df, He3_filter_parameters)
        energies = calculate_He3_energy(df_red, He3_distance_offset)
        hist, bins = get_hist(energies, number_bins, start, end)
        data = [duration, energies, hist, bins]
        if i < 1:
            # If it is a beam measurement, extract peaks
            peaks = get_peaks(hist, heights_He3, number_bins)
            #heights = get_heights_He3(bins)
            #peaks, *_ = find_peaks(hist, height=heights)
            widths, *_ = peak_widths(hist, peaks)
            data.extend([peaks, widths])
        full_data.append(data)
    return full_data


# =============================================================================
#                                HELPER FUNCTIONS
# =============================================================================

def get_peak_area(energies, x0, sigma, bin_width):
    # Get background range
    if x0 < 70:
        back_start = -30
        back_stop = -25
    else:
        back_start = -3
        back_stop = -2
    # Extract number of counts from regions of interest
    peak_counts = energies[(energies >= (x0 - sigma)) & (energies <= (x0 + sigma))]
    background_counts = energies[(energies >= (x0 + back_start*sigma)) & (energies <= (x0 + back_stop*sigma))]
    # Rename for easier calculation of uncertainties
    a = len(peak_counts)
    b = len(background_counts)
    background_range_in_meV = sigma*abs((back_start-back_stop))
    # Define normalization constants
    norm = (1/background_range_in_meV) * 2*sigma
    # Calculate FoM
    c = a - b * norm
    # Calculate uncertainites
    da = np.sqrt(a)
    db = np.sqrt(b)
    dc = np.sqrt(da ** 2 + (db*norm) ** 2)
    uncertainty = dc
    area = c
    # Plot background to cross-check calculation
    plt.axhline(y=b*(1/background_range_in_meV)*bin_width,
                color='black', linewidth=2, label=None)
    # Statistics for background
    plt.axvline(x=x0 + back_start*sigma, color='black', linewidth=2, label='Background')
    plt.axvline(x=x0 + back_stop*sigma, color='black', linewidth=2, label=None)
    # Statistics for peak area
    plt.axvline(x=x0 - sigma, color='orange', linewidth=2, label='-σ')
    plt.axvline(x=x0 + sigma, color='orange', linewidth=2, label='σ')
    return area, uncertainty


def get_FoM(beam, x0, sigma, bin_width):
    # Extract number of counts from regions of interest
    peak_beam_counts = beam[(beam >= (x0 - 3*sigma)) & (beam <= (x0 + 3*sigma))]
    shoulder_beam_counts = beam[(beam >= (x0 - 5*sigma)) & (beam <= (x0 - 3*sigma))]
    background_estimation = beam[(beam >= (x0 - 30*sigma)) & (beam <= (x0 - 25*sigma))]
    # Rename for easier calculation of uncertainties
    a = len(peak_beam_counts)
    b = len(shoulder_beam_counts)
    c = len(background_estimation)
    d_meV = 5*sigma
    # Subtract background
    norm_d = ((1/d_meV) * 6 * sigma)
    norm_e = ((1/d_meV) * 2 * sigma)
    d = a - c * norm_d
    e = b - c * norm_e
    # Calculate FoM
    f = e/d
    # Calculate uncertainites
    da = np.sqrt(a)
    db = np.sqrt(b)
    dc = np.sqrt(c)
    dd = np.sqrt(da ** 2 + (dc*norm_d) ** 2)
    de = np.sqrt(db ** 2 + (dc*norm_e) ** 2)
    df = np.sqrt((dd/d) ** 2 + (de/e) ** 2)
    uncertainty = df * f
    FoM = f
    # Plot background to cross-check calculation
    plt.axhline(y=c*(1/d_meV)*bin_width, color='black', linewidth=2, label=None)
    return FoM, uncertainty

def get_distance_FoM(beam, x0, sigma, start, end, bin_width):
    # Define norm from range
    norm_range = 1/(end-start)
    # Extract number of counts from regions of interest
    shoulder_counts = beam[(beam >= start) & (beam <= end)]
    background_counts = beam[(beam >= (x0 - 30*sigma)) & (beam <= (x0 - 25*sigma))]
    peak_counts = beam[(beam >= (x0 - sigma)) & (beam <= (x0 + sigma))]
    # Rename for easier calculation of uncertainties
    a = len(shoulder_counts)
    b = len(background_counts)
    c = len(peak_counts)
    background_range_in_meV = 5*sigma
    # Define normalization constants
    norm_shoulder = (1/background_range_in_meV) * (end-start)
    norm_peak = (1/background_range_in_meV) * (2 * sigma)
    # Calculate FoM
    d = a - b * norm_shoulder
    e = c - b * norm_peak
    f = d/e
    # Calculate uncertainites
    da = np.sqrt(a)
    db = np.sqrt(b)
    dc = np.sqrt(c)
    dd = np.sqrt(da ** 2 + (db*norm_shoulder) ** 2)
    de = np.sqrt(dc ** 2 + (db*norm_peak) ** 2)
    df = np.sqrt((dd/d) ** 2 + (de/e) ** 2)
    uncertainty = df * f * norm_range
    FoM = f * norm_range
    # Plot background to cross-check calculation
    #plt.axhline(y=b*(1/background_range_in_meV)*bin_width,
    #            color='black', linewidth=2, label=None)
    #plt.axvline(x=x0 - 30*sigma, color='orange', linewidth=2, label='Background')
    #plt.axvline(x=x0 - 25*sigma, color='orange', linewidth=2, label=None)
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
    # Get peaks
    peaks, *_ = find_peaks(hist, height=heights)
    return peaks

def plot_sigma_borders(x0, sigma):
    plt.axvline(x=x0 - 5*sigma, color='orange', linewidth=2, label='-5σ')
    plt.axvline(x=x0 - 3*sigma, color='purple', linewidth=2, label='-3σ')
    plt.axvline(x=x0 - sigma, color='green', linewidth=2, label='-σ')
    plt.axvline(x=x0 + sigma, color='green', linewidth=2, label='σ')

def calculate_distance_borders(bins, hist, d=28.413):
    def E_to_v(energy_in_meV):
        # Define constants
        JOULE_TO_meV = 6.24150913e18 * 1000
        meV_TO_JOULE = 1/JOULE_TO_meV
        NEUTRON_MASS = 1.674927351e-27
        # Calculate velocity of neutron
        v = np.sqrt((2*energy_in_meV*meV_TO_JOULE)/NEUTRON_MASS)
        return v

    def get_new_E(d, ToF, ToF_extra):
        # Define constants
        JOULE_TO_meV = 6.24150913e18 * 1000
        NEUTRON_MASS = 1.674927351e-27
        E_new = ((NEUTRON_MASS/2)*(d/(ToF+ToF_extra)) ** 2) * JOULE_TO_meV
        return E_new


    # Declare intervals, in m
    distances = np.array([5, 10, 20, 40]) * 1e-2
    # Extract average E
    average_E = bins[hist == max(hist)]
    average_v = E_to_v(average_E)
    # Calculate additional ToF for these distances
    ToF_extras = distances / average_v
    print('ToF extras')
    print(ToF_extras)
    # Calculate reduced energy from additional ToF (d is from closest voxel)
    ToF = d/average_v
    print('---')
    print('Average energy: %f' % average_E)
    print('Average v: %f' % average_v)
    print('Test to get back first energy: %f' % get_new_E(d, ToF, 0))
    linestyles = ['solid', 'dotted', 'dashed', 'dashdot']
    E_reduced = {0.05: 0, 0.10: 0, 0.20: 0, 0.40: 0}
    for ToF_extra, distance in zip(ToF_extras, distances):
        E_new = get_new_E(d, ToF, ToF_extra)
        E_reduced[distance] = E_new
        print('Reduced energy: %f' % E_new)
    return E_reduced, distances, linestyles








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

def get_heights_NonCoated(bin_centers):
    heights = np.zeros(len(bin_centers))
    heights[(bin_centers >= 0) & (bin_centers <= 1.1)] = 12000
    heights[(bin_centers >= 1.1) & (bin_centers <= 3.0)] = 7500
    heights[(bin_centers >= 3.0) & (bin_centers <= 50.0)] = 12000
    heights[(bin_centers >= 50.0) & (bin_centers <= 70.0)] = 1000
    heights[(bin_centers >= 70.0) & (bin_centers <= 100.0)] = 85
    return heights

def get_heights_Coated(bin_centers):
    heights = np.zeros(len(bin_centers))
    heights[(bin_centers >= 0) & (bin_centers <= 1.1)] = 20000
    heights[(bin_centers >= 1.1) & (bin_centers <= 3.0)] = 10000
    heights[(bin_centers >= 3.0) & (bin_centers <= 50.0)] = 20000
    heights[(bin_centers >= 50.0) & (bin_centers <= 70.0)] = 2000
    heights[(bin_centers >= 70.0) & (bin_centers <= 100.0)] = 105
    return heights

def get_heights_He3(bin_centers):
    heights = np.zeros(len(bin_centers))
    heights[(bin_centers >= 0) & (bin_centers <= 2.5)] = 2000
    heights[(bin_centers >= 2.5) & (bin_centers <= 50.0)] = 20000
    heights[(bin_centers >= 50.0) & (bin_centers <= 70.0)] = 2000
    heights[(bin_centers >= 70.0) & (bin_centers <= 100.0)] = 286
    return heights



def meV_to_A(energy):
    return np.sqrt(81.81/energy)


def A_to_meV(wavelength):
    return (81.81/(wavelength ** 2))

def get_duration(df):
    times = df.Time.values
    duration_in_seconds = (times[-1] - times[0]) * 62.5e-9
    return duration_in_seconds
