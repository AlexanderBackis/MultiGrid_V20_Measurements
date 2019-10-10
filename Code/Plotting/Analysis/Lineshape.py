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
from scipy.signal import find_peaks, peak_widths, peak_prominences
from scipy.optimize import curve_fit

from HeliumTubes.EnergyHe3 import calculate_He3_energy
from HeliumTubes.PlottingHe3 import energy_plot_He3

# =============================================================================
#                         LINESHAPE INVESTIGATION
# =============================================================================

def analyze_Lineshape(ce_MG, ce_He3, origin_voxel):
    """

    Non-Coated One Voxel: 800, 100
    Coated, Full Volume: 40000, 10000

    """
    def Gaussian(x, a, x0, sigma):
        return a*np.exp(-(x-x0)**2/(2*sigma**2))

    def meV_to_A(energy):
        return np.sqrt(81.81/energy)
    # Declare parameters
    number_bins = 5000
    plot_energy = True
    label_MG, label_He3 = 'Multi-Grid', 'He-3'
    useMaxNorm = True
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
        # Prepare new histograms within +/- 5 of our guessed sigma
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


# =============================================================================
#                                HELPER FUNCTIONS
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
