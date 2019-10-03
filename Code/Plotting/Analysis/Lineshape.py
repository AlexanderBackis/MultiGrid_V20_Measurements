#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Lineshape.py: Analyses the lineshape using our Figure-of-Merit
"""
import sys
from Plotting.Analysis.DeltaE import energy_plot
from HelperFunctions.EnergyTransfer import calculate_energy
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks, peak_widths, peak_prominences
from scipy.optimize import curve_fit

# =============================================================================
#                         LINESHAPE INVESTIGATION
# =============================================================================

def analyze_Lineshape(ce_MG, detector_type, origin_voxel):
    def Gaussian(x, a, x0, sigma):
        return a*np.exp(-(x-x0)**2/(2*sigma**2))
    number_bins = 5000
    plot_energy = True
    energies = calculate_energy(ce_MG, detector_type, origin_voxel)
    energy_hist, bin_centers = energy_plot(ce_MG, detector_type, origin_voxel,
                                           number_bins, 1, 10, plot_energy)
    # Find peaks
    bins_part_1 = number_bins - number_bins//3
    bins_part_2 = number_bins//3
    heights_part1 = np.ones(bins_part_1)*12000
    heights_part2 = np.ones(bins_part_2)*1000
    heights = np.append(heights_part1, heights_part2)
    plt.plot(bin_centers[:bins_part_1], heights_part1, color='purple')
    plt.plot(bin_centers[bins_part_1:], heights_part2, color='purple')
    peaks, *_ = find_peaks(energy_hist, height=heights)
    widths, *_ = peak_widths(energy_hist, peaks)
    plt.plot(bin_centers[peaks], energy_hist[peaks], color='red', zorder=5,
             linestyle='', marker='o')
    figs = []
    print('Number of peaks: %d' % len(peaks))
    for width, peak in zip(widths, peaks):
        left, right = bin_centers[peak]-width/20, bin_centers[peak]+width/20
        left_fit, right_fit = bin_centers[peak]-width/60, bin_centers[peak]+width/60
        # Perform new histogram within estimated peak limits
        peak_bins = 100
        peak_energies = energies[(energies >= left) & (energies <= right)]
        peak_hist, peak_edges = np.histogram(peak_energies, bins=peak_bins, range=[left, right])
        peak_bin_centers = 0.5 * (peak_edges[1:] + peak_edges[:-1])
        # Define points to to fit procedure on
        peaks_in_peak, *_ = find_peaks(peak_hist, height=max(peak_hist)/2)
        prominences, left_bases, right_bases = peak_prominences(peak_hist, peaks_in_peak)
        fig_new = plt.figure()
        # Try fit
        fit_start = find_nearest(peak_bin_centers, left_fit)
        fit_stop = find_nearest(peak_bin_centers, right_fit)
        p0 = [6.03437070e+04, 6.19310485e+00, 4.02349443e-03]
        try:
            print('fit start: %s' % fit_start)
            print('fit stop: %s' % fit_stop)
            # Prepare guesses
            maximum = max(peak_hist)
            maximum_idx = find_nearest(peak_hist, maximum)
            half_maximum = maximum/2
            half_maximum_idx_1 = find_nearest(peak_hist[:maximum_idx],
                                              half_maximum)
            half_maximum_idx_2 = (find_nearest(peak_hist[maximum_idx:],
                                               half_maximum) + maximum_idx)
            FWHM = peak_bin_centers[half_maximum_idx_2] - peak_bin_centers[half_maximum_idx_1]
            a_guess = maximum
            x0_guess = peak_bin_centers[maximum_idx]
            sigma_guess = FWHM/(2*np.sqrt(2*np.log(2)))
            # Fit
            popt, __ = curve_fit(Gaussian,
                                 peak_bin_centers[fit_start:fit_stop],
                                 peak_hist[fit_start:fit_stop],
                                 p0=[a_guess, x0_guess, sigma_guess])
            a, x0, sigma = popt[0], popt[1], abs(popt[2])
            # Plot Gaussian
            xx = np.linspace(left, right, 1000)
            plt.plot(xx, Gaussian(xx, a, x0, sigma), color='purple')
        except:
            print("Unexpected error:", sys.exc_info())
        # Print fit edges
        plt.axvline(x=left_fit, color='green', linewidth=0.5)
        plt.axvline(x=right_fit, color='green', linewidth=0.5)
        # Plot parameters from scipy to double-check
        plt.axvline(x=left, color='orange', linewidth=0.5)
        plt.axvline(x=right, color='orange', linewidth=0.5)
        plt.plot(peak_bin_centers, peak_hist, color='blue', marker='.', linestyle='')
        plt.plot(peak_bin_centers[fit_start:fit_stop], peak_hist[fit_start:fit_stop],
                 color='red', marker='x', linestyle='')
        plt.title('Peak at: %s meV' % bin_centers[peak])
        plt.xlabel('Energy [meV]')
        plt.ylabel('Counts')
        figs.append(fig_new)

    for fig_temp in figs:
        fig_temp.show()


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
