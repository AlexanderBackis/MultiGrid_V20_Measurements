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

def analyze_Lineshape(ce_MG, detector_type, origin_voxel, number_bins):
    def Gaussian(x, a, x0, sigma):
        return a*np.exp(-(x-x0)**2/(2*sigma**2))

    plot_energy = True
    energies = calculate_energy(ce_MG, detector_type, origin_voxel)
    energy_hist, bin_centers = energy_plot(ce_MG, detector_type, origin_voxel,
                                           number_bins, 1, 10, plot_energy)
    # Find peaks
    threshold = 2000
    peaks, *_ = find_peaks(energy_hist, threshold=threshold)
    widths, *_ = peak_widths(energy_hist, peaks)
    #prominences, left_bases, right_bases = peak_prominences(energy_hist, peaks)
    plt.plot(bin_centers[left_bases], energy_hist[left_bases], color='green',
             zorder=5, linestyle='', marker='o')
    #plt.plot(bin_centers[right_bases], energy_hist[right_bases], color='blue',
    #         zorder=5, linestyle='', marker='x')
    plt.plot(bin_centers[peaks], energy_hist[peaks], color='red', zorder=5,
             linestyle='', marker='o')
    figs = []
    for width, peak, left_base, right_base, peak in zip(widths, peaks, left_bases, right_bases, peaks):
        left, right = bin_centers[peak]-width/2, bin_centers[peak]+width/2
        try:
            # Perform new histogram within estimated peak limits
            peak_energies = energies[(energies >= left) & (energies <= right)]
            peak_hist, peak_edges = np.histogram(peak_energies, bins=300, range=[left, right])
            peak_bin_centers = 0.5 * (peak_edges[1:] + peak_edges[:-1])
            # Check if there are more than one peak within, if so, split them
            popt, __ = curve_fit(Gaussian, peak_bin_centers, peak_hist)
            a, x0, sigma = popt[0], popt[1], abs(popt[2])
            xx = np.linspace(left, right, 100)
            fig_new = plt.figure()
            # Plot parameters from scipy to double-check
            plt.axvline(x=left, color='orange', linewidth=0.5)
            plt.axvline(x=right, color='orange', linewidth=0.5)
            # Plot Gaussian
            plt.plot(xx, Gaussian(xx, a, x0, sigma), color='purple')
            plt.plot(peak_bin_centers, peak_hist)
            plt.title('Peak at: %s meV' % bin_centers[peak])
            plt.xlabel('Energy [meV]')
            plt.ylabel('Counts')
            figs.append(fig_new)
        except:
            print("Unexpected error:", sys.exc_info()[0])

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
