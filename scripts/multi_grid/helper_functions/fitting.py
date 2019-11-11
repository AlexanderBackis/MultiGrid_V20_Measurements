#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PathsAndFolders.py: Helper functions for handling of paths and folders.
"""

import numpy as np
from scipy.optimize import curve_fit

from multi_grid.helper_functions.misc import find_nearest

# =============================================================================
#                                 GAUSSIAN FIT
# =============================================================================

def fit_data(hist, bins, a_guess, x0_guess, sigma_guess):
    def Gaussian(x, a, x0, sigma):
        return a*np.exp(-(x-x0)**2/(2*sigma**2))

    popt, pcov = curve_fit(Gaussian, bins, hist, p0=[a_guess, x0_guess, sigma_guess])
    a, x0, sigma = popt[0], popt[1], abs(popt[2])
    xx = np.linspace(bins[0], bins[-1], 1000)
    return a, x0, sigma, xx, Gaussian(xx, a, x0, sigma), pcov


# =============================================================================
#                          GET ESTIMATION OF PARAMETERS
# =============================================================================

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

# =============================================================================
#                               HELPER FUNCTIONS
# =============================================================================

def get_hist(energies, number_bins, start, stop, weights=None):
    hist, bin_edges = np.histogram(energies, bins=number_bins,
                                   range=[start, stop], weights=weights)
    bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
    return hist, bin_centers
