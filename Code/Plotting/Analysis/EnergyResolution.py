#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EnergyResolution.py: Calculates the energy resolution of the detector, using
                     the FWHM as a measure.
"""

import numpy as np
from scipy.optimize import curve_fit

# =============================================================================
#                                ENERGY RESOLUTION
# =============================================================================

def calculate_energy_resolution(dE_values, Ei):
    """
    Calculates the FWHM of the elastic peak area. Does this in three steps:

    1. Histogram dE values
    2. Fit peak and extract relevant fit parameters
    3. FWHM, according to http://mathworld.wolfram.com/GaussianFunction.html,
       is calculated using:

                       FWHM = 2 * np.sqrt(2 * log(2)) * sigma

    Args:
        dE_values (numpy array): Energy transfer values
        Ei (float): Incident energy in meV

    Returns:
        FWHM (float): FWHM of the elastic peak
    """
    def Gaussian(x, a, x0, sigma):
        return a*np.exp(-(x-x0)**2/(2*sigma**2))
    # Histogram data, only use parts just around the peak: -Ei/3 -> Ei/3
    dE_hist, bin_edges = np.histogram(dE_values, bins=300, range=[-Ei/3, Ei/3])
    bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
    # Fit data, extract relevant fit parameters
    popt, __ = curve_fit(Gaussian, bin_centers, dE_hist)
    x0, sigma = popt[1], abs(popt[2])
    # Calculate FWHM
    FWHM = 2 * np.sqrt(2*log(2)) * sigma
    return FWHM
