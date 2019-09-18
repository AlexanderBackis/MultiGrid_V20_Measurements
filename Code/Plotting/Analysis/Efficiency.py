#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Efficiency.py: Calculates the efficiency through analysis of energy transfer
               spectra.
"""

from scipy.optimize import curve_fit
from HelperFunctions.EnergyTransfer import calculate_energy_transfer
from HelperFunctions.AreaAndSolidAngle import get_multi_grid_area_and_solid_angle


# =============================================================================
#                                EFFICIENCY
# =============================================================================

def calculate_efficiency(MG_dE_values, He3_dE_values, Ei, parameters):
    """
    Calculates the efficiency of the Multi-Grid detector at energy 'Ei'. Does
    through analysis in energy transfer spectra in three steps:

    1. Calculate number of counts in elastic peak, removing the background
    2. Get normalization on solid angle and, in the case of He-3, efficiency
    3. Normalize peak data and take fraction between Multi-Grid and He-3

    Args:
        MG_dE_values (numpy array): Energy transfer values from Multi-Grid
        He3_dE_values (numpy array): Energy transfer values from He-3 tubes
        Ei (float): Incident energy in meV
        parameters (dict): Dictionary containing the parameters on how the data
                           is reduced. Here we are only interested in the
                           filters which affects the total surface area.

    Returns:
        MG_efficiency (float): Efficiency of the Multi-Grid at energy Ei
    """
    # Get counts from peak in energy transfer spectra from MG and He-3
    MG_peak_counts = calculate_peak_area(MG_dE_values, Ei)
    He3_peak_counts = calculate_peak_area(He3_dE_values, Ei)
    # Get solid angle normalization
    __, MG_solid_angle = get_multi_grid_area_and_solid_angle(parameters)
    __, He3_solid_angle = get_He3_area_and_solid_angle()
    # Get efficiency normalization
    He3_efficiency = get_He3_efficiency(Ei)
    # Calculate Multi-Grid efficiency
    He3_normalized = (He3_peak_counts/He3_solid_angle) * (1/He3_efficiency)
    MG_normalized = MG_peak_counts/MG_solid_angle
    MG_efficiency = MG_normalized/He3_normalized
    return MG_efficiency


# =============================================================================
#                            CALCULATE PEAK AREA
# =============================================================================

def calculate_peak_area(dE_values, Ei):
    """
    Calculates the elastic peak area. Does this in three steps:

    1. Histogram dE values
    2. Fit peak and extract relevant fit parameters
    3. Sum counts within ±3 sigma of peak (peak counts)
    4. Get an estimate of background/meV at 5->7 sigma (THIS CAN DEFINITELY BE IMPROVED)
    5. Subtract background from peak counts

    Args:
        dE_values (numpy array): Energy transfer values
        Ei (float): Incident energy in meV

    Returns:
        signal_counts (float): Counts in elastic peak, background reduced
    """
    def Gaussian(x, a, x0, sigma):
        return a*np.exp(-(x-x0)**2/(2*sigma**2))
    # Histogram data, only use parts just around the peak: -Ei/3 -> Ei/3
    dE_hist, bin_edges = np.histogram(dE_values, bins=300, range=[-Ei/3, Ei/3])
    bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
    # Fit data, extract relevant fit parameters
    popt, __ = curve_fit(Gaussian, bin_centers, dE_hist)
    x0, sigma = popt[1], abs(popt[2])
    # Extract counts in peak
    left_edge, right_edge = (x0 - 3 * sigma), (x0 + 3 * sigma)
    dE_values_within_peak = dE_values[(dE_values >= left_edge) & (dE_values <= right_edge)]
    peak_counts = len(dE_values_within_peak)
    # Extract counts per meV at the background, between 5 -> 7 sigma
    background_left, background_right = (x0 + 5 * sigma), (x0 + 7 * sigma)
    dE_values_background = dE_values[(dE_values >= background_left) & (dE_values <= background_right)]
    background_counts = len(dE_values_background)
    background_counts_per_meV = background_counts / (background_right - background_left)
    # Calculate counts from signal
    signal_counts = peak_counts - (background_counts_per_meV * (right_edge - left_edge))
    return signal_counts
