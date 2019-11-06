#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Efficiency.py: Calculates the efficiency through analysis of energy transfer
               spectra.
"""

import os
import matplotlib.pyplot as plt
import numpy as np

from HelperFunctions.Misc import find_nearest, A_to_meV, meV_to_A


# =============================================================================
#                                EFFICIENCY
# =============================================================================


def plot_efficiency(He3_energies, MG_energies,
                    He3_areas, MG_areas,
                    He3_err, MG_err,
                    monitor_norm_He3, monitor_norm_MG):
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
    # Iterate through energies to find matching efficiency from calculation to our measured data points
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



# =============================================================================
#                            CALCULATE PEAK AREA
# =============================================================================

def get_peak_area(energies, x0, sigma, bin_width):
    """
    Calculates the elastic peak area. Does this in five steps:

    1. Histogram dE values
    2. Fit peak and extract relevant fit parameters
    3. Sum counts within ±σ of peak (peak counts)
    4. Get an estimate of background/meV between "back_start" and "back_stop"
    5. Subtract background from peak counts

    Args:
        dE_values (numpy array): Energy transfer values
        Ei (float): Incident energy in meV

    Returns:
        signal_counts (float): Counts in elastic peak, background reduced
    """
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
    # Calculate area
    c = a - b * norm
    # Calculate uncertainites
    da = np.sqrt(a)
    db = np.sqrt(b)
    dc = np.sqrt(da ** 2 + (db*norm) ** 2)
    uncertainty = dc
    area = c
    # Plot background to cross-check calculation
    plt.axhline(y=b*(1/background_range_in_meV)*bin_width, color='black', linewidth=2, label=None)
    # Statistics for background
    plt.axvline(x=x0 + back_start*sigma, color='black', linewidth=2, label='Background')
    plt.axvline(x=x0 + back_stop*sigma, color='black', linewidth=2, label=None)
    # Statistics for peak area
    plt.axvline(x=x0 - sigma, color='orange', linewidth=2, label='-σ')
    plt.axvline(x=x0 + sigma, color='orange', linewidth=2, label='σ')
    return area, uncertainty
