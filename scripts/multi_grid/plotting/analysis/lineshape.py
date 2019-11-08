#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Lineshape.py: Analyses the lineshape using our Figure-of-Merit
"""

import os
import matplotlib.pyplot as plt
import numpy as np


# =============================================================================
#                             PLOT FIGURE-OF-MERIT
# =============================================================================

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


# =============================================================================
#                           EXTRACT FIGURE-OF-MERIT
# =============================================================================

def get_FoM(beam, x0, sigma, start, end, bin_width):
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


# =============================================================================
#                               FIND SHOULDER
# =============================================================================

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
    # Calculate reduced energy from additional ToF (d is from closest voxel)
    ToF = d/average_v
    linestyles = ['solid', 'dotted', 'dashed', 'dashdot']
    E_reduced = {0.05: 0, 0.10: 0, 0.20: 0, 0.40: 0}
    for ToF_extra, distance in zip(ToF_extras, distances):
        E_new = get_new_E(d, ToF, ToF_extra)
        E_reduced[distance] = E_new
    return E_reduced, distances, linestyles
