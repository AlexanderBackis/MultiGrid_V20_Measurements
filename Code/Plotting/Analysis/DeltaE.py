#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
DeltaE.py: Function which histograms energy transfer data
"""

from HelperFunctions.EnergyTransfer import calculate_energy_transfer


# =============================================================================
#                             ENERGY TRANSFER
# =============================================================================


def energy_transfer_plot(df, Ei, number_bins):
    """
    Histograms the energy transfer values from a measurement

    Args:
        df (DataFrame): Clustered events
        Ei (float): Incident energy in meV
        number_bins (int): Number of bins to histogram energy transfer data

    Returns:
        fig (Figure): Figure containing nine 2D coincidences histograms, one
                      for each bus.
        dE_hist (numpy array): Numpy array containing the histogram data
        bin_centers (numpy array): Numpy array containing the bin centers
    """
    # Calculate DeltaE
    frame_shift = get_frame_shift(Ei)
    dE = calculate_energy_transfer(df, Ei, frame_shift)
    # Histogram DeltaE
    dE_hist, bin_edges = np.histogram(dE, bins=number_bins, range=[-Ei, Ei])
    bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
    # Plot data
    fig = plt.figure()
    plt.plot(bin_centers, dE_hist, '.-', color='black', zorder=5)
    plt.grid(True, which='major', linestyle='--', zorder=0)
    plt.grid(True, which='minor', linestyle='--', zorder=0)
    plt.xlabel('$E_i$ - $E_f$ [meV]')
    plt.ylabel('Counts')
    plt.yscale('log')
    plt.title('Energy transfer')
    return fig


# =============================================================================
#                             HELPER FUNCTIONS
# =============================================================================

def get_frame_shift(Ei):
    pass
