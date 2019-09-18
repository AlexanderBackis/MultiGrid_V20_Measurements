#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PathsAndFolders.py: Helper functions for handling of paths and folders.
"""

from HelperFunctions.EnergyTransfer import calculate_energy_transfer


# =============================================================================
#                             ENERGY TRANSFER
# =============================================================================


def energy_transfer_plot(df, Ei, number_bins):
    # Calculate DeltaE
    frame_shift = get_frame_shift(Ei)
    dE = calculate_energy_transfer(df, Ei, frame_shift)
    # Histogram DeltaE
    dE_hist, bins = np.histogram(dE, bins=number_bins, range=[-Ei, Ei])
    bin_centers = 0.5 * (bins[1:] + bins[:-1])
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
