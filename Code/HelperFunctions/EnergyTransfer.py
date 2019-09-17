#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EnergyTransfer.py: Tools to calculate energy transfer.
"""

import numpy as np
from HelperFunctions.CreateMapping import create_full_mapping

NEUTRON_MASS = 1.674927351e-27  # Neutron mass [kg]
meV_TO_JOULE = 1.60218e-19 * 0.001  # Convert meV to Joule
JOULE_TO_meV = 6.24150913e18 * 1000  # Convert Joule to meV
CHOPPER_TO_SAMPLE = 20.01 # Chopper to detector distance [m]

# =============================================================================
#                        CALCULATE ENERGY TRANSFER
# =============================================================================

def calculate_energy_transfer(df, Ei, frame_shift):
    """
    Calculates the energy transfer of a data set.

    Args:
        df (DataFrame): Clustered events
        Ei (float): Incident neutron energy in meV
        frame_shift (float): Frameshift in seconds

    Returns:
        DeltaE (float): Energy transfer in meV
    """

    # Get distance matric for each voxel
    distance_mapping = get_distances()
    # Extract necessary values from dataframe
    wChs = df.wCh
    gChs = df.gCh
    buses = df.Bus
    ToFs = df.ToF
    distances = distance_mapping[buses, gChs, wChs]
    # Calculate time spend travelling chopper-sample
    vi = np.sqrt((E_i * meV_TO_JOULE * 2)/NEUTRON_MASS)
    ti = (CHOPPER_TO_SAMPLE/ vi)
    # Calculate time spend travelling sample-detector
    tf = ToF * 62.5e-9 + frame_shift - ti
    # Filter away events with negative travel time
    tf = tf[tf > 0]
    # Calculate energy Ef of neutron after scattering from sample
    Ef = ((NEUTRON_MASS/2) * ((d/tf) ** 2)) * JOULE_TO_meV
    # Calculate energy transfer
    DeltaE = E_i - E_f
    return DeltaE



# =============================================================================
#                           VOXEL DISTANCES 3D MATRIX
# =============================================================================

def get_distances():
    """
    Calculates the distances to each of the voxels from the sample.

    Returns:
        distances (numpy array): 3D matrix where the axises are bus, gCh and wCh.
                                 The elements are the distances in m from the
                                 sample position.
    """
    detector_mappings = create_full_mapping()
    distances = np.zeros((9, 120, 80), dtype='float')
    for bus in range(0, 9):
        detector_mapping = detector_mappings(bus//3)
        for gCh in range(80, 120):
            for wCh in range(0, 80):
                coordinate = detector[bus%3, grid, wire]
                x = vox_coordinate['x']
                y = vox_coordinate['y']
                z = vox_coordinate['z']
                distances[bus, gCh, wCh] = np.sqrt(x ** 2 + y ** 2 + z ** 2)
    return distances
