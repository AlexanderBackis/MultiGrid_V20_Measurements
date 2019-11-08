#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EnergyHe3.py: Tools to calculate energy for the He-3 tube.
"""

import numpy as np

# =============================================================================
#                        CALCULATE ENERGY TRANSFER
# =============================================================================

def calculate_He3_energy(df, distance_offset=0):
    """
    Calculates the energy transfer of a data set. Distance is from source
    chopper to the center of the tube.

    Args:
        df (DataFrame): He-3 events

    Returns:
        energy (float): Energy in meV
    """
    # Declare necessary constants, neutron mass in [kg], DISTANCE in [m]
    NEUTRON_MASS = 1.674927351e-27 # [kg]
    JOULE_TO_meV = 6.24150913e18 * 1000
    DISTANCE = 28.239 + distance_offset # [m]
    time_offset = 0.6e-3 # [s]
    period_time = (1/14) # [s]
    ToF = (df.ToF * 8e-9 + time_offset) % period_time
    # Calculate energy Ef of detected neutron
    energy = ((NEUTRON_MASS/2) * ((DISTANCE/ToF) ** 2)) * JOULE_TO_meV
    return energy.values
