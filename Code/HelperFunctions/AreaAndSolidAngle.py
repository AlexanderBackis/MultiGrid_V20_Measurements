#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
AreaAndSolidAngle.py: Calculates the area and solid angle of the part of the
                      Multi-Grid currently being considered, as defined in the
                      filtering parameters.
"""

import numpy as np
from HelperFunctions.CreateMapping import create_mapping

# =============================================================================
#                      AREA AND SOLID ANGLE CALCULATIONS
# =============================================================================

def get_multi_grid_area_and_solid_angle(parameters):
    """
    Calculates the area and solid angle of the part of the Multi-Grid currently
    being considered, as defined in the filtering parameters. Solid angle is
    calculated assuming that detectors are at a circumference around the sample,
    and that only the up and down directions changes the solid angle. It is
    calculated by iterating through all surface voxels, and calculating the
    distance to the sample and angle theta according to:

                      theta = np.arctan(abs(y_vox/z_vox))
                      d = np.sqrt(x_vox ** 2 + y_vox ** 2 + z_vox ** 2)

    The projected area is then calculated according to:

                    projected_area = VOXEL_AREA * np.cos(theta)

    and the corresponding solid angle:

                    MG_solid_angle = (projected_area / (d ** 2))

    This calculation is repeated for all the surface voxels.



    Args:
        parameters (dict): Dictionary containing the parameters on how the data
                           is reduced. Here we are only interested in the
                           filters which affects the total surface area.

    Returns:
        MG_area (float): Active surface area of detector (i.e. active surface
                         sensitive to neutrons)

        MG_solid_angle (float): Solid angle corresponding to the active
                                detector surface area.
    """
    # Declare parameters (voxel area is in m^2)
    VOXEL_AREA = 0.0005434375
    mapping = create_mapping()
    # Get amount of surface to include
    modules_to_include = get_modules(parameters)
    grids = get_grids(parameters)
    # Iterate through surface and calculate area and solid angle
    MG_area = 0
    MG_projected_area = 0
    MG_solid_angle = 0
    for module in modules_to_include:
        wires = get_wires(module)
        for grid in grids:
            for wire in wires:
                # Extract coordinates
                vox_coordinate = mapping[module, grid, wire]
                x_vox = vox_coordinate['x']
                y_vox = vox_coordinate['y']
                z_vox = vox_coordinate['z']
                # Do calculations
                theta = np.arctan(abs(y_vox/z_vox))
                d = np.sqrt(x_vox ** 2 + y_vox ** 2 + z_vox ** 2)
                projected_area = VOXEL_AREA * np.cos(theta)
                MG_area += VOXEL_AREA
                MG_projected_area += projected_area
                MG_solid_angle += (projected_area / (d ** 2))
    return MG_area, MG_solid_angle


# =============================================================================
#                               HELPER FUNCTIONS
# =============================================================================

def get_modules(parameters):
    """
    Returns the buses which will be included in the area and solid angle
    calculations.

    Args:
        parameters (dict): Dictionary containing the parameters on how the data
                           is reduced. Here we are only interested in the
                           filters which affects the total surface area.
    Returns:
        buses (numpy array): Containing sequence of buses to be included in the
                             area and solid angle calculation

    """
    buses = [0, 1, 2]
    return buses


def get_wires(module):
    """
    Returns the wires which will be included in the area and solid angle
    calculations. Checks if any wire columns are disconnected, and if so,
    discards the area and solid angle contribution from these.

    Args:
        parameters (dict): Dictionary containing the parameters on how the data
                           is reduced. Here we are only interested in the
                           filters which affects the total surface area.

    Returns:
        wires (numpy array): Containing sequence of surface wires to be included
                             in the area and solid angle calculation
    """
    if module == 2:
        wires = [0, 20, 40]
    elif module == 3:
        wires = [0, 20, 60]
    elif module == 8:
        wires = [0, 40, 60]
    else:
        wires = [0, 20, 40, 60]
    return wires


def get_grids(parameters):
    """
    Returns the grids which will be included in the area and solid angle
    calculations.

    Args:
        parameters (dict): Dictionary containing the parameters on how the data
                           is reduced. Here we are only interested in the
                           filters which affects the total surface area.
    Returns:
        grids (numpy array): Containing sequence of grids to be included in the
                             area and solid angle calculation
    """
    if parameters['gCh'][2] is True:
        first_grid = parameters['gCh'][0] + 80 - 1
        last_grid = parameters['gCh'][1] + 80 - 1
        grids = np.arange(first_grid, last_grid+1)
    else:
        grids = np.arange(80, 120)
    return grids
