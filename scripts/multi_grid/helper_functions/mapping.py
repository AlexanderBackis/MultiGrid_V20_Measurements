#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CreateMapping.py: Creates the channel->coordinate mapping using Isaaks CAD
                  mapping. Please consult the 'Isaaks_coordinate_explanation' in
                  the 'Information'-folder. Isaak start from lower right corner
                  in a different coordinate system than the one we use as global
                  system in instrument hall. After the mapping import has been
                  completed, the mapping is then translated into the correct
                  ordering, because buses are flipped and wire rows are flipped.
"""

import os
import numpy as np
import pandas as pd

# =============================================================================
#                                CREATE MAPPING
# =============================================================================

def create_mapping(origin_voxel, distance_offset=0):
    """
    Creates a channel->coordinate mapping for the detector. We place our origin
    at the beam center, where z=0 is at the location of the slit.

    Returns:
        detector_mapping (dict): Dictionary conta
    """

    # Calculate offset of detectors in (x, y, z)
    theta = 0
    bus_origin, gCh_origin, wCh_origin = origin_voxel
    temp_origin = {'x': 0, 'y': 0, 'z': 0}
    # Get offshift from corner of first voxel
    mapping_temp = create_channel_to_coordinate_map(theta, temp_origin)
    hit_coordinate = mapping_temp[bus_origin, gCh_origin, wCh_origin]
    x_offset = -hit_coordinate['x']
    y_offset = -hit_coordinate['y']
    z_offset = -hit_coordinate['z'] + 46.514e-3 + 28.366 + distance_offset
    offset = {'x': x_offset, 'y': y_offset, 'z': z_offset}
    mapping = create_channel_to_coordinate_map(theta, offset)
    return mapping


# =============================================================================
#                        INTERNAL DETECTOR MAPPING
# =============================================================================

def create_channel_to_coordinate_map(theta, offset):
    """
    Creates a channel->coordinate mapping for the ESS-type detectors. It needs
    a 'xlsx'-file containing the coordinates, using lower right outer detector
    corner as origin.

    Args:
        theta (float): Inclination of vessel in (x, z)-plane (instrument system)
        offset (float): Lower right corner deviation from sample origin

    Returns:
        ch_to_coord (dict): Dictionary conataining mapping for ESS detector,
                            according to (Bus, gCh, wCh)->(x, y, z)
    """
    # Import voxel coordinates
    dirname = os.path.dirname(__file__)
    file_path = os.path.join(dirname, '../../../tables/Mapping_MG_SEQ_ESS.xlsx')
    matrix = pd.read_excel(file_path).values
    coordinates = matrix[2:802]
    # Create empty dictionary to store mapping
    ch_to_coord = np.empty((3, 120, 80), dtype='object')
    coordinate = {'x': -1, 'y': -1, 'z': -1}
    axises =  ['x','y','z']
    # Iterate through all coordinates and start mapping procedure
    for i, row in enumerate(coordinates):
        grid_ch = i // 20 + 80
        for j, col in enumerate(row):
            module = j // 12
            layer = (j // 3) % 4
            wire_ch = (19 - (i % 20)) + (layer * 20)
            coordinate_count = j % 3
            coordinate[axises[coordinate_count]] = col
            if coordinate_count == 2:
                x = coordinate['x']
                y = coordinate['y']
                z = coordinate['z']
                # Convert from [mm] to [m]
                x = x/1000
                y = y/1000
                z = z/1000
                # Shift to match internal and external coordinate system
                z, x, y = x, y, z
                # Apply rotation
                x, z = get_new_x(x, z, theta), get_new_y(x, z, theta)
                # Apply translation
                x += offset['x']
                y += offset['y']
                z += offset['z']
                # Insert result in dictionary and reset temporary variable.
                # Flip bus and wire row locations, as we are parsing the
                # coordinates backwards.
                ch_to_coord[flip_bus(module), grid_ch, flip_wire(wire_ch)] = {'x': x, 'y': y, 'z': z}
                coordinate = {'x': -1, 'y': -1, 'z': -1}
    return ch_to_coord


# =============================================================================
#                            Helper Functions
# =============================================================================

def flip_bus(bus):
    """
    Flips the bus into correct position in detector, which is on the
    mirrored side (detector surface projections).

    Args:
        bus (int): Bus

    Returns:
        mirrored_bus (int): Mirrored bus (on detector surface projection)
    """
    flip_bus_dict = {0: 2, 1: 1, 2: 0}
    mirrored_bus = flip_bus_dict[bus]
    return mirrored_bus

def flip_wire(wCh):
    """
    Flips the wire channel into correct location in bus, which is on the
    mirrored side (detector surface projection).

    Args:
        wCh (int): Wire channel

    Returns:
        mirrored_wCh (int): Mirrored wire channel position (on detector
                            surface projection)
    """
    if 0 <= wCh <= 19:
        mirrored_wCh = wCh + 60
    elif 20 <= wCh <= 39:
        mirrored_wCh = wCh + 20
    elif 40 <= wCh <= 59:
        mirrored_wCh = wCh - 20
    elif 60 <= wCh <= 79:
        mirrored_wCh = wCh - 60
    return mirrored_wCh

def get_new_x(x, y, theta):
    """
    This works because we are working in Isaaks coordinate system, where
    each consecutive voxel center is placed further out in the x-y plane.
    We then check the angle for the non-rotated system (i.e. np.arctan(y/x)),
    followed by a an addition to theta (corresponding to the rotation, which
    is "up" from the x-axis, which makes it easy for us to add the angles).
    The summation of these two angles are then used to translate the coordinate
    into polar form and calculate the final position of the voxel.

    Args:
        x (float): x-position
        y (float): y-position
        theta (float): Detector inclination [Radians]

    Returns:
        new_x (float): New x in rotated geometry
    """
    new_x = np.cos(np.arctan(y/x)+theta)*np.sqrt(x ** 2 + y ** 2)
    return new_x

def get_new_y(x, y, theta):
    """
    This works because we are working in Isaaks coordinate system, where
    each consecutive voxel center is placed further out in the x-y plane.
    We then check the angle for the non-rotated system (i.e. np.arctan(y/x)),
    followed by a an addition to theta (corresponding to the rotation, which
    is "up" from the x-axis, which makes it easy for us to add the angles).
    The summation of these two angles are then used to translate the coordinate
    into polar form and calculate the final position of the voxel.

    Args:
        x (float): x-position
        y (float): y-position
        theta (float): Detector inclination [Radians]

    Returns:
        new_y (float): New y in rotated geometry
    """
    new_y = np.sin(np.arctan(y/x)+theta)*np.sqrt(x ** 2 + y ** 2)
    return new_y
