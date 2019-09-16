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
#                        FULL MAPPING, ALL DETECTORS
# =============================================================================

def create_full_mapping():
    """
    Creates a channel->coordinate mapping for all three detectors.

    Returns:
        detector_mappings (list): List containing the mappings for all three
                                  detectors.
    """
    # Import measurement of corners
    ill_corners, ess_crb_corners, ess_pa_corners = get_detector_corners()
    # Calculate inclination of detectors in (x, z)-plane
    theta_ill = calculate_theta(ill_corners['left']['x'],
                                ill_corners['right']['x'],
                                ill_corners['left']['z'],
                                ill_corners['right']['z'])
    theta_crb = calculate_theta(ess_crb_corners['left']['x'],
                                ess_crb_corners['right']['x'],
                                ess_crb_corners['left']['z'],
                                ess_crb_corners['right']['z'])
    theta_pa = calculate_theta(ess_pa_corners['left']['x'],
                               ess_pa_corners['right']['x'],
                               ess_pa_corners['left']['z'],
                               ess_pa_corners['right']['z'])
    # Calculate offset of detectors in (x, y, z)
    offset_ill = ill_corners['right']
    offset_crb = ess_crb_corners['right']
    offset_pa = ess_pa_corners['right']
    # Calculate mappings
    ill_mapping = create_ill_channel_to_coordinate_map(theta_ill, offset_ill)
    crb_mapping = create_ess_channel_to_coordinate_map(theta_crb, offset_crb)
    pa_mapping = create_ess_channel_to_coordinate_map(theta_pa, offset_pa)
    detector_mappings = [ill_mapping, crb_mapping, pa_mapping]
    return detector_mappings


# =============================================================================
#                            ESS DETECTOR MAPPING
# =============================================================================

def create_ess_channel_to_coordinate_map(theta, offset):
    """
    Creates a channel->coordinate mapping for the ESS-type detectors. It needs
    a 'xlsx'-file containing the coordinates, using lower right outer detector
    corner as origin.

    Args:
        theta (float): Inclination of vessel in (x, z)-plane (instrument system)
        offset (float): Lower right corner deviation from sample origin

    Returns:
        ess_ch_to_coord (dict): Dictionary conataining mapping for ESS detector
    """
    # Import voxel coordinates
    dirname = os.path.dirname(__file__)
    file_path = os.path.join(dirname, '../../Tables/Mapping_MG_SEQ_ESS.xlsx')
    matrix = pd.read_excel(file_path).values
    coordinates = matrix[2:802]
    # Create empty dictionary to store mapping
    ess_ch_to_coord = np.empty((3, 120, 80), dtype='object')
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
                ess_ch_to_coord[flip_bus(module),
                                grid_ch,
                                flip_wire(wire_ch)] = {'x': x, 'y': y, 'z': z}
                coordinate = {'x': -1, 'y': -1, 'z': -1}
    return ess_ch_to_coord


# =============================================================================
#                            ILL DETECTOR MAPPING
# =============================================================================

def create_ill_channel_to_coordinate_map(theta, offset):
    """
    Creates a channel->coordinate mapping for the ESS-type detectors. It needs
    the offset from the lower right detector corner and the distances between
    consecutive voxels.

    Args:
        theta (float): Inclination of vessel in (x, z)-plane (instrument system)
        offset (float): Lower right corner deviation from sample origin

    Returns:
        ill_ch_to_coord (dict): Dictionary conataining mapping for ILL detector
    """
    # Declare distances between consecutive voxel centers in [mm]
    WireSpacing  = 10
    LayerSpacing = 23.5
    GridSpacing  = 23.5
    # Declare first voxel's offset
    x_offset = 46.514
    y_offset = 37.912
    z_offset = 37.95
    # Create dictionary to store data and generate full mapping
    ill_ch_to_coord = np.empty((3, 120, 80), dtype='object')
    for Bus in range(0, 3):
        for GridChannel in range(80, 120):
            for WireChannel in range(0, 80):
                    x = (WireChannel % 20)*WireSpacing + x_offset
                    y = ((WireChannel // 20)*LayerSpacing + (Bus*4*LayerSpacing) + y_offset)
                    z = ((GridChannel-80)*GridSpacing) + z_offset
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
                    # Insert result in dictionary, flip bus and wire row
                    # locations, as we are parsing the coordinates backwards.
                    ill_ch_to_coord[flip_bus(Bus),
                                    GridChannel,
                                    flip_wire(WireChannel)] = {'x': x, 'y': y, 'z': z}
    return ill_ch_to_coord

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

def calculate_theta(x1, x2, z1, z2):
    """
    Calculates the angle in a triangle with sides the length of (x2-x1) and
    (z2-z1).

    Args:
        x1 (float): First x-position
        x2 (float): Second x-position
        y1 (float): First y-position
        y2 (float): Second y-position

    Returns:
        theta (float): Angle theta in triangle with sides (x2-x1) and (z2-z1).
    """
    theta = np.arctan((z2-z1)/(x2-x1))
    return theta

def get_detector_corners():
    """
    Imports the measured detector corners of the three detectors and saves them
    in three dictionaries, which are stored in a list.

    Returns:
        theta (float): Angle theta in triangle with sides (x2-x1) and (z2-z1).
    """
    def create_corner_dict(corners):
        corners = {'left': {'x': corners[0],
                            'y': corners[1],
                            'z': corners[2]},
                   'right': {'x': corners[3],
                             'y': corners[4],
                             'z': corners[5]}}
        return corners
    # Import data from measurement of corners
    dirname = os.path.dirname(__file__)
    file_path = os.path.join(dirname, '../../Tables/Detector_Corners.xlsx')
    matrix = pd.read_excel(file_path).values
    # Create dictionary containing information about corners
    ill_corners = create_corner_dict(matrix[1][1:])
    ess_crb_corners = create_corner_dict(matrix[2][1:])
    ess_pa_corners = create_corner_dict(matrix[3][1:])
    return ill_corners, ess_crb_corners, ess_pa_corners
