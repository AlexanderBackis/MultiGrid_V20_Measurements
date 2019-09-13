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
    coordinates = matrix[1:801]
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
    flip_bus_dict = {0: 2, 1: 1, 2: 0}
    return flip_bus_dict[bus]

def flip_wire(wCh):
        if 0 <= wCh <= 19:
            wCh += 60
        elif 20 <= wCh <= 39:
            wCh += 20
        elif 40 <= wCh <= 59:
            wCh -= 20
        elif 60 <= wCh <= 79:
            wCh -= 60
        return wCh

def get_new_x(x, y, theta):
    return np.cos(np.arctan(y/x)+theta)*np.sqrt(x ** 2 + y ** 2)

def get_new_y(x, y, theta):
    return np.sin(np.arctan(y/x)+theta)*np.sqrt(x ** 2 + y ** 2)
