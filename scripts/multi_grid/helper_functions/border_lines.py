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

import plotly.graph_objs as go

# =============================================================================
#                          Detector border lines
# =============================================================================


def initiate_detector_border_lines(mapping):
    """
    Produces a 3D hit-position histogram in (x, y, z)-coordinates, where the
    colorbar indicates number of counts at that specific coordinate.

    Args:
        detector_vec (list): List containing three elements, and each element is
                             a dictionary conatining the (wCh, gCh, bus)->(x, y, z)
                             mapping for a specific detector.

    Returns:
        b_traces (list): List containing the border traces for each of the three
                         detectors.
    """

    # Initiate all pairs of corners were lines will go between
    pairs = [[[80, 0], [80, 60]],
             [[80, 0], [80, 19]],
             [[80, 79], [80, 60]],
             [[80, 79], [80, 19]],
             [[119, 0], [119, 60]],
             [[119, 0], [119, 19]],
             [[119, 79], [119, 60]],
             [[119, 79], [119, 19]],
             [[80, 0], [119, 0]],
             [[80, 19], [119, 19]],
             [[80, 60], [119, 60]],
             [[80, 79], [119, 79]]
             ]
    # For each of the pairs, create a plotly trace with the plot
    b_traces = []
    for bus in range(0, 3):
        for pair in pairs:
            x_vec = []
            y_vec = []
            z_vec = []
            for loc in pair:
                gCh = loc[0]
                wCh = loc[1]
                coord = mapping[bus, gCh, wCh]
                x_vec.append(coord['x'])
                y_vec.append(coord['y'])
                z_vec.append(coord['z'])
            b_trace = go.Scatter3d(x=x_vec,
                                   y=y_vec,
                                   z=z_vec,
                                   mode='lines',
                                   line=dict(color='rgba(0, 0, 0, 0.5)',
                                             width=5))
            b_traces.append(b_trace)
    # Append a line showing how the incident neutrons travel
    neutron_line = go.Scatter3d(x=[0, 0],
                                y=[0, 0],
                                z=[28, 28.5],
                                mode='lines',
                                line=dict(color='rgb(0, 0, 500)',
                                          width=15),
                                text='Neutron Beam')
    #b_traces.append(neutron_line)
    return b_traces
