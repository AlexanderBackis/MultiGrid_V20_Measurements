#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Coincidences_3D.py: Produces a 3D hit-position histogram in
                    (x, y, z)-coordinates, where the colorbar indicates number
                    of counts at that specific coordinate.
"""

import plotly as py
import numpy as np
import plotly.graph_objs as go
import plotly.io as pio
import os

from HelperFunctions.CreateMapping import create_mapping

# =============================================================================
#                           Coincidence Histogram (3D)
# =============================================================================

def coincidences_3D_plot(df, detector_type, origin_voxel):
    """
    Produces a 3D hit-position histogram in (x, y, z)-coordinates, where the
    colorbar indicates number of counts at that specific coordinate.

    Args:
        df (DataFrame): Clustered events

    Yields:
        HTML file containing the 3D-histogram plot, automatically opened in the
        default browser.
    """

    # Declare max and min count
    min_count = 0
    max_count = np.inf
    # Perform initial filters
    df = df[(df.wCh != -1) & (df.gCh != -1)]
    # Initiate 'voxel_id -> (x, y, z)'-mapping
    mapping = create_mapping(detector_type, origin_voxel)
    # Initiate border lines
    b_traces = initiate_detector_border_lines(mapping, detector_type)
    # Calculate 3D histogram
    H, edges = np.histogramdd(df[['wCh', 'gCh', 'Bus']].values,
                              bins=(80, 40, 3),
                              range=((0, 80), (80, 120), (0, 3))
                              )
    # Insert results into an array
    hist = [[], [], [], []]
    loc = 0
    labels = []
    for wCh in range(0, 80):
        for gCh in range(80, 120):
            for bus in range(0, 3):
                over_min = H[wCh, gCh-80, bus] > min_count
                under_max = H[wCh, gCh-80, bus] <= max_count
                if over_min and under_max:
                    coord = mapping[bus, gCh, wCh]
                    hist[0].append(coord['x'])
                    hist[1].append(coord['y'])
                    hist[2].append(coord['z'])
                    hist[3].append(H[wCh, gCh-80, bus])
                    loc += 1
                    labels.append('Detector: ' + detector_type
                                  + '<br>'
                                  + 'Module: %d<br>' % bus
                                  + 'WireChannel: %d<br>' % wCh
                                  + 'GridChannel: %d<br>' % gCh
                                  + 'Counts: %d' % H[wCh, gCh-80, bus])
    # Produce 3D histogram plot
    MG_3D_trace = go.Scatter3d(x=hist[0],
                               y=hist[1],
                               z=hist[2],
                               mode='markers',
                               marker=dict(size=5,
                                           color=np.log10(hist[3]),
                                           colorscale='Jet',
                                           opacity=1,
                                           colorbar=dict(thickness=20,
                                                         title='log10(counts)'
                                                         ),
                                           #cmin=0,
                                           #cmax=2.5
                                           ),
                               text=labels,
                               name='Multi-Grid',
                               scene='scene1'
                               )
    # Introduce figure and put everything together
    fig = py.tools.make_subplots(rows=1, cols=1, specs=[[{'is_3d': True}]])
    # Insert histogram
    fig.append_trace(MG_3D_trace, 1, 1)
    # Insert vessel borders
    for b_trace in b_traces:
        fig.append_trace(b_trace, 1, 1)
    # Assign layout with axis labels, title and camera angle
    a = 1
    camera = dict(up=dict(x=0, y=0, z=1), center=dict(x=0, y=0, z=0),
                  eye=dict(x=-2*a, y=-0.5*a, z=1.3*a))
    scene=dict(camera=dict(eye=dict(x=1.15, y=1.15, z=0.8)), #the default values are 1.25, 1.25, 1.25
           xaxis=dict(),
           yaxis=dict(),
           zaxis=dict(),
           aspectmode='data', #this string can be 'data', 'cube', 'auto', 'manual'
           #a custom aspectratio is defined as follows:
           #aspectratio=dict(x=1, y=1, z=1)
           )
    fig['layout']['scene1']['xaxis'].update(title='x [m]')
    fig['layout']['scene1']['yaxis'].update(title='y [m]')
    fig['layout']['scene1']['zaxis'].update(title='z [m]')
    fig['layout'].update(title='Coincidences (3D)')
    fig['layout']['scene1'].update(scene)
    fig.layout.showlegend = False
    # If in plot He3-tubes histogram, return traces, else save HTML and plot
    py.offline.plot(fig, filename='../Output/Ce3Dhistogram.html', auto_open=True)
    pio.write_image(fig, '../Output/Ce3Dhistogram.pdf')


# =============================================================================
#                          Detector border lines
# =============================================================================


def initiate_detector_border_lines(mapping, detector_type):
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
    pairs_ESS = [[[80, 0], [80, 60]],
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
    pairs_ILL = [[[80, 0, 0], [80, 60, 2]],
                 [[80, 0, 0], [80, 19, 0]],
                 [[80, 79, 2], [80, 60, 2]],
                 [[80, 79, 2], [80, 19, 0]],
                 [[119, 0, 0], [119, 60, 2]],
                 [[119, 0, 0], [119, 19, 0]],
                 [[119, 79, 2], [119, 60, 2]],
                 [[119, 79, 2], [119, 19, 0]],
                 [[80, 0, 0], [119, 0, 0]],
                 [[80, 19, 0], [119, 19, 0]],
                 [[80, 60, 2], [119, 60, 2]],
                 [[80, 79, 2], [119, 79, 2]]
                 ]
    # For each of the pairs, create a plotly trace with the plot
    b_traces = []
    if detector_type == 'ESS':
        for bus in range(0, 3):
            for pair in pairs_ESS:
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

    else:
        for pair in pairs_ILL:
            x_vec = []
            y_vec = []
            z_vec = []
            for loc in pair:
                gCh = loc[0]
                wCh = loc[1]
                bus = loc[2]
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
    neutron_line = go.Scatter3d(x=[0, 0],
                                y=[0, 0],
                                z=[28, 28.7],
                                mode='lines',
                                line=dict(color='rgb(0, 0, 500)',
                                          width=10),
                                text='Neutron Beam')
    b_traces.append(neutron_line)
    return b_traces
