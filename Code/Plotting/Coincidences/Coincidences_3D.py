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
from HelperFunctions.Borderlines import initiate_detector_border_lines

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
                                           #cmin=np.log10(min(hist[3])+1),
                                           #cmax=np.log10(max(hist[3])+2)
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
    a = 2
    camera = dict(up=dict(x=0, y=1, z=0), center=dict(x=0, y=0, z=0),
                  eye=dict(x=1*a, y=1*a, z=1*a))
    scene=dict(camera=camera, #the default values are 1.25, 1.25, 1.25
           xaxis=dict(),
           yaxis=dict(),
           zaxis=dict(),
           aspectmode='data')
    fig['layout']['scene1']['xaxis'].update(title='x [m]')
    fig['layout']['scene1']['yaxis'].update(title='y [m]')
    fig['layout']['scene1']['zaxis'].update(title='z [m]')
    fig['layout'].update(title='Coincidences (3D)')
    fig['layout']['scene1'].update(scene)
    fig.layout.showlegend = False
    py.offline.plot(fig, filename='../Output/Ce3Dhistogram.html', auto_open=True)
    pio.write_image(fig, '../Output/Ce3Dhistogram.pdf')
