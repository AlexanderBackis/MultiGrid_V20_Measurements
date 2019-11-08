#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Animation3D.py: Helper functions for handling of paths and folders.
"""

import plotly as py
import numpy as np
import plotly.graph_objs as go
import plotly.io as pio
import os
import imageio
import shutil

from multi_grid.helper_functions.mapping import create_mapping
from multi_grid.helper_functions.border_lines import initiate_detector_border_lines
from multi_grid.helper_functions.misc import mkdir_p

# =============================================================================
#                         3D Animation - Time sweep
# =============================================================================

def Animation_3D_plot(ce_full, origin_voxel):
    # Filters
    ce = ce_full[ce_full.Bus == 1]
    ce = ce[(ce.wCh != -1) & (ce.gCh != -1)]
    # Time intervals
    iter_start = 0
    iter_stop = 60
    step = 1
    conversion_factor = 1/(62.5e-9)
    # Mapping and Borders
    mapping = create_mapping(origin_voxel)
    borderline_traces = initiate_detector_border_lines(mapping)
    # Plotting
    min_count = 0
    max_count = np.inf
    cmin = 1
    cmax = 10
    # Storage
    dir_name = os.path.dirname(__file__)
    temp_animation_folder = os.path.join(dir_name, '../temp_animation_folder/')
    output_path = os.path.join(dir_name, '../3D_animation.gif')
    mkdir_p(temp_animation_folder)
    # Iteration
    delimiters = np.arange(iter_start, iter_stop, step)
    for i, delimiter in enumerate(delimiters):
        # Get current limits and filter data
        start = delimiter * conversion_factor
        stop = (delimiter + step) * conversion_factor
        ce_step = ce[(ce.Time >= start) & (ce.Time <= stop)]
        # Plot
        traces = []
        hist = get_3D_histogram(ce_step, mapping, min_count, max_count)
        MG_3D_trace = get_MG_3D_trace(hist, cmin, cmax)
        traces.append(MG_3D_trace)
        traces.extend(borderline_traces)
        fig = create_fig(traces, start, stop)
        # Save
        path = temp_animation_folder + '%d.png' % i
        pio.write_image(fig, path)
    # Animate
    images = []
    files = os.listdir(temp_animation_folder)
    files = [file[:-4] for file in files if file[-9:] != '.DS_Store' and file != '.gitignore']
    for file in sorted(files, key=int):
        images.append(imageio.imread(temp_animation_folder + file + '.png'))
    imageio.mimsave(output_path, images)
    shutil.rmtree(temp_animation_folder, ignore_errors=True)



# =============================================================================
#                                Helper Functions
# =============================================================================

def get_3D_histogram(ce, mapping, min_count, max_count):
    # Calculate 3D histogram
    H, edges = np.histogramdd(ce[['wCh', 'gCh', 'Bus']].values,
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
    return hist

def get_MG_3D_trace(hist, cmin=None, cmax=None):
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
                                           cmin=cmin,
                                           cmax=cmax
                                           ),
                               name='Multi-Grid',
                               scene='scene1'
                               )
    return MG_3D_trace

def create_fig(traces, start, stop):
    # Introduce figure and put everything together
    fig = py.tools.make_subplots(rows=1, cols=1, specs=[[{'is_3d': True}]])
    # Insert vessel borders
    for trace in traces:
        fig.append_trace(trace, 1, 1)
    # Assign layout with axis labels, title and camera angle
    a = 0.6
    camera = dict(up=dict(x=0, y=1, z=0), center=dict(x=-3, y=-6, z=-2),
                  eye=dict(x=1*a, y=1*a, z=1*a))
    scene=dict(camera=camera, #the default values are 1.25, 1.25, 1.25
           xaxis=dict(),
           yaxis=dict(),
           zaxis=dict(),
           aspectmode='data')
    fig['layout']['scene1']['xaxis'].update(title='x [m]')
    fig['layout']['scene1']['yaxis'].update(title='y [m]')
    fig['layout']['scene1']['zaxis'].update(title='z [m]')
    fig['layout'].update(title='Start: %.1f, Stop: %.1f' % (start, stop))
    fig['layout']['scene1'].update(scene)
    fig.layout.showlegend = False
    return fig
