#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Coincidences_3D.py: Helper functions for handling of paths and folders.
"""

# =============================================================================
#                           Coincidence Histogram (3D)
# =============================================================================

def Coincidences_3D_plot(df, data_sets, window):
    # Declare max and min count
    min_count = 0
    max_count = np.inf
    # Perform initial filters
    df = df[(df.wCh != -1) & (df.gCh != -1)]
    df = filter_ce_clusters(window, df)
    if data_sets == 'mvmelst_039.zip':
        df = df[df.Time < 1.5e12]
    # Initiate 'voxel_id -> (x, y, z)'-mapping
    detector_vec = get_detector_mappings()
    # Initiate border lines
    b_traces = initiate_detector_border_lines(detector_vec)
    # Calculate 3D histogram
    H, edges = np.histogramdd(df[['wCh', 'gCh', 'Bus']].values,
                              bins=(80, 40, 9),
                              range=((0, 80), (80, 120), (0, 9))
                              )
    # Insert results into an array
    hist = [[], [], [], []]
    loc = 0
    labels = []
    detector_names = ['ILL', 'ESS_CLB', 'ESS_PA']
    for wCh in range(0, 80):
        for gCh in range(80, 120):
            for bus in range(0, 9):
                detector = detector_vec[bus//3]
                over_min = H[wCh, gCh-80, bus] > min_count
                under_max = H[wCh, gCh-80, bus] <= max_count
                if over_min and under_max:
                    coord = detector[flip_bus(bus % 3), gCh, flip_wire(wCh)]
                    hist[0].append(coord['x'])
                    hist[1].append(coord['y'])
                    hist[2].append(coord['z'])
                    hist[3].append(H[wCh, gCh-80, bus])
                    loc += 1
                    labels.append('Detector: ' + detector_names[(bus//3)]
                                  + '<br>'
                                  + 'Module: ' + str(bus) + '<br>'
                                  + 'WireChannel: ' + str(wCh) + '<br>'
                                  + 'GridChannel: ' + str(gCh) + '<br>'
                                  + 'Counts: ' + str(H[wCh, gCh-80, bus])
                                  )
    # Produce 3D histogram plot
    MG_3D_trace = go.Scatter3d(x=hist[2],
                               y=hist[0],
                               z=hist[1],
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
    fig = py.tools.make_subplots(rows=1, cols=1,
                                 specs=[[{'is_3d': True}]]
                                 )
    # Insert histogram
    fig.append_trace(MG_3D_trace, 1, 1)
    # Insert vessel borders
    for b_trace in b_traces:
        fig.append_trace(b_trace, 1, 1)
    # Assign layout with axis labels, title and camera angle
    a = 1
    camera = dict(up=dict(x=0, y=0, z=1),
                  center=dict(x=0, y=0, z=0),
                  eye=dict(x=-2*a, y=-0.5*a, z=1.3*a)
                  )
    fig['layout']['scene1']['xaxis'].update(title='z [m]')
    fig['layout']['scene1']['yaxis'].update(title='x [m]')
    fig['layout']['scene1']['zaxis'].update(title='y [m]')
    fig['layout'].update(title='Coincidences (3D)<br>' + str(data_sets))
    fig['layout']['scene1']['camera'].update(camera)
    fig.layout.showlegend = False
    # If in plot He3-tubes histogram, return traces, else save HTML and plot
    if data_sets == '':
        return b_traces, hist[0], hist[1], hist[2], np.log10(hist[3])
    else:
        py.offline.plot(fig,
                        filename='../Results/HTML_files/Ce3Dhistogram.html',
                        auto_open=True)
        pio.write_image(fig, '../Results/HTML_files/Ce3Dhistogram.pdf')
