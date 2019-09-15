#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PHS_1D.py: Helper functions for handling of paths and folders.
"""

import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
from matplotlib.colors import LogNorm
from Plotting.HelperFunctions import (stylize, set_figure_properties,
                                      filter_ce_clusters)

# ============================================================================
#                                   PHS (1D)
# ============================================================================


def PHS_1D_plot(events, data_sets, window):
    # Declare parameters
    number_bins = int(window.phsBins.text())
    Channel = window.Channel.value()
    Bus = window.Module.value()
    # Plot
    fig = plt.figure()
    title = ('PHS (1D) - Channel: %d, Bus: %d\nData set(s): %s'
             % (Channel, Bus, data_sets))
    xlabel = 'Counts'
    ylabel = 'Collected charge [ADC channels]'
    fig = stylize(fig, xlabel, ylabel, title, yscale='log', grid=True)
    plt.hist(events[(events.Channel == Channel) & (events.Bus == Bus)].ADC,
             bins=number_bins, range=[0, 4400], histtype='step',
             color='black', zorder=5)
    return fig
