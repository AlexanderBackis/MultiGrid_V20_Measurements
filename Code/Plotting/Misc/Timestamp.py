#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Timestamp.py: Helper functions for handling of paths and folders.
"""

import matplotlib.pyplot as plt
import numpy as np

# =============================================================================
#                                   Timestamp
# =============================================================================


def timestamp_plot(df):
    # Generate list of event numbers
    event_number = np.arange(0, df.shape[0], 1)
    # Prepare figure
    fig = plt.figure()
    plt.title('Timestamp')
    plt.xlabel('Cluster [index]')
    plt.ylabel('Timestamp [TDC channels]')
    plt.grid(True, which='major', zorder=0)
    plt.grid(True, which='minor', linestyle='--', zorder=0)
    plt.scatter(event_number, df.Time, c=df.wM+df.gM, cmap='jet', zorder=5)
    cbar = plt.colorbar()
    cbar.set_label('Multiplicity summation [wM + gM]')
    plt.tight_layout()
    return fig
