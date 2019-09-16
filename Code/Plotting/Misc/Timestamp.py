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
    # Filter so we only plot every 100:th element, so it is not so slow
    df_temp = df[df.index % 100 == 0]
    # Prepare figure
    fig = plt.figure()
    plt.title('Timestamp (every 100:th cluster)')
    plt.xlabel('Cluster [index]')
    plt.ylabel('Timestamp [TDC channels]')
    plt.grid(True, which='major', zorder=0)
    plt.grid(True, which='minor', linestyle='--', zorder=0)
    plt.scatter(df_temp.index, df_temp.Time, c=df_temp.wM+df_temp.gM,
                cmap='jet', zorder=5)
    cbar = plt.colorbar()
    cbar.set_label('Multiplicity summation [wM + gM]')
    plt.tight_layout()
    return fig
