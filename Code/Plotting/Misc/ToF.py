#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ToF.py: Helper functions for handling of paths and folders.
"""

import matplotlib.pyplot as plt

# =============================================================================
#                               ToF - MG
# =============================================================================


def ToF_histogram(df, number_bins):
    fig = plt.figure()
    min_val = 0
    max_val = 16667
    plt.hist(df.ToF * 62.5e-9 * 1e6, bins=number_bins, log=True, color='black',
             zorder=4, range=[min_val, max_val], histtype='step')
    plt.title('ToF')
    plt.xlabel('ToF [$\mu$s]')
    plt.ylabel('Counts')
    plt.yscale('log')
    plt.grid(True, which='major', linestyle='--', zorder=0)
    plt.grid(True, which='minor', linestyle='--', zorder=0)
    return fig
