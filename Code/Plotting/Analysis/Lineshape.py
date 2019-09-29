#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Lineshape.py: Analyses the lineshape using our Figure-of-Merit
"""

from Plotting.Analysis.DeltaE import energy_plot
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks

def analyze_Lineshape(ce_MG, detector_type, origin_voxel, number_bins):
    energy_hist, bin_centers = energy_plot(ce_MG, detector_type, origin_voxel,
                                          number_bins, 1, 10, True)
    # Find peaks
    threshold = 2000
    peaks, *_ = find_peaks(energy_hist, threshold=threshold)
    plt.plot(bin_centers[peaks], energy_hist[peaks], color='red', zorder=5,
             linestyle='', marker='o')
