#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ImportHe3.py: Imports He3 data taken using the MCA4 Multichannel Analyzer
"""

import os
import struct
import shutil
import zipfile
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# =============================================================================
#                               PHS - HELIUM-3
# =============================================================================

def He3_PHS_plot(df, number_bins):
    plt.hist(df['ADC'], histtype='step', color='blue', zorder=5,
             bins=number_bins)
    plt.grid(True, which='major', linestyle='--', zorder=0)
    plt.grid(True, which='minor', linestyle='--', zorder=0)
    plt.xlabel('Charge [ADC Channels]')
    plt.ylabel('Counts')
    plt.title('PHS')

# =============================================================================
#                         CHANNEL HISTOGRAM - HELIUM-3
# =============================================================================

def He3_Ch_plot(df):
    plt.hist(df['Ch'], histtype='step', color='red', zorder=5, bins=20)
    plt.grid(True, which='major', linestyle='--', zorder=0)
    plt.grid(True, which='minor', linestyle='--', zorder=0)
    plt.xlabel('Channel')
    plt.ylabel('Counts')
    plt.title('Channel')

# =============================================================================
#                         TOF HISTOGRAM - HELIUM-3
# =============================================================================

def He3_ToF_plot(df, number_bins):
    plt.hist(df['ToF']*(8e-9)*1e6, histtype='step', color='green', zorder=5,
             bins=number_bins)
    plt.xlabel('ToF [µs]')
    plt.ylabel('Counts')
    plt.yscale('log')
    plt.title('ToF')
    plt.grid(True, which='major', linestyle='--', zorder=0)
    plt.grid(True, which='minor', linestyle='--', zorder=0)
