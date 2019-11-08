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
#                                 Filtering
# =============================================================================

def filter_He3(df, parameters):
    """
    Filters clusters based on preferences set on GUI.

    Args:
        ce (DataFrame): Clustered events
        parameters (dict): Dictionary containing information about which
                           parameters to filter on, and to what extent

    Returns:
        ce_red (DataFrame): DataFrame containing the reduced data according to
                            the specifications set on the GUI.
    """

    df_red = df
    for parameter, (min_val, max_val, filter_on) in parameters.items():
        if filter_on:
            df_red = df_red[(df_red[parameter] >= min_val) &
                            (df_red[parameter] <= max_val)]
    return df_red


# =============================================================================
#                            Helper Functions
# =============================================================================

def get_He3_filter_parameters(window):
    """
    Extracts filtering information from GUI window and saves it in a dictionary.

    Args:
        window (MainWindow): GUI window

    Returns:
        parameters (dict): Dictionary containing information about which
                           parameters to filter on, and to what extent. For each
                           key, the values are:

                           0. min_value (int/float),
                           1. max_value (int/float),
                           2. isFilterOn (bool)

    """

    parameters = {'ADC': [float(window.He3_ADC_min.text()),
                           float(window.He3_ADC_max.text()),
                           window.He3_ADC_filter.isChecked()],
                  'ToF': [float(window.He3_ToF_min.text()) / (8e-9 * 1e6),
                          float(window.He3_ToF_max.text()) / (8e-9 * 1e6),
                          window.ToF_filter.isChecked()]
                  }
    return parameters
