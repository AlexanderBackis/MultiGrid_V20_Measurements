#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Storage.py:
"""

import pandas as pd


# =============================================================================
#                                LOAD DATA
# =============================================================================

def load_data(path):
    ce = pd.read_hdf(path, 'ce')
    e = pd.read_hdf(path, 'e')
    data_sets = pd.read_hdf(path, 'data_sets')['data_sets'].iloc[0]
    adc_threshold = pd.read_hdf(path, 'adc_threshold')['adc_threshold'].iloc[0]
    ILL_buses = pd.read_hdf(path, 'ILL_buses')['ILL_buses'].values[0]
    return ce, e, data_sets, adc_threshold, ILL_buses

# =============================================================================
#                                SAVE DATA
# =============================================================================

def save_data(path, ce, e, data_sets, adc_threshold, ILL_buses):
    # Convert values to DataFrame for easy storage
    data_sets_df = pd.DataFrame({'data_sets': [data_sets]})
    adc_threshold_df = pd.DataFrame({'adc_threshold': [adc_threshold]})
    ILL_buses_df = pd.DataFrame({'ILL_buses': [ILL_buses]})
    # Export to HDF5, use highest compression level
    ce.to_hdf(path, 'ce', complevel=9)
    e.to_hdf(path, 'e', complevel=9)
    data_sets_df.to_hdf(path, 'data_sets', complevel=9)
    adc_threshold_df.to_hdf(path, 'adc_threshold', complevel=9)
    ILL_buses_df.to_hdf(path, 'ILL_buses', complevel=9)
