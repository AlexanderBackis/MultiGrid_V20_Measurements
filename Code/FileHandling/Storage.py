#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Storage.py: Saves and Loads clustered data via HDF5 files.
"""

import pandas as pd


# =============================================================================
#                                LOAD DATA
# =============================================================================

def load_data(path):
    """
    Loads clustered data from specified HDF5-path.

    Args:
        path (str): Path to HDF5 containing the saved clusters and events

    Returns:
        ce (DataFrame): Clusters
        e (DataFrame): Events
        data_sets (str): Data used for clustering
        adc_threshold (int): ADC threshold used in clustering procedure
        ILL_buses (list): List containing all ILL buses
    """
    ce = pd.read_hdf(path, 'ce')
    e = pd.read_hdf(path, 'e')
    data_sets = pd.read_hdf(path, 'data_sets')['data_sets'].iloc[0]
    return ce, e, data_sets

# =============================================================================
#                                SAVE DATA
# =============================================================================

def save_data(path, ce, e, data_sets):
    """
    Saves clustered data to specified HDF5-path.

    Args:
        path (str): Path to HDF5 containing the saved clusters and events
        ce (DataFrame): Clusters
        e (DataFrame): Events
        data_sets (str): Data used for clustering

    Yields:
        Clustered data is saved to specified path.
    """
    # Convert values to DataFrame for easy storage
    data_sets_df = pd.DataFrame({'data_sets': [data_sets]})
    # Export to HDF5, use highest compression level
    ce.to_hdf(path, 'ce', complevel=9)
    e.to_hdf(path, 'e', complevel=9)
    data_sets_df.to_hdf(path, 'data_sets', complevel=9)
