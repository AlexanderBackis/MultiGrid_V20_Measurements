#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Misc.py: Helper functions for handling of paths and folders.
"""

from errno import EEXIST
from os import makedirs,path
import numpy as np


# =============================================================================
#                            CREATE DIRECTORY
# =============================================================================

def mkdir_p(my_path):
    """
    Creates a directory, equivalent to using mkdir -p on the command line.

    Args:
        my_path (str): Path to where the new folder should be created.

    Yields:
        A new folder at the requested path.
    """
    try:
        makedirs(my_path)
    except OSError as exc: # Python >2.5
        if exc.errno == EEXIST and path.isdir(my_path):
            pass
        else: raise


# =============================================================================
#                       FIND NEAREST ELEMENT IN ARRAY
# =============================================================================

def find_nearest(array, value):
    """
    Returns the index of the element in 'array' which is closest to 'value'.

    Args:
        array (numpy array): Numpy array with elements
        value (float): Value which we want to find the closest element to in
                       arrray

    Returns:
        idx (int): index of the element in 'array' which is closest to 'value'
    """
    idx = (np.abs(array - value)).argmin()
    return idx

# =============================================================================
#                          GET MEASUREMENT DURATION
# =============================================================================

def get_duration(df):
    times = df.Time.values
    diff = np.diff(times)
    resets = np.where(diff < 0)
    duration_in_TDC_channels = sum(times[resets]) + times[-1]
    duration_in_seconds = duration_in_TDC_channels * 62.5e-9
    return duration_in_seconds


# =============================================================================
#                       APPEND FOLDER AND FILES
# =============================================================================

def append_folder_and_files(folder, files):
    folder_vec = np.array(len(files)*[folder])
    return np.core.defchararray.add(folder_vec, files)

# =============================================================================
#                       ENERGY-WAVELENGTH CONVERSION
# =============================================================================

def meV_to_A(energy):
    return np.sqrt(81.81/energy)


def A_to_meV(wavelength):
    return (81.81/(wavelength ** 2))
