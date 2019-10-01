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
from binascii import unhexlify, hexlify

# =============================================================================
#                             UNZIP HELIUM-3 DATA
# =============================================================================

def unzip_He3_data(zip_source):
    """ Unzips MCA4 '.zip'-file, extracts data into '.lst'-file, and returns
        the path to the file.

    Args:
        zip_source (str): Path to '.zip'-file that contains the data

    Returns:
        file_path (str): Path to the extracted '.lst'-file

    """
    dirname = os.path.dirname(__file__)
    destination = os.path.join(dirname, '../../')
    with zipfile.ZipFile(zip_source, "r") as zip_ref:
        zip_ref.extractall(destination)
    file_name = zip_source.rsplit('/', 1)[-1][:-3] + 'lst'
    file_path = destination + file_name
    return file_path


# =============================================================================
#                           IMPORT HELIUM-3 DATA
# =============================================================================

def import_He3_data(file_path):
    print(file_path)
    """ Imports MCA4 data. This is binary (or ascii?) and is in 64 bit "words".
        The bits in the word, starting from least significant:

        0, 1  : Channel number 0..3 for ADC1..ADC4.
        2     : Flag for pile-up detected
        3     : Indicates scope mode. If it is on, more data words follow.
        4->47 : 44 bit event time, in units of 8 [ns]
        48->63: 16 bit ADC value, size can be 4095, 8191, 16383 or 32767.

    Args:
        file_path (str): Path to '.mesytec'-file that contains the data

    Returns:
    """
    # 1->4 = amplitude
    # 5->15 = Time
    # 16 Channel

    # Masks
    ChannelMask = 0x0000000000000003 
    PileUpMask  = 0x000000000000000C
    TimeMask    = 0x0000FFFFFFFFFFF0
    ADCMask     = 0xFFFF000000000000
    BreakMask   = 0xFFF0000000000000
    # Bit shifts
    ChannelShift = 0
    TimeShift    = 4
    ADCShift     = 48
    # Import data
    data = np.loadtxt(file_path, dtype='str', delimiter='\n')
    start_idx = np.where(data == '[DATA]')[0][0]
    size = len(data)
    # Declare dictionary to store data
    He3_dict = {'Ch':  np.empty([size], dtype=int),
                'ToF': np.empty([size], dtype=int),
                'ADC': np.empty([size], dtype=int)}
    count = 0
    # Extracts information from data
    for i, row in enumerate(data[start_idx+1:]):
        # Convert ASCII encoded HEX to int (shouldn't it be uint?)
        word = int(row, 16)
        # Check if we should save data
        if (word & BreakMask) != 0:
            # Extract values using masks
            He3_dict['Ch'][count] = (word & ChannelMask)
            He3_dict['ToF'][count] = (word & TimeMask) >> TimeShift
            He3_dict['ADC'][count] = (word & ADCMask) >> ADCShift
            count += 1
    # Only save the events, cut unused rows
    for key in He3_dict:
        He3_dict[key] = He3_dict[key][0:count]
    He3_df = pd.DataFrame(He3_dict)
    return He3_df
