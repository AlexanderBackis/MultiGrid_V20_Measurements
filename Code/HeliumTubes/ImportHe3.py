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
import binascii

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
    """ Imports MCA4 data. This is binary (or ascii?) and is in 64 bit "words".
        The bits in the word, starting from least significant:

        0, 1  : Channel number 0..3 for ADC1..ADC4.
        2     : Flag for pile-up detected
        3     : Indicates scope mode. If it is on, more data words follow.
        4->47 : 44 bit event time. If the RTC option is enabled, it is in units
                of 8 ns, otherwise in units of 1 msec.
        48->63: 16 bit ADC value or in scope mode the length of the waveform
                data minus one in units of 16 bit words. This size value can be
                4095, 8191, 16383 or 32767. The waveform data follow then next.

    Args:
        file_path (str): Path to '.mesytec'-file that contains the data

    Returns:


    """
    test = np.loadtxt(file_path, dtype='str', delimiter='\n')
    for row in test:
        print(binascii.unhexlify(row))

    return data
