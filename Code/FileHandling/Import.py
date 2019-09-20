#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Import.py: Function which imports data obtained using the Mesytec VMMR 8/16.
"""

import os
import struct
import shutil
import zipfile
import re

from HelperFunctions.PathsAndFolders import mkdir_p

# =============================================================================
#                                UNZIP DATA
# =============================================================================

def unzip_data(zip_source):
    """ Unzips mestyec .zip-files and extracts data-files:

            1. Extracts data in temporary folder
            2. Selects the relevant file and moves it to a temporary location
            3. Removes the temporary fodler where the rest of zipped data is

    Args:
        zip_source (str): Path to '.zip'-file that contains the data

    Returns:
        '.mesytec'-file path (str): Path to the extracted '.mesytec'-file

    """
    dirname = os.path.dirname(__file__)
    zip_temp_folder = os.path.join(dirname, '../zip_temp_folder/')
    mkdir_p(zip_temp_folder)
    file_temp_folder = os.path.join(dirname, '../')
    destination = ''
    with zipfile.ZipFile(zip_source, "r") as zip_ref:
        zip_ref.extractall(zip_temp_folder)
        temp_list = os.listdir(zip_temp_folder)
        source_file = None
        for temp_file in temp_list:
            if temp_file[-8:] == '.mvmelst':
                source_file = temp_file
        source = zip_temp_folder + source_file
        destination = file_temp_folder + source_file
        shutil.move(source, destination)
        shutil.rmtree(zip_temp_folder, ignore_errors=True)
    return destination


# =============================================================================
#                                IMPORT DATA
# =============================================================================

def import_data(file_path, maximum_file_size_in_mb):
    """ Imports mestyec data in three steps:

            1. Reads file as binary and saves data in 'content'
            2. Finds the end of the configuration text, i.e. '}\n}\n' followed
               by 0 to n spaces, then saves everything after this to
               'reduced_content'.
            3. Groups data into 'uint'-words of 4 bytes (32 bits) length

    Args:
        file_path (str): Path to '.mesytec'-file that contains the data

    Returns:
        data (tuple): A tuple where each element is a 32 bit mesytec word

    """

    # Get maximum file size in [bytes]
    ONE_MB_IN_BYTES = (1 << 20)
    maximum_file_size_in_bytes = maximum_file_size_in_mb * ONE_MB_IN_BYTES
    # Assign piece size in [bytes]
    piece_size = 1000 * ONE_MB_IN_BYTES
    # Import data
    with open(file_path, mode='rb') as bin_file:
        # Get first piece of data
        content = bin_file.read(piece_size)
        # Skip configuration text
        match = re.search(b'}\n}\n[ ]*', content)
        start = match.end()
        content = content[start:]
        print(content)
        # Split first piece of data into groups of 4 bytes
        data = struct.unpack('I' * (len(content)//4), content)
        # Repeat for the rest of data
        moreData = True
        imported_data = piece_size
        while moreData and imported_data <= maximum_file_size_in_bytes:
            imported_data += piece_size
            piece = bin_file.read(piece_size)
            if not piece:  # Reached end of file
                moreData = False
            else:
                data += struct.unpack('I' * (len(piece)//4), piece)
    return data
