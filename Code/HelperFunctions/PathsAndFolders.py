#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PathsAndFolders.py: Helper functions for handling of paths and folders.
"""

from errno import EEXIST
from os import makedirs,path

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
