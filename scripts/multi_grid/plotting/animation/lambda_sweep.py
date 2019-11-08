#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Animation3D.py: Helper functions for handling of paths and folders.
"""

import plotly as py
import numpy as np
import matplotlib.pyplot as plt
import os
import imageio
import shutil

from multi_grid.plotting.coincidences.coincidences_projections import plot_front, plot_top, plot_side
from multi_grid.plotting.analysis.energy_and_wavelength import energy_plot
from multi_grid.helper_functions.energy_calculation import calculate_energy
from multi_grid.helper_functions.misc import mkdir_p


# =============================================================================
#                                LAMBDA SWEEP
# =============================================================================

def Lambda_Sweep_Animation(ce, number_bins, origin_voxel, bus_start, bus_stop):
    def meV_to_A(energy):
        return np.sqrt(81.81/energy)

    def A_to_meV(wavelength):
        return (81.81/(wavelength ** 2))
    # Filters
    ce = ce[(ce.wCh != -1) & (ce.gCh != -1)]
    # Storage
    dir_name = os.path.dirname(__file__)
    temp_folder = os.path.join(dir_name, '../temp_folder/')
    output_path = os.path.join(dir_name, '../../../../output/lambda_sweep.gif')
    mkdir_p(temp_folder)
    # Calculate energies
    energy = calculate_energy(ce, origin_voxel)
    landa = meV_to_A(energy)
    # Define intervals
    iter_start = 4.3
    iter_stop = 4.75
    step = 0.001
    delimiters = np.arange(iter_start, iter_stop, step)
    # Iteration
    vmin = 1
    vmax = 1e3
    for i, delimiter in enumerate(delimiters):
        start = delimiter
        stop = delimiter + step
        ce_temp = ce[(landa >= start) & (landa <= stop)]
        wChs, gChs, Buses = ce_temp.wCh, ce_temp.gCh, ce_temp.Bus
        fig = plt.figure()
        fig.set_figheight(8)
        fig.set_figwidth(14)
        plt.subplot2grid((2, 3), (0, 0), colspan=1)
        plot_front(wChs, gChs, Buses, bus_start, bus_stop, vmin, vmax)
        plt.subplot2grid((2, 3), (0, 1), colspan=1)
        plot_top(wChs, gChs, Buses, bus_start, bus_stop, vmin, vmax)
        plt.subplot2grid((2, 3), (0, 2), colspan=1)
        plot_side(wChs, gChs, Buses, vmin, vmax)
        plt.subplot2grid((2, 3), (1, 0), colspan=3)
        energy_plot(ce, origin_voxel, number_bins, iter_start, iter_stop)
        plt.axvline(x=(start+stop)/2, color='red', linewidth=2, alpha=0.8, zorder=5)
        plt.tight_layout()
        fig.savefig(temp_folder + '%d.png' % i)
        plt.close(fig)
    # Animate
    images = []
    files = os.listdir(temp_folder)
    files = [file[:-4] for file in files if file[-9:] != '.DS_Store' and file != '.gitignore']
    for file in sorted(files, key=int):
        images.append(imageio.imread(temp_folder + file + '.png'))
    imageio.mimsave(output_path, images)
    shutil.rmtree(temp_folder, ignore_errors=True)
