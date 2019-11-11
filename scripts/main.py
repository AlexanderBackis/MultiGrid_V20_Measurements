#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" main.py: Module containing the GUI used for data analysis.
"""

# Standard library
import os
import sys
import matplotlib.pyplot as plt
import pandas as pd
# QT
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5 import uic
# Data Analysis
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
# File handling
from multi_grid.file_handling.read import unzip_data, import_data
from multi_grid.file_handling.cluster import cluster_data
from multi_grid.file_handling.storage import save_data, load_data
from multi_grid.file_handling.data_preparation import prepare_data, plot_all_peaks
# He-3 tube
from helium_tube.file_handling_he3 import unzip_He3_data, import_He3_data, save_He3_data, load_He3_data
from helium_tube.plotting_he3 import He3_PHS_plot, He3_ToF_plot, He3_Ch_plot, energy_plot_He3, he3_pileup_plot
from helium_tube.filtering_he3 import get_He3_filter_parameters, filter_He3
from helium_tube.energy_he3 import calculate_He3_energy
# Helper functions
from multi_grid.helper_functions.mapping import create_mapping
from multi_grid.helper_functions.filtering import filter_clusters, get_filter_parameters
from multi_grid.helper_functions.energy_calculation import calculate_energy
from multi_grid.helper_functions.misc import append_folder_and_files, meV_to_A, A_to_meV, get_duration
# PHS
from multi_grid.plotting.phs.phs_1d import PHS_1D_plot
from multi_grid.plotting.phs.phs_2d import PHS_2D_plot
from multi_grid.plotting.phs.phs_wires_vs_grids import PHS_wires_vs_grids_plot
# Coincidences
from multi_grid.plotting.coincidences.coincidences_2d import coincidences_2D_plot
from multi_grid.plotting.coincidences.coincidences_3d import coincidences_3D_plot
from multi_grid.plotting.coincidences.coincidences_projections import coincidences_projections_plot
# Misc
from multi_grid.plotting.misc.multiplicity import multiplicity_plot
from multi_grid.plotting.misc.tof import ToF_histogram
from multi_grid.plotting.misc.timestamp import timestamp_plot
# Analysis
from multi_grid.plotting.analysis.energy_and_wavelength import energy_plot
from multi_grid.plotting.analysis.count_rate import calculate_count_rate
from multi_grid.plotting.analysis.efficiency import plot_efficiency
from multi_grid.plotting.analysis.energy_resolution import calculate_energy_resolution
from multi_grid.plotting.analysis.lineshape import plot_FoM
from multi_grid.plotting.analysis.layers import (investigate_layers_ToF,
                                                 investigate_layers_FWHM,
                                                 investigate_layers_delta_ToF,
                                                 investigate_layers_counts)
# Animation
from multi_grid.plotting.animation.animation_3d import Animation_3D_plot
from multi_grid.plotting.animation.lambda_sweep import Lambda_Sweep_Animation
from multi_grid.plotting.animation.time_sweep import Time_Sweep_Animation
from multi_grid.plotting.animation.tof_sweep import ToF_Sweep_Animation




# =============================================================================
#                                   Windows
# =============================================================================

class MainWindow(QMainWindow):
    def __init__(self, app, parent=None):
        super(MainWindow, self).__init__(parent)
        dir_name = os.path.dirname(__file__)
        title_screen_path = os.path.join(dir_name, '../windows/mainwindow.ui')
        self.ui = uic.loadUi(title_screen_path, self)
        self.app = app
        # Clustering attributes
        self.data_sets = ''
        self.maximum_file_size_in_mb = 3000
        # Cluster properties
        self.measurement_time = 0
        self.ce = pd.DataFrame()
        self.e = pd.DataFrame()
        self.fill_MG_information_window()
        # He3 attributes
        self.He3_data_sets = ''
        self.He3_counts = 0
        self.He3_df = pd.DataFrame()
        self.fill_He3_information_window()
        self.show()
        self.refresh_window()
        # Beam Monitor attributes
        self.BM_counts_dict = {}
        #self.import_beam_monitor_data()
        self.initialize_bm_dict()

    # =========================================================================
    # File handling
    # =========================================================================

    def cluster_action(self):
        zip_paths = QFileDialog.getOpenFileNames(self, "", "../data")[0]
        if len(zip_paths) > 0:
            # Import
            data = ()
            data_sets_temp = '<br/>'
            for i, zip_path in enumerate(zip_paths):
                file_name = zip_path.rsplit('/', 1)[-1]
                data_sets_temp += file_name + '<br/>'
                file_path = unzip_data(zip_path)
                data += import_data(file_path, self.maximum_file_size_in_mb)
                os.remove(file_path)
                print('%d/%d' % (i+1, len(zip_paths)))
            # Cluster
            clusters, events = cluster_data(data)
            # Write or append
            if self.write.isChecked():
                self.ce = clusters
                self.e = events
                self.data_sets = data_sets_temp
            else:
                self.ce = self.ce.append(clusters)
                self.e = self.e.append(events)
                self.data_sets += data_sets_temp
                # Reset index
                self.ce.reset_index(drop=True, inplace=True)
                self.e.reset_index(drop=True, inplace=True)
            # Update window
            self.measurement_time = get_duration(self.ce)
            self.fill_MG_information_window()
            self.refresh_window()
            # Print statements for debugging purposes
            print(self.ce)
            print(self.e)

    def save_action(self):
        path = QFileDialog.getSaveFileName()[0]
        if path != '':
            save_data(path, self.ce, self.e, self.data_sets)

    def load_action(self):
        path = QFileDialog.getOpenFileName(self, "", "../data")[0]
        if path != '':
            clusters, events, data_sets_temp = load_data(path)
            # Write or append
            if self.write.isChecked():
                self.ce = clusters
                self.e = events
                self.data_sets = data_sets_temp
            else:
                self.ce = self.ce.append(clusters)
                self.e = self.e.append(events)
                self.data_sets += data_sets_temp
                # Reset index
                self.ce.reset_index(drop=True, inplace=True)
                self.e.reset_index(drop=True, inplace=True)
            # Update window
            self.measurement_time = get_duration(self.ce)
            self.fill_MG_information_window()
            self.refresh_window()
        print(self.ce)


    # =========================================================================
    # Plotting
    # =========================================================================

    # ==== PHS ==== #

    def PHS_1D_action(self):
        if self.data_sets != '':
            filter_parameters = get_filter_parameters(self)
            ce_filtered = filter_clusters(self.ce, filter_parameters)
            number_bins = int(self.phsBins.text())
            fig = plt.figure()
            fig.set_figheight(5)
            fig.set_figwidth(10)
            PHS_1D_plot(ce_filtered, number_bins)
            fig.show()


    def PHS_2D_action(self):
        if self.data_sets != '':
            bus_start = self.module_min.value()
            bus_stop = self.module_max.value()
            fig = PHS_2D_plot(self.e, bus_start, bus_stop)
            fig.show()

    def PHS_wires_vs_grids_action(self):
        if (self.data_sets != ''):
            filter_parameters = get_filter_parameters(self)
            ce_filtered = filter_clusters(self.ce, filter_parameters)
            bus_start = self.module_min.value()
            bus_stop = self.module_max.value()
            fig = PHS_wires_vs_grids_plot(ce_filtered, bus_start, bus_stop)
            fig.show()

    def PHS_comparison_action(self):
        paths = QFileDialog.getOpenFileNames(self, "", "../data")[0]
        if len(paths) > 0:
            # Declare parameters
            filter_parameters = get_filter_parameters(self)
            number_bins = int(self.phsBins.text())
            ylabel = '(Normalized to beam monitor)'
            intervals = [[0, 4095], [0, 10000]]
            # Extract histogram data
            fig = plt.figure()
            fig.set_figheight(5)
            fig.set_figwidth(10)
            labels = ['(Beam)', '(Background)']
            hists = []
            bins_vec = []
            errors = []
            norms = []
            # Plot data
            for path, label in zip(paths, labels):
                # Import data
                data_set = path.rsplit('/', 1)[-1][:-3] + '.zip'
                ce = pd.read_hdf(path, 'ce')
                ce_filtered = filter_clusters(ce, filter_parameters)
                # Extract normalization
                norm = 1/self.BM_counts_dict[data_set]
                bins, hist = PHS_1D_plot(ce_filtered, number_bins, label, norm,
                                         ylabel, intervals)
                hists.append(np.array(hist))
                errors.append(np.array(hist)/norm)
                norms.append(norm)
                bins_vec.append(bins)
            # Plot difference
            for i in range(0, 2):
                plt.subplot(1, 2, i+1)
                plt.plot(bins_vec[0][i], hists[0][i] - hists[1][i],
                         #np.sqrt((errors[0][i]*norms[0]) ** 2 + (errors[1][i]*norms[1]) ** 2),
                         #fmt='.-', capsize=5,
                         zorder=3, label='Difference')
                plt.legend(loc=1)
            fig.show()

    def PHS_comparison_region_action(self):
        if (self.data_sets != ''):
            # Declare parameters
            number_bins = int(self.phsBins.text())
            filter_parameters = get_filter_parameters(self)
            df = filter_clusters(self.ce, filter_parameters)
            # Extract data from different regions
            df_cross = df[(((df.Bus * 4) + df.wCh//20) >= 4) &
                          (((df.Bus * 4) + df.wCh//20) <= 7) &
                          (df.gCh >= 98) &
                          (df.gCh <= 101)]
            df_neutron = df[(((df.Bus * 4) + df.wCh//20) == 6) &
                          (df.gCh >= 87) &
                          (df.gCh <= 89)]
            norm = 1/self.measurement_time
            fig = plt.figure()
            fig.set_figheight(5)
            fig.set_figwidth(10)
            PHS_1D_plot(df_cross, number_bins, '(Cross talk)', norm)
            PHS_1D_plot(df_neutron, number_bins, '(Neutrons)', norm)
            plt.legend()
            fig.show()

    # ==== Coincidences ==== #

    def Coincidences_2D_action(self):
        if self.data_sets != '':
            filter_parameters = get_filter_parameters(self)
            ce_filtered = filter_clusters(self.ce, filter_parameters)
            bus_start = self.module_min.value()
            bus_stop = self.module_max.value()
            fig, histograms = coincidences_2D_plot(ce_filtered, self.measurement_time, bus_start, bus_stop)
            # Export histograms to text
            dir_name = os.path.dirname(__file__)
            output_path = os.path.join(dir_name, '../output/')
            for bus, histogram in enumerate(histograms):
                path = output_path + '2D_Coincidences_Bus_%d.txt' % bus
                np.savetxt(path, histogram, fmt="%d", delimiter=",")
            # Plot figure
            fig.show()


    def Coincidences_3D_action(self):
        if self.data_sets != '':
            filter_parameters = get_filter_parameters(self)
            ce_filtered = filter_clusters(self.ce, filter_parameters)
            origin_voxel = [int(self.bus_origin.text()),
                            int(self.gCh_origin.text()),
                            int(self.wCh_origin.text())]
            coincidences_3D_plot(ce_filtered, origin_voxel)

    def Coincidences_Projections_action(self):
        if self.data_sets != '':
            filter_parameters = get_filter_parameters(self)
            ce_filtered = filter_clusters(self.ce, filter_parameters)
            bus_start = self.module_min.value()
            bus_stop = self.module_max.value()
            # Get beam monitor data
            file_name = self.data_sets[5:-5]
            norm = 1/get_duration(ce_filtered)#self.BM_counts_dict[file_name]
            fig, histograms = coincidences_projections_plot(ce_filtered, bus_start, bus_stop, norm)
            # Export histograms to text
            dir_name = os.path.dirname(__file__)
            output_path = os.path.join(dir_name, '../output/')
            file_names = ['Front', 'Top', 'Side']
            for file_name, histogram in zip(file_names, histograms):
                path = output_path + '2D_Coincidences_Projections_%s.txt' % file_name
                np.savetxt(path, histogram, fmt="%d", delimiter=",")
            # Plot figure
            fig.show()


    def CE_2D_comparison_action(self):
        paths = QFileDialog.getOpenFileNames(self, "", "../data")[0]
        if len(paths) > 0:
            # Declare parameters
            filter_parameters = get_filter_parameters(self)
            front_hists = []
            top_hists = []
            side_hists = []
            # Extract histogram data
            for i, path in enumerate(paths):
                # Import data
                data_set = path.rsplit('/', 1)[-1][:-3] + '.zip'
                ce = pd.read_hdf(path, 'ce')
                ce_filtered = filter_clusters(ce, filter_parameters)
                # Extract histograms
                bus_start = self.module_min.value()
                bus_stop = self.module_max.value()
                norm = 1/self.BM_counts_dict[data_set]
                __, histograms = coincidences_projections_plot(ce_filtered, bus_start, bus_stop, norm)
                front_hists.append(np.transpose(histograms[0]))
                top_hists.append(np.transpose(histograms[1]))
                side_hists.append(np.transpose(histograms[2]))
            # Take difference between histograms
            projections_hists = [front_hists, top_hists, side_hists]
            diffs = []
            for hists in projections_hists:
                diffs.append(hists[0] - hists[1])
            # Define figure properties
            labels_vec = [['Front view', 'Row', 'Grid'],
                          ['Top view', 'Row', 'Layer'],
                          ['Side view', 'Layer', 'Grid']]
            limits_vec = [[3e-6, 2e-3], [3e-4, 4e-3], [3e-5, 2e-3]]
            # Plot difference
            fig = plt.figure()
            fig.set_figheight(5)
            fig.set_figwidth(17)
            for i, (diff, labels, limits) in enumerate(zip(diffs, labels_vec, limits_vec)):
                vmin, vmax = limits
                plt.subplot(1, 3, 1+i)
                plt.imshow(diff, cmap='jet', norm=LogNorm(), origin='lower',
                           interpolation='nearest', aspect='auto',
                           vmin=vmin, vmax=vmax)
                # Stylize plot
                title, xlabel, ylabel = labels
                plt.title(title)
                plt.xlabel(xlabel)
                plt.ylabel(ylabel)
                cbar = plt.colorbar()
                cbar.set_label('Counts (Normalized to beam monitor)')
            plt.tight_layout()
            fig.show()








    # ==== Misc ==== #

    def Multiplicity_action(self):
        if (self.data_sets != ''):
            filter_parameters = get_filter_parameters(self)
            ce_filtered = filter_clusters(self.ce, filter_parameters)
            bus_start = self.module_min.value()
            bus_stop = self.module_max.value()
            fig = multiplicity_plot(ce_filtered, bus_start, bus_stop)
            fig.show()

    def ToF_action(self):
        if self.data_sets != '':
            # Filter data
            filter_parameters = get_filter_parameters(self)
            ce_filtered = filter_clusters(self.ce, filter_parameters)
            number_bins = int(self.tofBins.text())
            # Plot data
            fig = plt.figure()
            hist, bins = ToF_histogram(ce_filtered, number_bins, range=[0, 71429])
            # Save histogram in ASCII-format
            dir_name = os.path.dirname(__file__)
            output_path = os.path.join(dir_name, '../output/ToF_%s.txt' % self.data_sets[5:-5])
            np.savetxt(output_path,
                       np.transpose(np.array([bins, hist])),
                       delimiter=",",
                       header='ToF (µs), Counts')
            fig.show()

    def ToF_Overlay_action(self):
        paths = QFileDialog.getOpenFileNames(self, "", "../data")[0]
        if len(paths) > 0:
            # Declare parameters
            filter_parameters = get_filter_parameters(self)
            number_bins = int(self.tofBins.text())
            labels = ['Beam', 'Background']
            # Plot
            fig = plt.figure()
            hists = []
            bin_centers_vec = []
            range = [0, 71e3]
            for path, label in zip(paths, labels):
                # Get data
                data_set = path.rsplit('/', 1)[-1][:-3] + '.zip'
                ce = pd.read_hdf(path, 'ce')
                ce_filtered = filter_clusters(ce, filter_parameters)
                # Get normalization
                norm = 1/self.BM_counts_dict[data_set]
                # Plot
                hist, bins = ToF_histogram(ce_filtered, number_bins, label, norm, range)
                hists.append(hist)
                bin_centers_vec.append(bins)
            plt.plot(bin_centers_vec[0], hists[0]-hists[1], label='Difference', zorder=5)
            plt.title('ToF')
            plt.xlabel('ToF [$\mu$s]')
            plt.ylabel('Counts (Normalized by beam monitor)')
            plt.legend()
            plt.grid(True, which='major', linestyle='--', zorder=0)
            plt.grid(True, which='minor', linestyle='--', zorder=0)
            fig.show()

    def Timestamp_action(self):
        if (self.data_sets != ''):
            filter_parameters = get_filter_parameters(self)
            ce_filtered = filter_clusters(self.ce, filter_parameters)
            number_bins = int(self.tofBins.text())
            unit = 'hours'
            fig = plt.figure()
            bus_start = self.module_min.value()
            bus_stop = self.module_max.value()
            for Bus in np.arange(bus_start, bus_stop+1, 1):
                df_temp = ce_filtered[ce_filtered.Bus == Bus]
                if df_temp.shape[0] > 0:
                    Time = (df_temp.Time*62.5e-9)/(60 ** 2)
                    label = 'Bus %d' % Bus
                    timestamp_plot(Time, number_bins, unit, label)
            plt.legend(loc=2)
            fig.show()

    # ==== Analysis ==== #

    def Energy_action(self):
        if (self.data_sets != ''):
            filter_parameters = get_filter_parameters(self)
            ce_filtered = filter_clusters(self.ce, filter_parameters)
            number_bins = int(self.dE_bins.text())
            origin_voxel = [int(self.bus_origin.text()),
                            int(self.gCh_origin.text()),
                            int(self.wCh_origin.text())]
            fig = plt.figure()
            start = 0.8
            stop = 10
            plot_energy = True
            energy_plot(ce_filtered, origin_voxel, number_bins, start, stop, plot_energy)
            fig.show()

    def Wavelength_action(self):
        if (self.data_sets != ''):
            filter_parameters = get_filter_parameters(self)
            ce_filtered = filter_clusters(self.ce, filter_parameters)
            number_bins = int(self.dE_bins.text())
            origin_voxel = [int(self.bus_origin.text()),
                            int(self.gCh_origin.text()),
                            int(self.wCh_origin.text())]
            fig = plt.figure()
            start = 0
            stop = 10
            plot_energy = False
            energy_plot(ce_filtered, origin_voxel, number_bins, start, stop, plot_energy)
            fig.show()

    def Wavelength_overlay_action(self):
        paths = QFileDialog.getOpenFileNames(self, "", "../data")[0]
        number_bins = int(self.dE_bins.text())
        origin_voxel = [int(self.bus_origin.text()),
                        int(self.gCh_origin.text()),
                        int(self.wCh_origin.text())]
        MG_filter_parameters = get_filter_parameters(self)
        He3_filter_parameters = get_He3_filter_parameters(self)
        fig = plt.figure()
        # Multi-Grid
        for i, path in enumerate(paths):
            label = path.rsplit('/', 1)[-1]
            ce = pd.read_hdf(path, 'ce')
            ce_filtered = filter_clusters(ce, MG_filter_parameters)
            duration = get_duration(ce)
            norm = 1/duration
            energy = calculate_energy(ce_filtered, origin_voxel)
            plt.hist(meV_to_A(energy), bins=number_bins, range=[0, 10], zorder=5,
                     histtype='step', label=label, weights=norm*np.ones(len(energy)))
            print('Progress: %d/%d' % (i+1, len(paths)))
        # He-3
        He3_df_red = filter_He3(self.He3_df, He3_filter_parameters)
        energy_He3 = calculate_He3_energy(He3_df_red)
        norm_He3 = 1/54304
        plt.hist(meV_to_A(energy_He3), bins=number_bins, range=[0, 10], zorder=5,
                 histtype='step',
                 label='2019_09_HZB_He3InBeam54304s_overnight.lst',
                 weights=norm*np.ones(len(energy_He3)))
        plt.grid(True, which='major', linestyle='--', zorder=0)
        plt.grid(True, which='minor', linestyle='--', zorder=0)
        plt.ylabel('Counts (Normalized to duration)')
        plt.xlabel('Wavelength [Å]')
        plt.title('Wavelength Distribution')
        plt.legend(loc=1)
        fig.show()

    def Count_Rate_action(self):
        if (self.data_sets != ''):
            # Declare parameters
            number_bins = int(self.tofBins.text())
            time_offset = (0.6e-3)
            period_time = (1/14)
            filter_parameters = get_filter_parameters(self)
            ce_filtered = filter_clusters(self.ce, filter_parameters)
            ToF_values = (ce_filtered.ToF * 62.5e-9 + time_offset) % period_time
            fig = plt.figure()
            count_rate = calculate_count_rate(ToF_values,
                                              self.measurement_time,
                                              number_bins)
            fig.show()
            print('Count rate: %.1f [Hz]' % count_rate)

    def Efficiency_action(self):
        dirname = os.path.dirname(__file__)
        He3_efficiency_path = os.path.join(dirname, '../tables/He3_efficiency.txt')
        MG_efficiency_path = os.path.join(dirname, '../tables/MG_efficiency.txt')
        He3_efficiency = np.loadtxt(He3_efficiency_path, delimiter=",", unpack=True)
        MG_efficiency = np.loadtxt(MG_efficiency_path, delimiter=",", unpack=True)[[0, 2]]
        fig = plt.figure()
        start = 0.01
        end = 12
        fig.suptitle('Efficiency, Multi-Grid and He-3')
        fig.set_figheight(5)
        fig.set_figwidth(10)
        plt.subplot(1, 2, 1)
        plt.plot(He3_efficiency[0], He3_efficiency[1], color='red',
                 label='He-3 (Average, Metal Container)', zorder=5)
        #plt.plot(MG_efficiency[0], MG_efficiency[1], color='blue',
        #         label='MG (90° incident angle)', zorder=5)
        plt.title('Wavelength')
        plt.xlabel('Wavelength (Å)')
        plt.ylabel('Efficiency')
        plt.grid(True, which='major', linestyle='--', zorder=0)
        plt.grid(True, which='minor', linestyle='--', zorder=0)
        plt.xlim(start, end)
        plt.legend()
        plt.subplot(1, 2, 2)
        plt.plot(A_to_meV(He3_efficiency[0]), He3_efficiency[1], color='red',
                 label='He-3 (Average, Metal Container)', zorder=5)
        #plt.plot(A_to_meV(MG_efficiency[0]), MG_efficiency[1], color='blue',
        #         label='MG (90° incident angle)', zorder=5)
        plt.title('Energy')
        plt.xlabel('Energy (meV)')
        plt.ylabel('Efficiency')
        plt.grid(True, which='major', linestyle='--', zorder=0)
        plt.grid(True, which='minor', linestyle='--', zorder=0)
        plt.xlim(A_to_meV(end), A_to_meV(start))
        plt.xscale('log')
        plt.legend()
        fig.show()


    def Energy_Resolution_action(self):
        if (self.data_sets != ''):
            filter_parameters = get_filter_parameters(self)
            ce_filtered = filter_clusters(self.ce, filter_parameters)
            dE_values = calculate_energy_transfer(ce_filtered, Ei)
            # NEED A WAY TO GET HE3 VALUES
            FWHM = calculate_energy_resolution(dE_values, self.Ei, filter_parameters)
            print('FWHM: %.2f' % FWHM)

    def Lineshape_action(self):
        origin_voxel = [int(self.bus_origin.text()),
                        int(self.gCh_origin.text()),
                        int(self.wCh_origin.text())]
        MG_filter_parameters = get_filter_parameters(self)
        He3_filter_parameters = get_He3_filter_parameters(self)
        analyze_all_lineshapes(origin_voxel, MG_filter_parameters, He3_filter_parameters)

    def Layers_action(self):
        # Extract origin voxel
        origin_voxel = [int(self.bus_origin.text()),
                        int(self.gCh_origin.text()),
                        int(self.wCh_origin.text())]
        # Get MG data
        MG_parameters = get_filter_parameters(self)
        df_MG = filter_clusters(self.ce, MG_parameters)
        # Get He3 data
        #He3_parameters = get_He3_filter_parameters(self)
        #df_He3 = filter_He3(self.He3_df, He3_parameters)
        # Investigate ToF spread
        #investigate_layers_FWHM(df_MG, df_He3, origin_voxel)
        #investigate_layers_delta_ToF(df_MG, df_He3, origin_voxel)
        duration = get_duration(df_MG)
        investigate_layers_counts(df_MG, duration)

    def full_analysis_action(self):
        # Prepare filter parameters
        origin_voxel = [int(self.bus_origin.text()),
                        int(self.gCh_origin.text()),
                        int(self.wCh_origin.text())]
        MG_filter_parameters = get_filter_parameters(self)
        He3_filter_parameters = get_He3_filter_parameters(self)
        # Define colors and normalization
        colors = {'MG_Coated': 'blue', 'MG_Non_Coated': 'green', 'He3': 'red'}
        monitor_norm_coated = 1/11411036
        monitor_norm_non_coated = 1/9020907
        monitor_norm_He3 = 1/10723199
        # Prepare data
        full_data = prepare_data(origin_voxel, MG_filter_parameters, He3_filter_parameters)
        MG_coated_data, MG_non_coated_data, He3_data = full_data[0], full_data[1], full_data[4]
        MG_coated_background, MG_non_coated_background, He3_background = full_data[2], full_data[3], full_data[5]
        # Plot all peaks
        Coated_values = plot_all_peaks(MG_coated_data, 'MG_Coated', colors['MG_Coated'], 28.413)
        NonCoated_values =  plot_all_peaks(MG_non_coated_data, 'MG_Non_Coated', colors['MG_Non_Coated'], 28.413+1.5e-3)
        He3_values = plot_all_peaks(He3_data, 'He3', colors['He3'], 28.239+3e-3)
        # Extract key values
        energies_Coated, FoM_Coated, FoM_err_Coated, peak_areas_Coated, peak_err_Coated = Coated_values
        energies_NonCoated, FoM_NonCoated, FoM_err_NonCoated, peak_areas_NonCoated, peak_err_NonCoated = NonCoated_values
        energies_He3, FoM_He3, FoM_err_He3, peak_areas_He3, peak_err_He3 = He3_values
        # Store all values in vectors
        energies = [energies_Coated, energies_NonCoated, energies_He3]
        labels = ['MG_Coated', 'MG_Non_Coated', 'He3']
        FoMs = [FoM_Coated, FoM_NonCoated, FoM_He3]
        FoM_errors = [FoM_err_Coated, FoM_err_NonCoated, FoM_err_He3]
        # Plot efficiency
        fig = plt.figure()
        fig.set_figheight(5)
        fig.set_figwidth(15)
        plot_efficiency(np.array(energies_He3), np.array(energies_NonCoated),
                        np.array(peak_areas_He3), np.array(peak_areas_NonCoated),
                        np.array(peak_err_He3), np.array(peak_err_NonCoated),
                        monitor_norm_He3, monitor_norm_non_coated)
        fig.show()
        # Plot FoM
        fig = plt.figure()
        for energy, FoM, error, label in zip(energies, FoMs, FoM_errors, labels):
            plot_FoM(energy, FoM, error, label, colors[label])
        plt.legend()
        fig.show()


    # ==== Animation ==== #

    def Animation_3D_action(self):
        filter_parameters = get_filter_parameters(self)
        ce_filtered = filter_clusters(self.ce, filter_parameters)
        origin_voxel = [int(self.bus_origin.text()),
                        int(self.gCh_origin.text()),
                        int(self.wCh_origin.text())]
        Animation_3D_plot(ce_filtered, origin_voxel)

    def lambda_sweep_action(self):
        filter_parameters = get_filter_parameters(self)
        ce_filtered = filter_clusters(self.ce, filter_parameters)
        origin_voxel = [int(self.bus_origin.text()),
                        int(self.gCh_origin.text()),
                        int(self.wCh_origin.text())]
        number_bins = int(self.dE_bins.text())
        bus_start = self.module_min.value()
        bus_stop = self.module_max.value()
        Lambda_Sweep_Animation(ce_filtered, number_bins,
                               origin_voxel, bus_start, bus_stop)

    def time_sweep_action(self):
        filter_parameters = get_filter_parameters(self)
        ce_filtered = filter_clusters(self.ce, filter_parameters)
        origin_voxel = [int(self.bus_origin.text()),
                        int(self.gCh_origin.text()),
                        int(self.wCh_origin.text())]
        number_bins = int(self.tofBins.text())
        bus_start = self.module_min.value()
        bus_stop = self.module_max.value()
        Time_Sweep_Animation(ce_filtered, number_bins,
                             origin_voxel, bus_start, bus_stop)

    def ToF_sweep_action(self):
        filter_parameters = get_filter_parameters(self)
        ce_filtered = filter_clusters(self.ce, filter_parameters)
        number_bins = int(self.tofBins.text())
        bus_start = self.module_min.value()
        bus_stop = self.module_max.value()
        ToF_Sweep_Animation(ce_filtered, number_bins, bus_start, bus_stop)





    # =========================================================================
    # He-3 tubes
    # =========================================================================

    def Import_He3_action(self):
        file_path = QFileDialog.getOpenFileName(self, "", "../data")[0]
        if file_path != '':
            self.He3_df = import_He3_data(file_path)
            self.He3_data_sets = '<br/>' + file_path.rsplit('/', 1)[-1]
            self.He3_counts = self.He3_df.shape[0]
            self.fill_He3_information_window()
            self.show()
            self.refresh_window()

    def Save_He3_action(self):
        path = QFileDialog.getSaveFileName()[0]
        if path != '':
            save_He3_data(path, self.He3_df, self.He3_data_sets)

    def Load_He3_action(self):
        path = QFileDialog.getOpenFileName(self, "", "../data")[0]
        if path != '':
             self.He3_df, self.He3_data_sets = load_He3_data(path)
        self.fill_He3_information_window()
        self.show()
        self.refresh_window()

    def He3_PHS_action(self):
        number_bins = int(self.phsBins.text())
        parameters = get_He3_filter_parameters(self)
        df_red = filter_He3(self.He3_df, parameters)
        fig = plt.figure()
        He3_PHS_plot(df_red, number_bins)
        fig.show()

    def He3_ToF_action(self):
        number_bins = int(self.tofBins.text())
        parameters = get_He3_filter_parameters(self)
        df_red = filter_He3(self.He3_df, parameters)
        fig = plt.figure()
        hist, bins = He3_ToF_plot(df_red, number_bins, range=[0, 71429])
        # Save histogram in ASCII-format
        dir_name = os.path.dirname(__file__)
        output_path = os.path.join(dir_name, '../output/ToF_%s.txt' % self.He3_data_sets[5:-5])
        np.savetxt(output_path,
                   np.transpose(np.array([bins, hist])),
                   delimiter=",",
                   header='ToF (µs), Counts')
        plt.legend()
        fig.show()

    def He3_Ch_action(self):
        fig = plt.figure()
        parameters = get_He3_filter_parameters(self)
        df_red = filter_He3(self.He3_df, parameters)
        He3_Ch_plot(df_red)
        fig.show()

    def He3_Energy_action(self):
        parameters = get_He3_filter_parameters(self)
        df_red = filter_He3(self.He3_df, parameters)
        number_bins = int(self.dE_bins.text())
        plot_energy = True
        fig = plt.figure()
        plt.subplot(1, 2, 1)
        hist_full, bins_full = energy_plot_He3(df_red, number_bins, plot_energy, 'Full Data')
        hist_pileup, bins_pileup = energy_plot_He3(df_red[df_red.PileUp == 1], number_bins, plot_energy, 'Pile Up Events')
        plt.legend()
        plt.subplot(1, 2, 2)
        plt.grid(True, which='major', linestyle='--', zorder=0)
        plt.grid(True, which='minor', linestyle='--', zorder=0)
        plt.xlabel('Energy (meV)')
        plt.ylabel('Fraction of Counts (PileUp/Full)')
        plt.title('Investigation of Pileup')
        plt.plot(bins_full, hist_pileup/hist_full, color='black', zorder=5,
                 label='PileUpEvents/AllEvents')
        plt.xscale('log')
        plt.legend()
        fig.show()

    def He3_Wavelength_action(self):
        parameters = get_He3_filter_parameters(self)
        df_red = filter_He3(self.He3_df, parameters)
        number_bins = int(self.dE_bins.text())
        plot_energy = False
        fig = plt.figure()
        plt.subplot(1, 2, 1)
        hist_full, bins_full = energy_plot_He3(df_red, number_bins, plot_energy, 'Full Data')
        hist_pileup, bins_pileup = energy_plot_He3(df_red[df_red.PileUp == 1], number_bins, plot_energy, 'Pile Up Events')
        plt.legend()
        plt.subplot(1, 2, 2)
        plt.grid(True, which='major', linestyle='--', zorder=0)
        plt.grid(True, which='minor', linestyle='--', zorder=0)
        plt.xlabel('Energy (meV)')
        plt.ylabel('Fraction of Counts (PileUp/Full)')
        plt.title('Investigation of Pileup')
        plt.plot(bins_full, hist_pileup/hist_full, color='black', zorder=5)
        fig.show()

    def ToF_MG_vs_ToF_He3_action(self):
        # Get filter parameters for MG and He-3
        parameters_MG = get_filter_parameters(self)
        parameters_He3 = get_He3_filter_parameters(self)
        # Filter data
        MG_red = filter_clusters(self.ce, parameters_MG)
        He3_red = filter_He3(self.He3_df, parameters_He3)
        # Declare paremeters
        number_bins = int(self.tofBins.text())
        MG_label, He3_label = self.data_sets, self.He3_data_sets
        useMaxNorm = True
        # Plot data
        fig = plt.figure()
        He3_ToF_plot(He3_red, number_bins, He3_label)
        ToF_histogram(MG_red, number_bins, MG_label)
        plt.legend()
        fig.show()

    def He3_pileup_action(self):
        # Filter data
        parameters_He3 = get_He3_filter_parameters(self)
        He3_red = filter_He3(self.He3_df, parameters_He3)
        # Plot data
        fig = plt.figure()
        fig.set_figheight(15)
        fig.set_figwidth(15)
        he3_pileup_plot(He3_red)
        plt.tight_layout()
        fig.show()

    # =========================================================================
    # Beam Monitors
    # =========================================================================

    def import_beam_monitor_data(self):
        dir_name = os.path.dirname(__file__)
        folder = os.path.join(dir_name, '../tables/')
        files = os.listdir(folder)
        files = [file for file in files if file[-4:] == '.asc']
        fig = plt.figure()
        for file in files:
            path = folder + file
            data = np.transpose(np.loadtxt(path, delimiter="\t"))
            x = data[0] * 1000
            y = data[1]
            norm = 1/max(y)
            plt.errorbar(x, y*norm, np.sqrt(y)*norm, label=file, fmt='.', capsize=5, zorder=5)
            plt.title('ToF - Beam Monitor')
            plt.xlabel('ToF [µs]')
            plt.ylabel('Counts (Normalized by maximum)')
            plt.yscale('log')
            plt.grid(True, which='major', linestyle='--', zorder=0)
            plt.grid(True, which='minor', linestyle='--', zorder=0)
            print('---')
            print(file)
            print(sum(y))
            print('---')
        fig.show()
        plt.legend()

    def initialize_bm_dict(self):
        """
        In order of appearance:

        1. Coated Radial, beam
        2. NonCoated Radial, beam
        3. He-3, beam
        4. Coated Radial, background
        5. NonCoated Radial, background
        6. He-3, background

        """
        bm_dict = {'mvmelst_165_191002_111641_Det2_overnight3.zip': 11411036,
                   'mvmelst_135_190930_141618_Det1_overnight2_30x80_14x60.zip': 9020907,
                   '2019_09_HZB_He3InBeam54304s_overnight.zip': 10723199,
                   'mvmelst_169_191003_075039_Det2_He3InBeam_overnight4.zip': 14052542,
                   'mvmelst_141_191001_120405_He3InBeam_overnight3.zip': 10723199,
                   '2019_09_HZB_He_3_background_beam_blocked_by_boron_cadmium_9208s.zip': 1809302
                   }
        self.BM_counts_dict = bm_dict



    # ========================================================================
    # Helper Functions
    # ========================================================================

    def setup_buttons(self):
        # File handling
        self.cluster_button.clicked.connect(self.cluster_action)
        self.save_button.clicked.connect(self.save_action)
        self.load_button.clicked.connect(self.load_action)
        # PHS
        self.PHS_1D_button.clicked.connect(self.PHS_1D_action)
        self.PHS_2D_button.clicked.connect(self.PHS_2D_action)
        self.PHS_wires_vs_grids_button.clicked.connect(self.PHS_wires_vs_grids_action)
        self.PHS_comparison_button.clicked.connect(self.PHS_comparison_action)
        # Misc
        self.multiplicity_button.clicked.connect(self.Multiplicity_action)
        self.ToF_button.clicked.connect(self.ToF_action)
        self.timestamp_button.clicked.connect(self.Timestamp_action)
        # Coincidences
        self.Coincidences_2D_button.clicked.connect(self.Coincidences_2D_action)
        self.Coincidences_3D_button.clicked.connect(self.Coincidences_3D_action)
        self.Coincidences_Projections_button.clicked.connect(self.Coincidences_Projections_action)
        self.CE_2D_comparison_button.clicked.connect(self.CE_2D_comparison_action)
        # Analysis
        self.wavelength_button.clicked.connect(self.Wavelength_action)
        self.energy_button.clicked.connect(self.Energy_action)
        self.count_rate_button.clicked.connect(self.Count_Rate_action)
        self.efficiency_button.clicked.connect(self.Efficiency_action)
        self.ToF_Overlay_button.clicked.connect(self.ToF_Overlay_action)
        self.Wavelength_Overlay_button.clicked.connect(self.Wavelength_overlay_action)
        self.layers_button.clicked.connect(self.Layers_action)
        self.full_analysis_button.clicked.connect(self.full_analysis_action)
        # Animation
        self.time_sweep_button.clicked.connect(self.time_sweep_action)
        self.lambda_sweep_button.clicked.connect(self.lambda_sweep_action)
        self.ToF_sweep_button.clicked.connect(self.ToF_sweep_action)
        # He-3 tubes
        self.he3_import_button.clicked.connect(self.Import_He3_action)
        self.he3_save_button.clicked.connect(self.Save_He3_action)
        self.he3_load_button.clicked.connect(self.Load_He3_action)
        self.he3_PHS_button.clicked.connect(self.He3_PHS_action)
        self.he3_ToF_button.clicked.connect(self.He3_ToF_action)
        self.he3_ch_button.clicked.connect(self.He3_Ch_action)
        self.he3_energy_button.clicked.connect(self.He3_Energy_action)
        self.he3_wavelength_button.clicked.connect(self.He3_Wavelength_action)
        self.ToF_MG_vs_ToF_He3_button.clicked.connect(self.ToF_MG_vs_ToF_He3_action)
        self.PileUp_button.clicked.connect(self.He3_pileup_action)

    def refresh_window(self):
        self.app.processEvents()
        self.update()
        self.app.processEvents()
        self.update()
        self.app.processEvents()
        self.app.processEvents()
        self.app.processEvents()

    def fill_MG_information_window(self):
        information_text = '<b>Measurement time:</b> %d [s]' % int(self.measurement_time)
        information_text += '<br/><b>Data sets:</b> ' + self.data_sets
        self.information_window.setText(information_text)

    def fill_He3_information_window(self):
        information_text = "<b>Counts:</b> %d [Counts]" % self.He3_counts
        information_text += '<br/><b>Data sets:</b> ' + self.He3_data_sets
        self.He3_information_window.setText(information_text)


# =============================================================================
# Start GUI
# =============================================================================

app = QApplication(sys.argv)
main_window = MainWindow(app)
main_window.setAttribute(Qt.WA_DeleteOnClose, True)
main_window.setup_buttons()
sys.exit(app.exec_())
