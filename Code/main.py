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
# File handling
from FileHandling.Import import unzip_data, import_data
from FileHandling.Cluster import cluster_data
from FileHandling.Storage import save_data, load_data
# He-3 tubes
from HeliumTubes.FileHandlingHe3 import unzip_He3_data, import_He3_data, save_He3_data, load_He3_data
from HeliumTubes.PlottingHe3 import He3_PHS_plot, He3_ToF_plot, He3_Ch_plot, energy_plot_He3
from HeliumTubes.FilteringHe3 import get_He3_filter_parameters, filter_He3
# Helper functions
from HelperFunctions.CreateMapping import create_mapping
from HelperFunctions.Filtering import filter_clusters, get_filter_parameters
from HelperFunctions.EnergyTransfer import calculate_energy
# PHS
from Plotting.PHS.PHS_1D import PHS_1D_plot
from Plotting.PHS.PHS_2D import PHS_2D_plot
from Plotting.PHS.PHS_Wires_Vs_Grids import PHS_wires_vs_grids_plot
# Coincidences
from Plotting.Coincidences.Coincidences_2D import coincidences_2D_plot
from Plotting.Coincidences.Coincidences_3D import coincidences_3D_plot
from Plotting.Coincidences.Coincidences_Projections import coincidences_projections_plot
# Misc
from Plotting.Misc.Multiplicity import multiplicity_plot
from Plotting.Misc.ToF import ToF_histogram
from Plotting.Misc.Timestamp import timestamp_plot
# Analysis
from Plotting.Analysis.DeltaE import energy_plot
from Plotting.Analysis.CountRate import calculate_count_rate
from Plotting.Analysis.Efficiency import calculate_efficiency
from Plotting.Analysis.EnergyResolution import calculate_energy_resolution
from Plotting.Analysis.Lineshape import analyze_Lineshape
# Animation
from Plotting.Analysis.Animation3D import Animation_3D_plot
from Plotting.Analysis.LambdaSweep import Lambda_Sweep_Animation
from Plotting.Analysis.TimeSweep import Time_Sweep_Animation


# TEMP
import matplotlib.pyplot as plt

# =============================================================================
#                                   Windows
# =============================================================================

class MainWindow(QMainWindow):
    def __init__(self, app, parent=None):
        super(MainWindow, self).__init__(parent)
        dir_name = os.path.dirname(__file__)
        title_screen_path = os.path.join(dir_name, '../Windows/mainwindow.ui')
        self.ui = uic.loadUi(title_screen_path, self)
        self.app = app
        # Clustering attributes
        self.data_sets = ''
        self.ILL_buses = [-1, -1, -1]
        self.maximum_file_size_in_mb = 3000
        self.adc_threshold = 0
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

    # =========================================================================
    # File handling
    # =========================================================================

    def cluster_action(self):
        zip_paths = QFileDialog.getOpenFileNames(self, "", "../Data")[0]
        if len(zip_paths) > 0:
            self.set_clustering_parameters()
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
            clusters, events = cluster_data(data, self.ILL_buses, self.adc_threshold)
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
            save_data(path, self.ce, self.e, self.data_sets, self.adc_threshold, self.ILL_buses)

    def load_action(self):
        path = QFileDialog.getOpenFileName(self, "", "../Data")[0]
        if path != '':
            clusters, events, data_sets_temp, adc_threshold_temp, ILL_buses_temp = load_data(path)
            # Write or append
            if self.write.isChecked():
                self.ce = clusters
                self.e = events
                self.data_sets = data_sets_temp
                self.adc_threshold = adc_threshold_temp
                self.ILL_buses = ILL_buses_temp
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


    # =========================================================================
    # Plotting
    # =========================================================================

    # ==== PHS ==== #

    def PHS_1D_action(self):
        if self.data_sets != '':
            filter_parameters = get_filter_parameters(self)
            ce_filtered = filter_clusters(self.ce, filter_parameters)
            number_bins = int(self.phsBins.text())
            fig = PHS_1D_plot(self.e, ce_filtered, number_bins)
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
            output_path = os.path.join(dir_name, '../Output/')
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
            if self.ESS_button.isChecked():
                detector_type = 'ESS'
            else:
                detector_type = 'ILL'
            coincidences_3D_plot(ce_filtered, detector_type, origin_voxel)

    def Coincidences_Projections_action(self):
        if self.data_sets != '':
            filter_parameters = get_filter_parameters(self)
            ce_filtered = filter_clusters(self.ce, filter_parameters)
            bus_start = self.module_min.value()
            bus_stop = self.module_max.value()
            fig, histograms = coincidences_projections_plot(ce_filtered, bus_start, bus_stop)
            # Export histograms to text
            dir_name = os.path.dirname(__file__)
            output_path = os.path.join(dir_name, '../Output/')
            file_names = ['Front', 'Top', 'Side']
            for file_name, histogram in zip(file_names, histograms):
                path = output_path + '2D_Coincidences_Projections_%s.txt' % file_name
                np.savetxt(path, histogram, fmt="%d", delimiter=",")
            # Plot figure
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
            filter_parameters = get_filter_parameters(self)
            ce_filtered = filter_clusters(self.ce, filter_parameters)
            number_bins = int(self.tofBins.text())
            fig = plt.figure()
            ToF_histogram(ce_filtered, number_bins)
            fig.show()

    def ToF_Overlay_action(self):
        paths = QFileDialog.getOpenFileNames(self, "", "../Data")[0]
        if len(paths) > 0:
            # Declare parameters
            time_offset = (0.6e-3) * 1e6
            period_time = (1/14) * 1e6
            filter_parameters = get_filter_parameters(self)
            number_bins = int(self.tofBins.text())
            # Plot
            fig = plt.figure()
            for i, path in enumerate(paths):
                label = path.rsplit('/', 1)[-1]
                ce = pd.read_hdf(path, 'ce')
                ce_filtered = filter_clusters(ce, filter_parameters)
                duration = get_duration(ce)
                ToF_shifted = (ce_filtered.ToF * 62.5e-9 * 1e6 + time_offset) % period_time
                hist, *_ = np.histogram(ToF_shifted, bins=2000)
                norm = 1/duration
                plt.hist(ToF_shifted, bins=number_bins, zorder=4, histtype='step',
                         label=label,
                         weights=(norm)*np.ones(len(ToF_shifted)))
            plt.title('ToF')
            plt.xlabel('ToF [$\mu$s]')
            plt.ylabel('Counts (Normalized by duration)')
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
            Time = (ce_filtered.Time*62.5e-9)/(60 ** 2)
            timestamp_plot(Time, number_bins, unit)
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
            if self.ESS_button.isChecked():
                detector_type = 'ESS'
            else:
                detector_type = 'ILL'
            fig = plt.figure()
            energy_plot(ce_filtered, detector_type, origin_voxel, number_bins)
            fig.show()

    def Count_Rate_action(self):
        if (self.data_sets != ''):
            filter_parameters = get_filter_parameters(self)
            ce_filtered = filter_clusters(self.ce, filter_parameters)
            ToF_values = ce_filtered.ToF * 62.5e-9
            count_rate = calculate_count_rate(ToF_values, self.measurement_time)
            print('Count rate: %.1f [Hz]' % count_rate)

    def Efficiency_action(self):
        if (self.data_sets != ''):
            parameters_MG = get_filter_parameters(self)
            parameters_He3 = get_He3_filter_parameters(self)
            MG_red = filter_clusters(self.ce, parameters_MG)
            He3_red = filter_He3(self.He3_df, parameters_He3)
            number_bins = int(self.dE_bins.text())
            origin_voxel = [int(self.bus_origin.text()),
                            int(self.gCh_origin.text()),
                            int(self.wCh_origin.text())]
            if self.ESS_button.isChecked():
                detector_type = 'ESS'
            else:
                detector_type = 'ILL'
            MG_label, He3_label = self.data_sets, self.He3_data_sets
            start, stop = 1, 10
            useMaxNorm = True
            fig = plt.figure()
            fig.suptitle('Comparison He-3 and MG (normalized to max value)')
            plt.subplot(1, 2, 1)
            plot_energy = True
            energy_plot(MG_red, detector_type, origin_voxel, number_bins,
                        start, stop, plot_energy, MG_label, useMaxNorm)
            energy_plot_He3(He3_red, number_bins, plot_energy, He3_label, useMaxNorm)
            plt.legend()
            plt.subplot(1, 2, 2)
            plot_energy = False
            energy_plot(MG_red, detector_type, origin_voxel, number_bins,
                        start, stop, plot_energy, MG_label, useMaxNorm)
            energy_plot_He3(He3_red, number_bins, plot_energy, He3_label, useMaxNorm)
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
        paths = QFileDialog.getOpenFileNames(self, "", "../Data")[0]
        if len(paths) > 0:
            path = paths[0]
            # Declare parameters
            filter_parameters = get_filter_parameters(self)
            number_bins = int(self.dE_bins.text())
            origin_voxel = [int(self.bus_origin.text()),
                            int(self.gCh_origin.text()),
                            int(self.wCh_origin.text())]
            if self.ESS_button.isChecked():
                detector_type = 'ESS'
            else:
                detector_type = 'ILL'
            # Plot
            ce_MG = pd.read_hdf(path, 'ce')
            ce_MG_filtered = filter_clusters(ce_MG, filter_parameters)
            parameters = get_He3_filter_parameters(self)
            ce_He3_filtered = filter_He3(self.He3_df, parameters)
            fig = plt.figure()
            analyze_Lineshape(ce_MG_filtered, ce_He3_filtered, detector_type, origin_voxel)
            fig.show()

    # ==== Animation ==== #

    def Animation_3D_action(self):
        filter_parameters = get_filter_parameters(self)
        ce_filtered = filter_clusters(self.ce, filter_parameters)
        origin_voxel = [int(self.bus_origin.text()),
                        int(self.gCh_origin.text()),
                        int(self.wCh_origin.text())]
        if self.ESS_button.isChecked():
            detector_type = 'ESS'
        else:
            detector_type = 'ILL'
        Animation_3D_plot(ce_filtered, detector_type, origin_voxel)

    def lambda_sweep_action(self):
        filter_parameters = get_filter_parameters(self)
        ce_filtered = filter_clusters(self.ce, filter_parameters)
        origin_voxel = [int(self.bus_origin.text()),
                        int(self.gCh_origin.text()),
                        int(self.wCh_origin.text())]
        number_bins = int(self.dE_bins.text())
        bus_start = self.module_min.value()
        bus_stop = self.module_max.value()
        if self.ESS_button.isChecked():
            detector_type = 'ESS'
        else:
            detector_type = 'ILL'
        Lambda_Sweep_Animation(ce_filtered, number_bins, detector_type,
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
        if self.ESS_button.isChecked():
            detector_type = 'ESS'
        else:
            detector_type = 'ILL'
        Time_Sweep_Animation(ce_filtered, number_bins, detector_type,
                             origin_voxel, bus_start, bus_stop)





    # =========================================================================
    # He-3 tubes
    # =========================================================================

    def Import_He3_action(self):
        file_path = QFileDialog.getOpenFileName(self, "", "../Data")[0]
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
        path = QFileDialog.getOpenFileName(self, "", "../Data")[0]
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
        He3_ToF_plot(df_red, number_bins)
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
        energy_plot_He3(df_red, number_bins, plot_energy)
        fig.show()

    def He3_Wavelength_action(self):
        parameters = get_He3_filter_parameters(self)
        df_red = filter_He3(self.He3_df, parameters)
        number_bins = int(self.dE_bins.text())
        plot_energy = False
        fig = plt.figure()
        energy_plot_He3(df_red, number_bins, plot_energy)
        fig.show()

    def ToF_MG_vs_ToF_He3_action(self):
        # Get filter parameters for MG and He-3
        parameters_MG = get_filter_parameters(self)
        parameters_He3 = get_He3_filter_parameters(self)
        # Filter data
        MG_red = filter_clusters(self.ce, parameters_MG)
        He3_red = filter_He3(self.He3_df, parameters_He3)
        # Declare paremeters
        number_bins = int(self.dE_bins.text())
        MG_label, He3_label = self.data_sets, self.He3_data_sets
        # plot data
        fig = plt.figure()
        He3_ToF_plot(He3_red, number_bins, He3_label)
        ToF_histogram(MG_red, number_bins, MG_label)
        plt.legend()
        fig.show()



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
        # Misc
        self.multiplicity_button.clicked.connect(self.Multiplicity_action)
        self.ToF_button.clicked.connect(self.ToF_action)
        self.timestamp_button.clicked.connect(self.Timestamp_action)
        # Coincidences
        self.Coincidences_2D_button.clicked.connect(self.Coincidences_2D_action)
        self.Coincidences_3D_button.clicked.connect(self.Coincidences_3D_action)
        self.Coincidences_Projections_button.clicked.connect(self.Coincidences_Projections_action)
        # Analysis
        self.energy_button.clicked.connect(self.Energy_action)
        self.count_rate_button.clicked.connect(self.Count_Rate_action)
        self.efficiency_button.clicked.connect(self.Efficiency_action)
        self.ToF_Overlay_button.clicked.connect(self.ToF_Overlay_action)
        self.lineshape_button.clicked.connect(self.Lineshape_action)
        # Animation
        self.time_sweep_button.clicked.connect(self.time_sweep_action)
        self.lambda_sweep_button.clicked.connect(self.lambda_sweep_action)
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
        # Button toogle
        self.toogle_VMM_MG()

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
        information_text += "<br/><b>ADC Threshold:</b> %d [ADC Ch's]" % self.adc_threshold
        information_text += '<br/><b>ILL buses:</b> ' + str(self.ILL_buses)
        information_text += '<br/><b>Data sets:</b> ' + self.data_sets
        self.information_window.setText(information_text)

    def fill_He3_information_window(self):
        information_text = "<b>Counts:</b> %d [Counts]" % self.He3_counts
        information_text += '<br/><b>Data sets:</b> ' + self.He3_data_sets
        self.He3_information_window.setText(information_text)

    def set_clustering_parameters(self):
        self.ILL_buses = [self.ILL_bus_1.value(), self.ILL_bus_2.value(), self.ILL_bus_3.value()]
        self.maximum_file_size_in_mb = float(self.maximum_file_size_in_mb_value.text())
        self.adc_threshold = float(self.adc_threshold_value.text())

    def toogle_VMM_MG(self):
        self.ESS_button.toggled.connect(
            lambda checked: checked and self.ILL_button.setChecked(False))
        self.ILL_button.toggled.connect(
            lambda checked: checked and self.ESS_button.setChecked(False))



# =============================================================================
# Helper Functions
# =============================================================================

def append_folder_and_files(folder, files):
    folder_vec = np.array(len(files)*[folder])
    return np.core.defchararray.add(folder_vec, files)

def get_duration(df):
    times = df.Time.values
    diff = np.diff(times)
    resets = np.where(diff < 0)
    duration_in_TDC_channels = sum(times[resets]) + times[-1]
    duration_in_seconds = duration_in_TDC_channels * 62.5e-9
    return duration_in_seconds


# =============================================================================
# Start GUI
# =============================================================================

app = QApplication(sys.argv)
main_window = MainWindow(app)
main_window.setAttribute(Qt.WA_DeleteOnClose, True)
main_window.setup_buttons()
sys.exit(app.exec_())
