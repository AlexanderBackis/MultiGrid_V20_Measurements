#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" main.py: Module containing the GUI used for data analysis.
"""

# Standard library
import os
import sys
# QT
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5 import uic
# Data Analysis
import pandas as pd
# Local sources
from FileHandling.Import import unzip_data, import_data
from FileHandling.Cluster import cluster_data
from HelperFunctions.CreateMapping import create_full_mapping
from HelperFunctions.Filtering import filter_clusters
# PHS
from Plotting.PHS.PHS_1D import PHS_1D_plot
from Plotting.PHS.PHS_2D import PHS_2D_plot
from Plotting.PHS.PHS_Wires_Vs_Grids import PHS_wires_vs_grids_plot
# Misc
from Plotting.Misc.Multiplicity import multiplicity_plot
from Plotting.Misc.ToF import ToF_histogram

# =============================================================================
# Windows
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
        self.Ei = -1
        self.ce = pd.DataFrame()
        self.e = pd.DataFrame()
        self.fill_information_window()
        self.show()
        self.refresh_window()

    # =========================================================================
    # File handling
    # =========================================================================

    def cluster_action(self):
        zip_paths = QFileDialog.getOpenFileNames(self, "Select Directory", "../Data")[0]
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
            self.fill_information_window()
            self.refresh_window()
            print(self.ce)
            print(self.e)


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
            fig = PHS_2D_plot(self.e)
            fig.show()

    def PHS_wires_vs_grids_action(self):
        if (self.data_sets != ''):
            filter_parameters = get_filter_parameters(self)
            ce_filtered = filter_clusters(self.ce, filter_parameters)
            fig = PHS_wires_vs_grids_plot(ce_filtered)
            fig.show()

    # ==== Misc ==== #

    def Multiplicity_action(self):
        if (self.data_sets != ''):
            filter_parameters = get_filter_parameters(self)
            ce_filtered = filter_clusters(self.ce, filter_parameters)
            fig = multiplicity_plot(ce_filtered)
            fig.show()

    def ToF_action(self):
        if self.data_sets != '':
            filter_parameters = get_filter_parameters(self)
            ce_filtered = filter_clusters(self.ce, filter_parameters)
            number_bins = int(self.tofBins.text())
            fig = ToF_histogram(ce_filtered, number_bins)
            fig.show()

    def Coincidences_2D_action(self):
        if self.data_sets != '':
            pass

    def Coincidences_3D_action(self):
        if self.data_sets != '':
            pass

    def Coincidences_Projections_action(self):
        if self.data_sets != '':
            pass

    # ========================================================================
    # Helper Functions
    # ========================================================================

    def setup_buttons(self):
        # File handling
        self.cluster_button.clicked.connect(self.cluster_action)
        # PHS
        self.PHS_1D_button.clicked.connect(self.PHS_1D_action)
        self.PHS_2D_button.clicked.connect(self.PHS_2D_action)
        self.PHS_wires_vs_grids_button.clicked.connect(self.PHS_wires_vs_grids_action)
        # Misc
        self.multiplicity_button.clicked.connect(self.Multiplicity_action)
        self.ToF_button.clicked.connect(self.ToF_action)
        # Coincidences
        #self.Coincidences_2D_button.clicked.connect(self.Coincidences_2D_action)
        #self.Coincidences_3D_button.clicked.connect(self.Coincidences_3D_action)
        #self.Coincidences_Projections.clicked.connect(self.Projections_action)

    def refresh_window(self):
        self.app.processEvents()
        self.update()
        self.app.processEvents()
        self.update()
        self.app.processEvents()
        self.app.processEvents()
        self.app.processEvents()

    def fill_information_window(self):
        information_text = '<b>Measurement time:</b> %d [s]' % int(self.measurement_time)
        information_text += '<br/><b>Incident energy:</b> %.2f [meV]' % self.Ei
        information_text += "<br/><b>ADC Threshold:</b> %d [ADC Ch's]" % self.adc_threshold
        information_text += '<br/><b>ILL buses:</b> ' + str(self.ILL_buses)
        information_text += '<br/><b>Data sets:</b> ' + self.data_sets
        self.information_window.setText(information_text)

    def set_clustering_parameters(self):
        self.ILL_buses = [self.ILL_bus_1.value(), self.ILL_bus_2.value(), self.ILL_bus_3.value()]
        self.Ei = float(self.Ei_value.text())
        self.maximum_file_size_in_mb = float(self.maximum_file_size_in_mb_value.text())
        self.adc_threshold = float(self.adc_threshold_value.text())



# =============================================================================
# Helper Functions
# =============================================================================

def append_folder_and_files(folder, files):
    folder_vec = np.array(len(files)*[folder])
    return np.core.defchararray.add(folder_vec, files)

def get_filter_parameters(window):
    parameters = {'wM': [window.wM_min.value(),
                         window.wM_max.value(),
                         window.wM_filter.isChecked()],
                  'gM': [window.gM_min.value(),
                         window.gM_max.value(),
                         window.gM_filter.isChecked()],
                  'ceM': [window.ceM_min.value(),
                          window.ceM_max.value(),
                          window.ceM_filter.isChecked()],
                  'wADC': [float(window.wADC_min.text()),
                           float(window.wADC_max.text()),
                           window.wADC_filter.isChecked()],
                  'gADC': [float(window.gADC_min.text()),
                           float(window.gADC_max.text()),
                           window.gADC_filter.isChecked()],
                  'ToF': [float(window.ToF_min.text()) / (62.5e-9 * 1e6),
                          float(window.ToF_max.text()) / (62.5e-9 * 1e6),
                          window.ToF_filter.isChecked()],
                  'Time': [float(window.Time_min.text()),
                           float(window.Time_max.text()),
                           window.Time_filter.isChecked()],
                  'Bus': [window.module_min.value(),
                          window.module_max.value(),
                          window.module_filter.isChecked()],
                  'wire': [window.wire_min.value(),
                           window.wire_max.value(),
                           window.wire_filter.isChecked()],
                  'gCh': [window.grid_min.value() + 80 - 1,
                          window.grid_max.value() + 80 - 1,
                          window.grid_filter.isChecked()]
                  }
    return parameters

# =============================================================================
# Start GUI
# =============================================================================

app = QApplication(sys.argv)
main_window = MainWindow(app)
main_window.setAttribute(Qt.WA_DeleteOnClose, True)
main_window.setup_buttons()
sys.exit(app.exec_())
