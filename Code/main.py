from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5 import uic
import os
import sys
import pandas as pd

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
        self.measurement_time = 0
        self.data_sets = ''
        self.show()
        self.refresh_window()

    # =========================================================================
    # File handling
    # =========================================================================

    def cluster_action(self):
        pass


    # =========================================================================
    # Plotting
    # =========================================================================

    def PHS_1D_action(self):
        if self.data_sets != '':
            pass

    def PHS_2D_action(self):
        if self.data_sets != '':
            pass

    def ToF_action(self):
        if self.data_sets != '':
            pass

    def Coincidences_2D_action(self):
        if self.data_sets != '':
            pass

    def Coincidences_3D_action(self):
        if self.data_sets != '':
            pass

    def Coincidences_Front_Top_Side_action(self):
        if self.data_sets != '':
            pass

    # ========================================================================
    # Helper Functions
    # ========================================================================

    def setup_buttons(self):
        # File handling
        #self.cluster_button.clicked.connect(self.cluster_action)
        # Plotting
        #self.PHS_1D_button.clicked.connect(self.PHS_1D_action)
        #self.PHS_2D_button.clicked.connect(self.PHS_2D_action)
        #self.Coincidences_2D_button.clicked.connect(self.Coincidences_2D_action)
        #self.Coincidences_3D_button.clicked.connect(self.Coincidences_3D_action)
        #self.Coincidences_Front_Top_Side_button.clicked.connect(self.Coincidences_Front_Top_Side_action)
        #self.ToF_button.clicked.connect(self.ToF_action)
        pass

    def refresh_window(self):
        self.app.processEvents()
        self.update()
        self.app.processEvents()
        self.update()
        self.app.processEvents()
        self.app.processEvents()
        self.app.processEvents()


# =============================================================================
# Helper Functions
# =============================================================================

def append_folder_and_files(folder, files):
    folder_vec = np.array(len(files)*[folder])
    return np.core.defchararray.add(folder_vec, files)

# =============================================================================
# Start GUI
# =============================================================================

app = QApplication(sys.argv)
main_window = MainWindow(app)
main_window.setAttribute(Qt.WA_DeleteOnClose, True)
main_window.setup_buttons()
sys.exit(app.exec_())
