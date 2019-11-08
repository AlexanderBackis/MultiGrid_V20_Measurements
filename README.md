# Multi-Grid V20 Measurements, 23/9-2019 -> 6/10-2019

Application for analysis of Multi-Grid data taken at the HZB using the V20 instrument.
The program consists of a GUI Interface which allows the user to cluster Mesytec VMMR-8/16 output and analyse the data using different tools, such as:

- PHS
- Coincidences (2D and 3D)
- Time-of-Flight
- Energy and Wavelength

In addition to this, the program has a 'filtering'-feature which allows the user to gate the data on different parameters, such as the ones mentioned above. It also contains the relevant analysis code, written for the three main purposes of the measurement:

1. Efficiency at cold neutrons energies
2. Elastic line shape in energy spectra
3. Detector response at high instanteous flux

## Requisties
- Python3 (https://www.python.org/downloads/)
- Anaconda (https://www.anaconda.com/distribution/)

## Installation
Install dependencies:
```
conda install -c anaconda pyqt 
conda install -c plotly plotly
```

Clone the repository:
```
git clone https://github.com/AlexanderBackis/MultiGrid_V20_Measurements.git
```

## Execution
Navigate to MultiGrid_V20_Measurements->scripts and enter:
```
python main.py
```
