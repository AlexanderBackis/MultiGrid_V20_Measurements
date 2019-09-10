# Multi-Grid V20 Measurements, 23/9-2019 -> 6/10-2019

Application for analysis of Multi-Grid data taken at the HZB using the V20 instrument.
The program consists of a GUI Interface which allows the user to cluster Mesytec VMMR-8/16 output and analyse the data using different tools, such as:

- PHS (Cumulative or Individual Channel)
- Coincidences (2D and 3D)
- Time-of-Flight

The program has a 'filtering'-feature which allows the user to gate the data on different parameters, such as PHS, ToF and position. In addition to this, it conatains the relevant analysis features used at the measurement:

1. Efficiency at cold neutrons energies
2. Detector response at high instanteous flux
3. Elastic line shape in energy transfer spectra

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
Navigate to MultiGrid_V20_Measurements->Code and enter:
```
python main.py
```
## Notes

The code requires the following files:
- mapping.xlsx
- window.ui

These can be found under the 'Tables'- and 'Window'-folders in the repository, and the files can be manipulated according to the specific conditions of the measurement.
