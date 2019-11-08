from readFile import read
import os

# Get paths
dir_name = os.path.dirname(__file__)
folder = os.path.join(dir_name, 'data/')
files = os.listdir(folder)
files = [file for file in files if file[-4:] == '.bn4']
# Declare parameters
binSize = 400
channels = 8
maxLength = 720000	 # in 100ns unith
# Iterate through data
for file in files:
	path = folder + file
	read(path, binSize, channels, maxLength)
