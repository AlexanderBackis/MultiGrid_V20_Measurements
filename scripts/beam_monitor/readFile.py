################################################################################################################
# Converter for TSDAU bn4 files
# 
# creates single ASCII files for every channel in bn4 file
# not used chnnels are skipped
#
# dependecy on matplotlib ony for writing channel spectrum as image file
# 
# author: Christian Jacobsen (christian.jacobsen@hzg.de)
# date: 26-07-2017
# version: 0.5
#
################################################################################################################


import matplotlib.pyplot as plt
from struct import *
import codecs
import sys
import numpy as np

def read(filename, binSize, channels, maxLength):

	#print(sys.path)

	#filename = '/home/jacobsec/workspace/Mantid/tsdau/Spectrum12.bn4'
	#binSize = 2500
	#channels = 8
	#maxLength = 720000	# in 100ns unit

	print('read File')



	bins = [[0]*(binSize) for i in range(channels+1)]

	hasErrors = False
	cntGood = 0
	cntBad = 0


	with open(filename, "rb") as f:

		# search for end of headerblock
		while(True):


			b = f.read(4)

			if b == b'\x00\x00\x00\x00':	# typical pattern for first data
				break

		longtime_ts = 0	# holds high bytes of timestamp 1,6s unit
		trigger_ts = 0	# holds timestamp of last trigger event

		f.read(4)

		# read data block
		while(True):
			#b = f.read(4)	# 4 bytes contains single event, 1 byte channelnumber, 3 byte ts
			read = f.read(1)

			#if chNum == x'ff':
			#	print('test')

			if len(read) == 0:
				break

			chNum = int(codecs.encode(read, 'hex'), 16)

			read = f.read(3)

			if len(read) < 3:	# found end of file
				break

			if chNum >= 0 and chNum <= channels:	# select only data channels
				chTs = int(codecs.encode(read, 'hex'), 16)
				ts = (longtime_ts<<24) + chTs	# extracts ts from data chunk
				#print(ts)
				idx = int((float(ts-trigger_ts)/float(maxLength)*float(binSize)))+2	# calculate inex of bin in array
				#print(idx)
				if chNum == 4:		#Trigger

					trigger_ts = ts

				# skip saving ts if no trigger occure before, or channel is not valid data 
				if trigger_ts != 0 and chNum != 4 and chNum != 0:
					if idx >= 0 and idx < binSize:	
						bins[chNum][idx] += 1
						cntGood += 1
					else:
						hasErrors = True
						cntBad += 1
					#	print(chNum)
					#	print(idx)						

				#break

			else:
				#print(type(b[2]))
				chTs = int(codecs.encode(read, 'hex'), 16)
				#chTs = chTs & 0xffff
				longtime_ts = chTs#(b[2]<<8) + b[3]	# longtime part of ts
				#print(longtime_ts)
				
	c = 0
	summe = 0
	# save data to file
	print('convert ...')
	for n in bins:

		#print(bins)

		print('... channel ' + str(c+1))

		
		c=c+1
		summe = 0

		if all(v == 0 for v in n):		# skip this channel if array contains no data
			continue

		# plot and save to file as image
		plt.clf()
		fig = plt.figure()
		plt.plot(np.arange(0, 400, 1)*17000/500, n, zorder=5, color='black')
		plt.ylabel('Counts')
		plt.xlabel('ToF [µs]')
		plt.title(filename.rsplit('/', 1)[-1])
		plt.grid(True, which='major', linestyle='--', zorder=0)
		plt.grid(True, which='minor', linestyle='--', zorder=0)
		fig.show()

		#plt.show()
		plt.savefig(filename + 'ch' + str(c-1) + '_bin' + str(binSize)+'.pdf', bbox_inches='tight')

		# save to ascii file
		with open(filename + 'ch' + str(c-1) + '_bin' + str(binSize)+'.asc', "w") as f:

			for i in range(len(n)):
				#print(i)
				#print(binSize)
				#print(maxLength)
				file_ts = float(i)/binSize*(maxLength/10.0)/1000.0 + ((maxLength/10.0/binSize)/2.0)/1000.0
				#print(n[i])
				#print(file_ts)

				f.write('{:06.3f}\t{:4d}\r\n'.format(file_ts, n[i]))

				summe += n[i]

			#print(str(c-1) + " " + str(summe))
		
	if hasErrors:
		print('File may be (partial) corrupt, ' + str(cntBad) + ' events of ' + str(cntGood) + ' ignored')
	print('End')
