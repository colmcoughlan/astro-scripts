#! /usr/bin/env python
# Colm Coughlan. Write out FITS file flux as csv
import sys
import numpy as np
import pyfits
import matplotlib.pyplot as plt
# Check arguments
if(len(sys.argv)!=5):
	print("\tError: Takes 4 arguments.")
	print("\tUseage: fits_to_tsv <inname1> <inname2> <inname3> <outname>.")
	print("\tWARNING: Assumes CASA FITS image.")
	sys.exit()
else:
	inputname1 = str(sys.argv[1])
	inputname2 = str(sys.argv[2])
	inputname3 = str(sys.argv[3])
	outputname = str(sys.argv[4])

plt.rcParams.update({'font.size': 15})

f = pyfits.open(inputname1,mode='readonly')
data1 = (f[0].data[0,0,:,:]).flatten()
f.close()

f = pyfits.open(inputname2,mode='readonly')
data2 = (f[0].data[0,0,:,:]).flatten()
f.close()

f = pyfits.open(inputname3,mode='readonly')
data3 = (f[0].data[0,0,:,:]).flatten()
f.close()

data = np.concatenate((data1,data2,data3), axis=0) # combine fields
data = data[~np.isnan(data)]	# remove NANs
data = np.multiply(data, 1.0e6)	# convert to microJy/beam
xmax = np.max(data)
xmin = 0.5*np.min(data)
step = 400

hist, bin_edges = np.histogram(data, bins = np.linspace(xmin, xmax, num=step))
hist = np.divide(hist, float(len(data)))
cdf = np.cumsum(hist)
plt.bar(bin_edges[:-1], hist, width = np.abs((xmax-xmin)/step))
ax1 = plt.gca()
ax1.set_ylabel(r'$\mathrm{Fraction}$')
ax1.set_xlabel(r'$\mathrm{RMS\,Noise\,(}\mu\mathrm{Jy})$')

ax2 = ax1.twinx()
ax2.plot(bin_edges[:-1], cdf, color='red')
ax2.set_ylabel(r'$\mathrm{Cumulative\,Fraction}$')
ax2.set_ylim(0.0,1.0)

plt.xlim(xmin, 25.0*xmin)

plt.savefig(outputname)

print('Median = '+str(np.median(data))+' microJy.')
print('\tProcess complete.')

