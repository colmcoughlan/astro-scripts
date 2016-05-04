#! /usr/bin/env python
# PyBDSM starts counting islands at zero, while Duchamp starts at 1. Flagged regions =-1 in PyBDSM, 0 in Duchamp
# (Using PyBDSM process_image task)
# This script just adds one to the entire PyBDSM map, converting it to a Duchamp-style map
# CASA masks are binary and can also be converted to Duchamp style maps with option=1
# Colm Coughlan
# 12.2.15
import sys
import numpy as np
from scipy import ndimage
import pyfits
# Check arguments
if(len(sys.argv)!=3):
	print("\tError: Takes two arguments.")
	print("\tUseage: pybdsm_to_duchamp_mask <filename> <type>")
	print("\tWARNING: Developed for older version of PyBDSM masks (where they were more like duchamp). Optype 2 should still be useful.")
	print("\tUseage: Type = 0 for old PyBDSM to Duchamp mask. Type = 1 for CASA to Duchamp mask. Type = 2 for old PyBDSM to CASA conversion.")
	sys.exit()
else:
	inputname = str(sys.argv[1])
	type = int(sys.argv[2])
	print('\tConverting '+inputname+'.')
# Open fits file and add one to data area
# Update existing file and close it

if type==0:
	print('\tAssuming currently in PyBDSM format.')
	print('\tConverting to Duchamp format...')
	f = pyfits.open(inputname,mode='update')
	f[0].data = np.add(f[0].data,1)
	f.flush()
	f.close()
	print('\tProcess complete.')
elif type==1:
	print('\tAssuming currently in CASA format.')
	print('\tConverting to Duchamp format...')
	f = pyfits.open(inputname,mode='update')
	original_data = f[0].data[0,0,:,:]
	labelled_array, n_islands = ndimage.measurements.label(original_data)
	print('\tNumber of islands found = '+str(n_islands))
	f[0].data[0,0,:,:] = labelled_array
	f.flush()
	f.close()
	print('\tProcess complete.')
else:
	print('\tAssuming currently in old PyBDSM format.')
	print('\tConverting to CASA format...')
	f = pyfits.open(inputname,mode='update')
	f[0].data = np.add(f[0].data,1)
	f[0].data = np.clip(f[0].data, 0 , 1)
	f.flush()
	f.close()
	print('\tProcess complete.')
