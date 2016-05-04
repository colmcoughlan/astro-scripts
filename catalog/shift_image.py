#! /usr/bin/env python
# Colm Coughlan
# Shift a FITS file along a given axis
# 12.10.2015
import shutil
import sys
import numpy as np
import pyfits

# Check arguments
if(len(sys.argv)!=5):
	print('\tError: Takes three arguments.')
	print('\tUseage: shift_image <filename> <no. of pixels to shift> <axis> <outputname>')
	print('\tIn a normal FITS image, axis = 0 => DEC, axis = 1 => RA')
	print('\tA positive number of pixels shifts sources upwards.')
	sys.exit()
else:
	inputname = str(sys.argv[1])
	shift_pix = int(sys.argv[2])
	axis_val = int(sys.argv[3])
	if (axis_val!=0 and axis_val!=1):
		print('\tError:Axis val must be either 0 or 1')
		sys.exit()
	outputname = str(sys.argv[4])

print('\tShifting '+inputname+' by '+str(shift_pix)+' pixel along axis '+str(axis_val)+' and saving result in '+outputname)
print('\tReplacement pixels given a value of NaN (should be interpreted as blanked)')

# make new file

shutil.copyfile(inputname, outputname)

# access data image

f = pyfits.open(outputname,mode='update')

# Perform shift

tdata = f[0].data[0,0,:,:]
tdata = np.roll(tdata, shift_pix, axis = axis_val)

# Zero the appropriate data

if(axis_val == 0):
	if(shift_pix >= 0):
		tdata[0:shift_pix,:] = np.nan
	else:
		n = tdata.shape[0]
		tdata[n + shift_pix:n,:] = np.nan
else:
	if(shift_pix >= 0):
		tdata[:,0:shift_pix] = np.nan
	else:
		n = tdata.shape[1]
		tdata[:,n + shift_pix:n] = np.nan

f[0].data[0,0,:,:] = tdata

# Update history

f[0].header.add_history('shift_image.py: This file has been shifted '+str(shift_pix)+' pixels along axis '+str(axis_val))
f[0].header.add_history('shift_image.py: No changes have been made to the header information')

# Flush changes and close

f.flush()
f.close()

