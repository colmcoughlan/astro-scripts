#! /usr/bin/env python
# Colm Coughlan
# 25.9.15
import shutil
import sys
import numpy as np
import pyfits
# Check arguments
if(len(sys.argv)!=4):
	print("\tError: Takes three arguments.")
	print("\tUseage: clip_image_size_to_match <filename> <file_to_be_matched> <outputname>")
	print("\tNote: Matched file only used to determine image size.")
	sys.exit()
else:
	inputname = str(sys.argv[1])
	comparision_name = str(sys.argv[2])
	outputname = str(sys.argv[3])

# make new file

shutil.copyfile(inputname, outputname)

# get new size

fin = pyfits.open(comparision_name,mode='readonly')
new_shape = np.shape(fin[0].data)
fin.close()

# resize image

fin = pyfits.open(inputname,mode='readonly')
fout = pyfits.open(outputname,mode='update')
old_shape = np.shape(fin[0].data)
print('Converting to '+str(new_shape))
w_clip = int(0.5*(old_shape[2] - new_shape[2]))
h_clip = int(0.5*(old_shape[3] - new_shape[3]))
#print(str(np.shape(fin[0].data[:,:,(w_clip):(old_shape[2]-w_clip),(h_clip):(old_shape[3]-h_clip)])))
fout[0].data = fin[0].data[:,:,(w_clip):(old_shape[2]-w_clip),(h_clip):(old_shape[3]-h_clip)]

# Update header information to account for clipping

fout[0].header.update('NAXIS1',new_shape[2])
fout[0].header.update('NAXIS2',new_shape[3])

fout[0].header.update('CRPIX1', int(fout[0].header.get('CRPIX1'))-w_clip )
fout[0].header.update('CRPIX2', int(fout[0].header.get('CRPIX2'))-h_clip )

fout.flush()

fin.close()
fout.close()
