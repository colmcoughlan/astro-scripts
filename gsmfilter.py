# Colm Coughlan 20.10.2014
# Read LOFAR GSM and clip at specified flux
# Colm Coughlan 3.11.2014
# Add ability to clip sources outside specified region (works for Northern hem. sources only atm)

import sys
from math import sqrt

# String for DGTAU based off Rachael's data in thesis

dgtau_str='DGTAU, POINT, 04:27:04.7, +26.06.16.3, 1.3976, , , , , [-0.84]'

# function to convert RA and DEC strings into degree floats

def ra_to_degrees(ra):
	hours = float(ra[0:2])
	minutes = float(ra[2:6])
	return(hours*15 + minutes*0.25)

def dec_to_degrees(dec):
	degrees = float(dec[0:2])
	minutes = float(dec[2:4])
	return(minutes/60.0 + degrees)

# function to find the separtion between two points, given their RA and DEC in strings

def separation_degrees(point1, point2):
	string = point1.split('+')
	ra1 = ra_to_degrees(string[0])
	dec1 = dec_to_degrees(string[1])

	string = point2.split('+')
	ra2 = ra_to_degrees(string[0])
	dec2 = dec_to_degrees(string[1])

	return( sqrt(pow( ra1 - ra2 ,2) + pow( dec1 - dec2 ,2)) )
	

# Check arguments

if(len(sys.argv)!=7):
	print "\tError: Takes six arguments."
	print "\tUseage: gsmfilter <filename> <flux cutoff> <outputname> <addsource> <phase centre> <radius>"
	print "\taddsource 0: none"
	print "\taddsource 1: DGTAU"
	print "\tPhase centre: eg. 0427+2606, 0 ==> no clip"
	print "\tAllowed radius (degrees)"
	sys.exit()
else:
	inputname = str(sys.argv[1])
	cutoff = float(sys.argv[2])
	outputname = str(sys.argv[3])
	addsource = int(sys.argv[4])
	phase_centre = str(sys.argv[5])
	radius = float(sys.argv[6])
	print '\tReading:', inputname
	print '\tWriting:', outputname
	print '\tClipping at:', str(cutoff)
	if(addsource == 1):
		print '\tAdding DGTAU --> '+dgtau_str+' Jy.'
	if(phase_centre == '0'):
		print '\tNot clipping by region.'
	else:
		print '\tPhase centre: '+phase_centre
		print '\tAllowed radius: '+str(radius)+' degrees.'


# Attempt to open input file

try:
	f = open(inputname)
except:
	print "\tError opening ", inputname
	sys.exit()


# Read data and close input file

lines = f.readlines()
f.close()

# open output file
# Find start of data in file
f = open(outputname,'w')

nlines = len(lines)
ctr = 0
while(lines[ctr][0]!='#'):
	f.write(lines[ctr])
	ctr = ctr + 1

if(phase_centre != '0'):
	f.write('# Clipped from orignal by gsmfilter.py, with Stokes I cutoff = '+str(cutoff)+' Jy and clipping all points outside '+str(radius)+' degrees of '+phase_centre+'.\n')
else:
	f.write('# Clipped from orignal by gsmfilter.py with Stokes I cutoff = '+str(cutoff)+' Jy.\n')

ctr = ctr + 1
print "\tInput points = ", str(nlines - ctr)

# Split by comma to find flux, if flux > cutoff flux, include in output
# close output file


out_ctr=0
for ctr in range(ctr,nlines):
	string = lines[ctr].split(',')
	if (float(string[4]) > cutoff):
		if (phase_centre == '0'):
			f.write(lines[ctr])
			out_ctr = out_ctr + 1
		else:
			if ( separation_degrees(phase_centre, string[0]) < radius):
				f.write(lines[ctr])
				out_ctr = out_ctr + 1

# Write out additional source if requested (only DGTAU supported right now)

if(addsource == 1):
	f.write('# gsmfilter.py: Adding DGTAU\n')
	f.write(dgtau_str)
	out_ctr = out_ctr + 1

print "\tOutput points = ", str(out_ctr)
f.close()
