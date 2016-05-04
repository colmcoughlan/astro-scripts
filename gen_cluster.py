# Colm Coughlan 22.12.2014
# Read LOFAR bbs file and create cluster file clustering sources within a user-specified radius
# Note: Will not distinguish between positive and negative declination

import sys
from math import sqrt
from scipy.cluster.vq import kmeans,vq

min_length = 3	# ignore entry lines shorter than this
chunk_size = 1
bad_axis_cutoff = 1e-6

name_col = 0
type_col = 1
ra_col = 2
dec_col = 3
flux_col = 4
major_axis_col = 8
minor_axis_col = 9
bpa_col = 10


# function to convert RA and DEC strings into degree floats

def ra_to_degrees(ra):
	hours = float(ra[1:3])	# account for space
	minutes = float(ra[4:6])
	seconds = float(ra[7:14])
	return(hours*15 + minutes*0.25 + seconds/240.0)

def dec_to_degrees(dec):
	degrees = float(dec[2:4])	# account for space and minus sign
	minutes = float(dec[5:7])
	seconds = float(dec[8:15])
	return(minutes/60.0 + degrees + seconds/3600.0)

# function to find the separtion between two points, given their RA and DEC in strings

def separation_degrees(point1, point2):
	ra1 = ra_to_degrees(point1[ra_col])
	dec1 = dec_to_degrees(point1[dec_col])

	ra2 = ra_to_degrees(point2[ra_col])
	dec2 = dec_to_degrees(point2[dec_col])

	return( sqrt(pow( ra1 - ra2 ,2) + pow( dec1 - dec2 ,2)) )
	

# Check arguments

if(len(sys.argv)!=5):
	print "\tError: Takes four arguments."
	print "\tUseage: gen_cluster <filename> <flux cutoff> <radius> <outputname>"
	print "\tAllowed radius (degrees)"
	sys.exit()
else:
	inputname = str(sys.argv[1])
	cutoff = float(sys.argv[2])
	radius = float(sys.argv[3])
	outputname = str(sys.argv[4])
	print '\tReading:', inputname
	print '\tWriting:', outputname
	print '\tClipping at:', str(cutoff)
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


nlines = len(lines)
ctr = 0
nodata=True
while(nodata):
	if( len(lines[ctr])<min_length or (lines[ctr][0]=='#' or lines[ctr][0]==',' )):
		ctr = ctr + 1
	else:
		if( lines[ctr].split(" ")[0]=='format' ):
			ctr = ctr + 1
		else:
			nodata=False


# Split by comma to find flux, if flux > cutoff flux, include in output
# close output file

# Set up first cluster
cluster_index=[]
cluster_index.append([ctr])	# Add first part of model to first cluster


for ctr in range(ctr+1,nlines):
	if len(lines[ctr])<min_length:
		break
	string = lines[ctr].split(',')
	if (float(string[flux_col]) > cutoff):
		cluster_found = False
		for cluster_ctr in range(len(cluster_index)):
			if ( separation_degrees(string, lines[cluster_index[cluster_ctr][0]].split(',')) < radius):
				cluster_index[cluster_ctr].append(ctr)
				cluster_found = True
				break
		if cluster_found==False:
			cluster_index.append([ctr])

# Write out additional source if requested (only DGTAU supported right now)

print "Number of clusters found = "+str(len(cluster_index))

try:
	f = open(outputname,'w')
except:
	print "\tError opening ", outputname
	sys.exit()

nbad = 0

for cluster_ctr in range(len(cluster_index)):
	f.write(str(cluster_ctr)+" "+str(chunk_size))
	for ctr in range(len(cluster_index[cluster_ctr])):
		if(lines[cluster_index[cluster_ctr][ctr]].split(',')[type_col] == ' GAUSSIAN'):
			if((float(lines[cluster_index[cluster_ctr][ctr]].split(',')[major_axis_col]) > bad_axis_cutoff) and (float(lines[cluster_index[cluster_ctr][ctr]].split(',')[minor_axis_col]) > bad_axis_cutoff)):
				f.write(" G"+str(lines[cluster_index[cluster_ctr][ctr]].split(',')[name_col]))
			else:
				print "Bad Gaussian detected : "+str(lines[cluster_index[cluster_ctr][ctr]].split(',')[name_col])
				nbad = nbad + 1
		else:
			f.write(" P"+str(lines[cluster_index[cluster_ctr][ctr]].split(',')[name_col]))
	f.write("\n")
f.close()

print str(nbad)+" bad Gaussians ignored."
