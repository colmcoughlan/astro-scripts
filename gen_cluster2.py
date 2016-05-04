#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Colm Coughlan 9.3.2015
# Improved version of gen_cluster
# Automatically scans format line
# Runs an actual clustering algorithm
# Strangely, no increase in effeciency noted over first version (similar clusters per max separation)
# Should have better clustering though

import sys
import numpy as np
from scipy.cluster.vq import vq, kmeans2, whiten
from pylab import plot,show

min_length = 3	# ignore entry lines shorter than this
chunk_size = 1	# SAGECal chunk size (hybrid solution intervals only)


# function to convert RA and DEC strings into degree floats

def gen_coords(data, ra_col, dec_col):
	coords = np.zeros((len(data),2),dtype=float)
	for i in range(len(data)):
		ra = data[i,ra_col].split(':')
		coords[i,0] = float(ra[0])*15 + float(ra[1])*0.25 + float(ra[2])/240.0
		dec = data[i,dec_col].split('.',2)
		coords[i,1] = float(dec[0]) + float(dec[1])/60.0 + float(dec[2])/3600.0
	return(coords)

	

# Check arguments

if(len(sys.argv)!=5):
	print("\tError: Takes 4 arguments.")
	print("\tUseage: gen_cluster <filename> <nclusters> <outputname> <doplot>")
	print("\tImportant notice: Initial guesses for cluster centroids are partly random")
	print("\tPerformance may vary slightly - experiment with kmeans2 if interested (parameter minit)")
	sys.exit()
else:
	inputname = str(sys.argv[1])
	nclusters = int(sys.argv[2])
	outputname = str(sys.argv[3])
	doplot = int(sys.argv[4])
	print("\tReading: "+inputname)
	print("\tWriting: "+outputname)
	print("\tNumber of clusters: "+str(nclusters))


# Attempt to open input file

try:
	f = open(inputname)
except:
	print("\tError opening "+inputname)
	sys.exit()


# Read format string and get rid of any blank lines or comments

ctr = 0
nodata=True
while(nodata):
	line = f.readline()
	# Check for blank lines or comments
	if( len(line)<min_length or line[0]=='#' ):
		ctr = ctr + 1
	else:
		# read format string
		if( line.split(" ")[0]=='format' ):
			formatstr = line.split(" = ")[1].split(", ")
			print("\tDetected format = "+str(formatstr))
			for i in range(len(formatstr)):
				if formatstr[i]=='Name':
					name_col = i
				if formatstr[i]=='Type':
					type_col = i
				if formatstr[i]=='Ra':
					ra_col = i
				if formatstr[i]=='Dec':
					dec_col = i
				if formatstr[i]=='I':
					flux_col = i
				if formatstr[i]=='Q':
					q_col = i
				if formatstr[i]=='U':
					u_col = i
				if formatstr[i]=='V':
					v_col = i
				if formatstr[i]=='MajorAxis':
					major_axis_col = i
				if formatstr[i]=='MinorAxis':
					minor_axis_col = i
				if formatstr[i]=='Orientation':
					bpa_col = i
				if formatstr[i].split("=")[0]=='ReferenceFrequency':
					ref_freq = float(formatstr[i].split("\'")[1])
					ref_freq_col = i
				if formatstr[i].split("=")[0]=='SpectralIndex':
					spectral_col = i
			ctr = ctr + 1
		else:
			# presume we have reached the data
			nodata=False
			f.close()

print('\tReference frequency = '+str(ref_freq)+' Hz')

# Read in data as strings
data = np.genfromtxt(inputname,delimiter=', ',dtype=str,skip_header = ctr-1)

# Convert to coords in degrees
coords = gen_coords(data,ra_col,dec_col)

# Set up for clustering (normalize std dev)
whitened = whiten(coords)

#for i in range(len(data)):
#	print "at "+data[i,ra_col]+" "+data[i,dec_col]+' = '+str(coords[i,0])+" "+str(coords[i,1])+' = '+str(whitened[i,0])+" "+str(whitened[i,1])

# Clustering (note default is to start with random locations, minit = 'random')
# Available methods are ‘random’, ‘points’, ‘uniform’, and ‘matrix’:
#	‘random’: generate k centroids from a Gaussian with mean and variance estimated from the data.
#	‘points’: choose k observations (rows) at random from data for the initial centroids.
#	‘uniform’: generate k observations from the data from a uniform distribution defined by the data set (unsupported).
#	‘matrix’: interpret the k parameter as a k by M (or length k array for one-dimensional data) array of initial centroids.

centroids, closest_centroid = kmeans2(whitened,nclusters,minit='points')
indices, distance = vq(whitened,centroids)

max_sep = np.max(distance)

print("\tMaximum separation between centroid and source = "+str(max_sep)+" degrees.")

# Plot
if doplot==1:
	c1,c2,c3 = 0,0,0
	for i in range(nclusters):
		plot(whitened[indices==i,0],whitened[indices==i,1],'o',color = [c1,c2,c3])
		c1 = c1 + 0.1
		if c1>1:
			c1 = 0
		c2 = 1-c1
		c3 = (c1+c2)/2
	plot(centroids[:,0],centroids[:,1],'sg',markersize=8)
	show()

# Find cluster closest to centre (to apply correction from this cluster to residuals if required)
cen_coords = np.zeros((1,2))
ra = '04:21:59.4'.split(':')
cen_coords[0][0] = float(ra[0])*15 + float(ra[1])*0.25 + float(ra[2])/240.0
dec = '+19.32.06'.split('.',3)
cen_coords[0][1] = float(dec[0]) + float(dec[1])/60.0 + float(dec[2])/360.
closest_source=np.argmin(np.sum(np.abs(np.subtract(cen_coords,coords))**2,axis=-1)**(1./2))

print('\tCentre coords (TTAU) are closest to source '+str(closest_source)+' in cluster ID '+str(indices[closest_source]))

# Write out

try:
	f = open(outputname,'w')
except:
	print "\tError opening ", outputname
	sys.exit()

ngauss = 0
npoints = 0
flux = 0.0

f.write('# Cluster file generated by gen_cluster2.py using '+inputname+'.\n')
f.write('# Total number of clusters = '+str(nclusters)+'.\n')
f.write('# Max separation between cluster centroid and source = '+str(max_sep)+" degrees.\n")
f.write('# Residual cluster ID = '+str(indices[closest_source])+' for specified source with RA = '+str(ra)+', DEC = '+str(dec)+'.\n')
for i in range(nclusters):
	f.write(str(i)+" "+str(chunk_size))
	for j in range(len(data[indices==i])):
		if data[indices==i][j][type_col]=='GAUSSIAN':
			f.write(' G'+data[indices==i][j][name_col])
			ngauss = ngauss + 1
			flux = flux + float(data[indices==i][j][flux_col])
		else:
			f.write(' P'+data[indices==i][j][name_col])
			npoints = npoints + 1
			flux = flux + float(data[indices==i][j][flux_col])
	f.write('\n')
f.close()
print("\t"+str(ngauss)+" Gaussians and "+str(npoints)+" point sources with a total model flux of "+str(flux)+'Jy assigned to clusters.')

