#! /usr/bin/env python
# Plot clusters and give good imaging settings for a particular one
# Colm Coughlan. 18.8.15

import sys
import numpy as np

min_length = 3	# ignore BBS skymodel lines shorter than this

# function to convert RA and DEC strings into degree floats

def gen_coords(data, ra_col, dec_col):
	ra = data[ra_col].split(':')
	ra_deg = float(ra[0])*15 + float(ra[1])*0.25 + float(ra[2])/240.0
	dec = data[dec_col].split('.',2)
	dec_deg = float(dec[0]) + float(dec[1])/60.0 + float(dec[2])/360.
	return(ra_deg , dec_deg)

def deg_to_radec(data):
	ra = float(ra[0])*15 + float(ra[1])*0.25 + float(ra[2])/240.0
	dec = float(dec[0]) + float(dec[1])/60.0 + float(dec[2])/3600.0
	return(ra_deg , dec_deg)

def get_format(filename, format_dict):
	# Attempt to open input file

	try:
		f = open(filename)
	except:
		print("\t Error opening "+filename)
		sys.exit()


	# Read format string and get rid of any blank lines or comments

	ctr = 0
	nodata=True
	key_list = format_dict.keys()
	while(nodata):
		line = f.readline()
		# Check for blank lines or comments
		if( len(line)<min_length or line[0]=='#' ):
			ctr = ctr + 1
		else:
			# read format string
			if( line.split(" ")[0]=='format' ):
				saved_format = line
				formatstr = (line.rstrip('\n')).split(" = ")[1].split(", ")
				print("\t Detected format line from "+filename+" = "+str(formatstr))
				for i in range(len(formatstr)):
					for j in range(len(key_list)):
						if formatstr[i]==key_list[j]:
							format_dict[ key_list[j] ] = i
						else:
							if formatstr[i].split("=")[0]==key_list[j]:
								format_dict[ key_list[j] ] = i
				ctr = ctr + 1
			else:
				# presume we have reached the data
				nodata=False
				f.close()
	return(format_dict , ctr)


def get_skymodel_data(source , skymodel_data, format_dict):
	for i in range(len(skymodel_data)):
		if(skymodel_data[i][format_dict['Name']]==source):
			return(skymodel_data[i])
	return(1)


# Start of code

if(len(sys.argv)!=5):
	print("\t Error: Takes 4 arguments.")
	print("\t Useage: cluster_plotter <BBS skymodel> <SAGEcal cluster file> <cluster id> <doplot>")
	sys.exit()
else:
	skymodel_file = str(sys.argv[1])
	cluster_file = str(sys.argv[2])
	cluster_id = int(sys.argv[3])
	doplot = int(sys.argv[4])
	print("\t Reading skymodel (BBS format): "+skymodel_file)
	print("\t Reading cluster file: "+cluster_file)

# Obtain the column numbers corresponding to different quantities from the format line of the BBS file

key_list = ['Name' , 'Type', 'Ra', 'Dec', 'I', 'Q', 'U', 'V', 'MajorAxis', 'MinorAxis', 'Orientation', 'ReferenceFrequency', 'SpectralIndex']
num_keys = len(key_list)
init_vals = [-1] * num_keys
format_dict = dict(zip(key_list,init_vals))

format_dict , i = get_format(skymodel_file, format_dict)
skymodel_data = np.genfromtxt(skymodel_file,delimiter=', ',dtype=str,skip_header = i-1)

# Read cluster information from file

cluster_found = False
with open(cluster_file) as f:
	for line in f:
		# Check for blank lines or comments
		if( len(line)>min_length and line[0]!='#' ):
			if( line.split(" ")[0]==str(cluster_id) ):
				cluster_line = line.split("\n")[0]	# remove newline
				cluster_line = cluster_line.split(" ")
				cluster_found = True

if( not(cluster_found) ):
	print("\tError: cluster ID "+str(cluster_id)+" not found in "+cluster_file)
	sys.exit()

offset = 2	#offset for source name in cluster file (ignore cluster number and integration fraction)
nsources = len(cluster_line)-offset
print('\t'+str(nsources)+' sources detected')

coords = np.zeros((nsources,2))

for i in range(nsources):

#search data for string (remember to remove P, S or G character from LSM -> BBS format)

	newdata = get_skymodel_data( cluster_line[offset+i][1:] , skymodel_data , format_dict)
	if(newdata == 1):
		print("\tError: Cannot find source "+cluster_line[offset+i][1:]+" from cluster ID "+str(cluster_id)+" in "+skymodel_file)
		sys.exit()
	coords[i] = gen_coords(newdata, format_dict['Ra'], format_dict['Dec'])


print(str(coords))
cent = [np.mean(coords[:,0]),np.mean(coords[:,1])]
print(str(cent))

