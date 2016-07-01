#! /usr/bin/env python
# Colm Coughlan 20.8.2015
# Dublin Institute for Advanced Studies

import numpy as np
import scipy.spatial
import pandas as pd
import argparse

# Function to get frequency in Hz
def get_ref_freq(inputname):
	# Attempt to open input file

	try:
		f = open(inputname)
	except:
		print("\t Error opening "+inputname)
		exit()

	f.readline()	# Burn the first two lines
	f.readline()
	return(float(f.readline().split(" ")[8]))	# This assumes the very particular structure of the output from PyBDSM (csv mode)

def deg_to_sexagesimal(df,delim,doround=-1):
	coords = []

	rem_ra , hrs = np.modf(np.divide(df['RA'],15.0))
	seconds_ra , minutes_ra = np.modf(np.multiply(rem_ra,60.0))
	seconds_ra = np.multiply(seconds_ra,60.0)

	rem_dec , deg = np.modf(df['DEC'])
	seconds_dec , minutes_dec = np.modf(np.multiply(rem_dec,60.0))
	seconds_dec = np.multiply(seconds_dec,60.0)

	if(doround>0):	# doround specifies number of degrees of precision in DEC. Use one more in RA (as it is less precise).
		if doround==1:
			prec = ['%05.2f','%04.1f']
		if doround==2:
			prec = ['%06.3f','%05.2f']
	else:
		prec = ['%09f','%09f']

	for i in range(len(df['RA'])):
		if((prec[0]%seconds_ra.iloc[i]) == prec[0]%60):	# watch out for rounding 60sec, 60min etc.
			seconds_ra.iloc[i] = 0.0
			minutes_ra.iloc[i] = minutes_ra.iloc[i] + 1
		if((prec[0]%minutes_ra.iloc[i]) == prec[0]%60):
			minutes_ra.iloc[i] = 0.0
			hrs.iloc[i] = hrs.iloc[i] + 1
		if((prec[1]%seconds_dec.iloc[i]) == prec[1]%60):
			seconds_dec.iloc[i] = 0.0
			minutes_dec.iloc[i] = minutes_dec.iloc[i] + 1
		if((prec[1]%minutes_dec.iloc[i]) == prec[1]%60):
			minutes_dec.iloc[i] = 0.0
			deg.iloc[i] = deg.iloc[i] + 1
		coords.append(( ('%02d:%02d:'+prec[0])%(hrs.iloc[i],minutes_ra.iloc[i],seconds_ra.iloc[i]))+delim+('+%02d:%02d:'+prec[1])%(deg.iloc[i],minutes_dec.iloc[i],seconds_dec.iloc[i]))

	return(coords)

def sexagesimal_to_deg(ra_list, dec_list):
	
	ra = (float(ra_list[0]) + float(ra_list[1])/60.0 +float(ra_list[2])/3600.0)*15.0
	dec = (float(dec_list[0]) + float(dec_list[1])/60.0 +float(dec_list[2])/3600.0)
	
	return(ra, dec)
	
	# Calculte the spectral index and error as normal
def gen_spec_indx(flux1, flux2, flux1_err, flux2_err, freq1, freq2):
	freq_ratio = np.log10(freq2/freq1)
	spec_indx = np.divide(np.log10( np.divide(flux2 , flux1 ) ) , freq_ratio)
	spec_indx_err = np.sqrt( (flux1_err / flux1 )**2 + (flux2_err / flux2 )**2 ) / np.abs(freq_ratio)
	return(spec_indx , spec_indx_err)
	
	# Make a list of catalog IDs, given dataframes at each freq and the matching relation between the two (to give duplicate IDs)
def gen_cat_id(df1, df2):
	df1['Cat_ID'] = -1
	df2['Cat_ID'] = -1
	
	# Get IDs for the lower frequency. Copy them to matches in the second freq.
	id = 0
	for i in df1.index:
		df1.loc[i,'Cat_ID'] = id
		if( df1.loc[i,'Match'] == True):
			df2.loc[ df1.loc[i,'nearest_df2_index'],'Cat_ID'] = id
		id = id + 1
		
	# Finish the IDs in the second freq.	
	for i in df2.index:
		if( df2.loc[i,'Cat_ID'] == -1):
			df2.loc[i,'Cat_ID'] = id
			id = id + 1
	return(df1['Cat_ID'], df2['Cat_ID'])
	
	
	# Re-assign IDs after merging two catalogs, watching out for gaussian's belonging to the same source
def resort_cat_id(df):
	size = len(df)	
	cat_id = df['Cat_ID'].values
	cat_id_new = np.zeros((size,1),dtype=np.int)
	
	for i in range(size):
		cat_id_new[i] = i
		if(i +1 < size):
			if(cat_id[i] == cat_id[i+1]):
				i = i + 1
				cat_id_new[i] = i

	return(cat_id_new)
				

	
def make_kvis_file(filename, df):
	n = len(df)
	done = np.zeros((n,1),dtype=np.bool)
	
	with open(filename+'.ann', 'w') as f:
		f.write('color red\n')
		
		# make sure not to write out high and low detections for the same source. Write out high freq position if there is a match
		for i in range(n):
			if(~done[i]):
				done[i] = True
				j = i
				match = False
				if(dfc['Match'].iloc[i]):
					if(i +1 != n):
						for k in range(i+1,n):
							if(dfc['Cat_ID'].iloc[j] == dfc['Cat_ID'].iloc[k]):	# if there is a match, report the 610 data as k
								done[k] = True
								match = True
								if(dfc['Freq'].iloc[j] > dfc['Freq'].iloc[k]):
									temp = j
									j = k		# set k to be the high freq position, j to be the low freq
									k = temp
								break
							else:
								match = False
				if(match):
					f.write('CROSS W '+str(df['RA'].iloc[k])+' '+str(df['DEC'].iloc[k])+' '+str(df['E_RA'].iloc[k])+' '+str(df['E_DEC'].iloc[k])+' 90.0\n')
				else:
					f.write('CROSS W '+str(df['RA'].iloc[j])+' '+str(df['DEC'].iloc[j])+' '+str(df['E_RA'].iloc[j])+' '+str(df['E_DEC'].iloc[j])+' 90.0\n')
			
# Print out an entry in the latex table file			

def print_entry(entry,real):
	if(real):
		f.write('$'+'%.2f'%entry['Peak_flux'].values[0]+'\pm'+'%.2f'%entry['E_Peak_flux'].values[0]+'$&')
		f.write('$'+'%.2f'%entry['Total_flux'].values[0]+'\pm'+'%.2f'%entry['E_Total_flux'].values[0]+'$&')
		f.write('$'+str(entry['Resid_Isl_rms'].values[0])+'$&')		
		'''
		if(entry['Resolved'].values[0]):
			if(entry['S_Code'].values!='M'):
				f.write('$'+str(entry['DC_Maj'].values[0])+'\pm'+str(entry['DC_E_Maj'].values[0])+'$&')
				f.write('$'+str(entry['DC_Min'].values[0])+'\pm'+str(entry['DC_E_Min'].values[0])+'$&')
				f.write('$'+str(entry['DC_PA'].values[0])+'\pm'+str(entry['DC_E_PA'].values[0])+'$&')
			else:
				f.write('-&-&-&')
		else:
			f.write('-&-&-&')
		if(entry['Resolved'].values[0]):
			f.write('T&')
		else:
			f.write('F&')
		'''
		f.write(str(entry['S_Code'].values[0])+'&')
	else:
#		f.write('-&-&-&-&-&-&-&')
		f.write('-&-&-&-&')


# Write out a list of coords in the NVSS format (RA, DEC, radius in as)
def write_nvss(filename, list, dfc, radius):
	n = len(dfc)
	done = np.zeros((n,1),dtype=np.bool)
	list = [l.replace(':', ' ') for l in list]
	
	with open(filename, 'w') as f:
		if radius < 0.0:
			for i in range(n):
				f.write(str(list[i])+'\n')		# if a negative radius is given, ignore it
		else:
			for i in range(n):
				f.write(str(list[i])+' '+str(radius)+'\n')	

def read_nvss(filename, n):
	have_detection = np.zeros((n,1),dtype=np.bool)
	flux = np.zeros((n,1),dtype=np.float)
	flux_error = np.zeros((n,1),dtype=np.float)
	ra = np.zeros((n,1),dtype=np.float)
	dec = np.zeros((n,1),dtype=np.float)
	# note nmatches is the number of matches within the data (not to NVSS)

	i = 0
	n_matched = 0
	n_unmatched = 0
	trigger = False
	grab_next_line = False
	grab_error_line = False
	with open(filename, 'r') as f:
		for line in f:
			if(grab_error_line):
				if( ((line.find('NVSS')<0) and (line.find('RA')<0)) and ((line.find('deg')<0) and (len(line) > 1))):
					error_line = line
					grab_error_line = False
		
			if(grab_next_line):		# grab the source line following a :, watch out for new pages
				if( ((line.find('NVSS')<0) and (line.find('RA')<0)) and ((line.find('deg')<0) and (len(line) > 1))):
					source_line = line
					grab_next_line = False
					grab_error_line = True
				
			if(line[0]==':'):	# a colon indicates a search for a position in the file. Reaching two colons, without an "S" as the first letter
				grab_next_line = True
				if(trigger):	# is an indication of a detection. The trigger variable is set to detect this
					have_detection[i] = True
					split_line = source_line.split()
					flux[i] = float(split_line[7])
					flux_error[i] = float((error_line.split())[3])
					ra[i], dec[i] = sexagesimal_to_deg(split_line[0:3], split_line[3:6])
					i = i + 1
					n_matched = n_matched + 1
				else:
					trigger = True
					
					
			if(line[0]=='S'):
				have_detection[i] = False	# "SOURCE NOT FOUND" -> no detection at current source. reset trigger to find start of next source
				trigger = False
				i = i + 1
				n_unmatched = n_unmatched + 1
				
				
	if(trigger):
		have_detection[i] = True	# just in case there is a detection at the last entry
		split_line = source_line.split()
		flux[i] = float(split_line[7])
		flux_error[i] = float((error_line.split())[3])
		ra[i], dec[i] = sexagesimal_to_deg(split_line[0:3], split_line[3:6])
		i = i + 1
		n_matched = n_matched + 1
				
	if(n_matched + n_unmatched != n):
		print('Error in NVSS scanning!. NVSS column incorrect!')
		print('Found '+str(n_matched)+' matches')
		print('Found '+str(n_unmatched)+' unmatched')
		print('Giving a total of '+str(n_matched + n_unmatched))
		print('But the total number of sources processed should be '+str(n))

	return(have_detection, flux, flux_error, ra, dec)
	
def read_vizier(cat_name, nskip, col_list):
	df = pd.read_csv(cat_name,skiprows=nskip,delimiter='\t', engine='python', comment='#')
	return(df[col_list])

###############################
#
#	Main code starts here
#
###############################

# Check arguments

parser = argparse.ArgumentParser(description='Colm Coughlan. Dublin Institute for Advanced Studies. Make a survey map from 2 PyBDSM source lists.')
parser.add_argument('input_cat_1', type=str, help='Lower frequency source list filename')
parser.add_argument('input_cat_2', type=str, help='Higher frequency source list filename')
parser.add_argument('phase_centre', type=str, help='Phase centre. Form: \'12 34 56.78+12 34 56.78\'')
parser.add_argument('radius', type=float, help='Maximum separation assumed to be in as.')
parser.add_argument('output_stem', type=str, help='Stem for output files.')
parser.add_argument('--nvss', type=str, help='NVSS detection printout.')
parser.add_argument('--xest', type=str, help='XEST vizier file.')
parser.add_argument('--mass', type=str, help='2MASS vizier file.')
parser.add_argument('--spitzer', type=str, help='Spitzer C2E vizier file.')
parser.add_argument('--gbs', type=str, help='Gould Belt survey vizier file.')
parser.add_argument('--gbs_counterparts', type=str, help='Gould Belt survey contourparts vizier file.')
parser.add_argument('--aclass', type=str, help='Additional classification list.')

args = parser.parse_args()


inputname1 = args.input_cat_1
inputname2 = args.input_cat_2
radius = args.radius/3600.
outputname = args.output_stem

phase_centre = np.zeros((2,1))
if((args.phase_centre).find('+') > 0):
	pc_list = (args.phase_centre).split('+')
	ra_list = pc_list[0].split(' ')
	dec_list = pc_list[1].split(' ')
	phase_centre[0], phase_centre[1] = sexagesimal_to_deg(ra_list, dec_list)
#	phase_centre[0] = (float(ra_list[0]) + float(ra_list[1])/60.0 +float(ra_list[2])/3600.0)*15.0
#	phase_centre[1] = (float(dec_list[0]) + float(dec_list[1])/60.0 +float(dec_list[2])/3600.0)
else:
	if((args.phase_centre).find('-') > 0):
		pc_list = (args.phase_centre).split('-')
		ra_list = pc_list[0].split(' ')
		dec_list = pc_list[1].split(' ')
		phase_centre[0] = (float(ra_list[0]) + float(ra_list[1])/60.0 +float(ra_list[2])/3600.0)*15.0
		phase_centre[1] = -(float(dec_list[0]) + float(dec_list[1])/60.0 +float(dec_list[2])/3600.0)
	else:
		print('Error in phase center. Format should be : \'1234+4567\'')
		exit()
		
# Check for for optional arguments

# NVSS

if(args.nvss is not None):
	have_nvss = True
	nvss_file = args.nvss
else:
	have_nvss = False
	
# xmm (XEST)

if(args.xest is not None):
	have_xmm = True
	xmm_file = args.xest
else:
	have_xmm = False
	
# 2MASS

if(args.mass is not None):
	have_mass = True
	mass_file = args.mass
else:
	have_mass = False
	
# spitzer

if(args.spitzer is not None):
	have_spitzer = True
	spitzer_file = args.spitzer
else:
	have_spitzer = False
	
# GBS
	
if(args.gbs is not None):
	have_gbs = True
	gbs_file = args.gbs
else:
	have_gbs = False
	
# GBS counterparts
	
if(args.gbs_counterparts is not None):
	have_gbs_counterparts = True
	gbs_counterparts_file = args.gbs_counterparts
else:
	have_gbs_counterparts = False
	
# Classification file
	
if(args.aclass is not None):
	have_class = True
	class_list = args.aclass
else:
	have_class = False
print('\t Reading: '+inputname1+', '+inputname2)
print('\t Writing: '+outputname)
if have_nvss:
	print('\t Using '+nvss_file+' for NVSS detection information.')
if have_class:
	print('\t Using '+class_list+' for additional classification information.')

	
# Set the primary and synthesised beam information

pbeam1 = 81.0/60.0	# in degrees, low freq
pbeam2 = 43.0/60.0

	
beam1 = 10.0	# in arcsec
beam2 = 10.0
res_factor = 1.0

nvss_freq = 1.4E9	# in Hz
nvss_comp_radius = 5.0 # in as
gbs_radius = 5.0/3600.0	# comparison radius in degrees


# use pandas to read csv file. Only out of date pandas easily available on Ubuntu 12.04
# Assumed format: .csv file from PyBDSM

# Get the reference frequencies of each file

freq1 = get_ref_freq(inputname1)
freq2 = get_ref_freq(inputname2)

if freq1 > freq2:
	print('Please enter the lower frequency first')
	exit()


# Then read in column names (different delimiter...). Assume both files have the same structure

df1 = pd.read_csv(inputname1,skiprows=5,delimiter=' ',nrows=1, engine='python')

column_names = df1.columns.values.tolist()
column_names.pop(0)	# Get rid of "#" at the start of the file

# Now read in entire files, assuming a whitespace delimiter. Do a sort on Ra + Dec. Index by Source_id

df1 = (pd.read_csv(inputname1,skiprows=6,delimiter='\s+',names = column_names, skipinitialspace=True, engine='python', index_col='Source_id')).sort_values(by=['RA','DEC'],ascending=[1,1])
df2 = (pd.read_csv(inputname2,skiprows=6,delimiter='\s+',names = column_names, skipinitialspace=True, engine='python', index_col='Source_id')).sort_values(by=['RA','DEC'],ascending=[1,1])

# Drop sources outside the primary beam

distances = scipy.spatial.distance.cdist(np.array([df1['RA'],df1['DEC']]).T , np.array([phase_centre[0],phase_centre[1]]).T,'euclidean')
df1 = df1[distances < pbeam1]

distances = scipy.spatial.distance.cdist(np.array([df2['RA'],df2['DEC']]).T , np.array([phase_centre[0],phase_centre[1]]).T,'euclidean')
df2 = df2[distances < pbeam2]

if( (len(df1)==0) or (len(df2)==0) ):
	print('Error, no sources inside primary beam at at least one frequency')
	exit()

# Change units to mJy, arcsec etc. Make sure to change errors too

df1['Total_flux'] = df1['Total_flux']*1000.0
df1['E_Total_flux'] = df1['E_Total_flux']*1000.0
df1['Peak_flux'] = df1['Peak_flux']*1000.0
df1['E_Peak_flux'] = df1['E_Peak_flux']*1000.0
df1['Maj'] = df1['Maj']*3600.0
df1['E_Maj'] = df1['E_Maj']*3600.0
df1['Min'] = df1['Min']*3600.0
df1['E_Min'] = df1['E_Min']*3600.0
df1['DC_Maj'] = df1['Maj']*3600.0
df1['E_DC_Maj'] = df1['E_Maj']*3600.0
df1['DC_Min'] = df1['Min']*3600.0
df1['E_DC_Min'] = df1['E_Min']*3600.0

df2['Total_flux'] = df2['Total_flux']*1000.0
df2['E_Total_flux'] = df2['E_Total_flux']*1000.0
df2['Peak_flux'] = df2['Peak_flux']*1000.0
df2['E_Peak_flux'] = df2['E_Peak_flux']*1000.0
df2['Maj'] = df2['Maj']*3600.0
df2['E_Maj'] = df2['E_Maj']*3600.0
df2['Min'] = df2['Min']*3600.0
df2['E_Min'] = df2['E_Min']*3600.0
df2['DC_Maj'] = df2['Maj']*3600.0
df2['E_DC_Maj'] = df2['E_Maj']*3600.0
df2['DC_Min'] = df2['Min']*3600.0
df2['E_DC_Min'] = df2['E_Min']*3600.0


# Get better flux errors first - using method from Rachael's thesis without the additional rms term

df1['E_Total_flux'] = np.sqrt(df1['E_Total_flux'].values**2 + (0.05*df1['Total_flux'].values)**2)
df1['E_Peak_flux'] = np.sqrt(df1['E_Peak_flux'].values**2 + (0.05*df1['Peak_flux'].values)**2)

df2['E_Total_flux'] = np.sqrt(df2['E_Total_flux'].values**2 + (0.05*df2['Total_flux'].values)**2)
df2['E_Peak_flux'] = np.sqrt(df2['E_Peak_flux'].values**2 + (0.05*df2['Peak_flux'].values)**2)


# Insert frequency columns

df1['Freq'] = freq1
df2['Freq'] = freq2

# Classify type S and C Gaussians as resolved or not. Note all type M sources are by definition resolved

df1['Resolved'] = np.greater(df1['Maj'].values , beam1*res_factor)
df2['Resolved'] = np.greater(df2['Maj'].values , beam2*res_factor)
df1.loc[df1['S_Code']=='M','Resolved'] = True
df2.loc[df2['S_Code']=='M','Resolved'] = True

# Now find the distances between each source (N.B. Class M sources merged already)

distances = scipy.spatial.distance.cdist(np.array([df1['RA'],df1['DEC']]).T , np.array([df2['RA'],df2['DEC']]).T,'euclidean')
min_dist = np.amin(distances,axis=1)
df1['nearest_df2_index'] = df2.index[np.argmin(distances,axis=1)] # Would fail if the entire array was NaN

# Identify sources across frequencies

df1['Match'] = False
df2['Match'] = False

df1.loc[min_dist < radius, 'Match'] = True	# if min_distance < radius, we have a match
df2.loc[df1.loc[df1['Match'], 'nearest_df2_index'], 'Match'] = True # set the corresponding source at freq2 to True also
nmatches = np.sum(df1['Match'].values)
nmatches_test = np.sum(df2['Match'].values)

if(nmatches!=nmatches_test):
	print('Uneven number of matches!')
else:
	print(str(nmatches)+' matches detected.')
	print('Mean separation = '+str(3600.0*np.mean(min_dist[min_dist < radius]))+' as.')

'''
# Sanity check for matching
print('RA1 = '+str(df1.loc[df1['Match'], 'RA'].values[0]))
print('DEC1 = '+str(df1.loc[df1['Match'], 'DEC'].values[0]))
print('RA2 = '+str(df2.loc[df1.loc[df1['Match'], 'nearest_df2_index'], 'RA'].values[0]))
print('DEC2 = '+str(df2.loc[df1.loc[df1['Match'], 'nearest_df2_index'], 'DEC'].values[0]))
print(radius)
'''

# Find Spectral indices if sources are close enough to be considered the same. Assume freq2 > freq1
# Use total flux if resolved, and peak flux otherwise

df1['SI'] = -999.0
df2['SI'] = -999.0
df1['E_SI'] = -999.0
df2['E_SI'] = -999.0

# Collect fluxes and indices of sources in df1 with a match in df2

f1_tf = df1.loc[df1['Match'],'Total_flux'].values
f1_pf = df1.loc[df1['Match'],'Peak_flux'].values
f1_e_tf = df1.loc[df1['Match'],'E_Total_flux'].values
f1_e_pf = df1.loc[df1['Match'],'E_Peak_flux'].values
f1_index, f1_res = df1.loc[df1['Match'],'Resolved'].index , df1.loc[df1['Match'],'Resolved'].values

# The corresponding df2 information

f2_tf = df2.loc[df1.loc[f1_index, 'nearest_df2_index'], 'Total_flux'].values
f2_pf = df2.loc[df1.loc[f1_index, 'nearest_df2_index'], 'Peak_flux'].values
f2_e_tf = df2.loc[df1.loc[f1_index, 'nearest_df2_index'], 'E_Total_flux'].values
f2_e_pf = df2.loc[df1.loc[f1_index, 'nearest_df2_index'], 'E_Peak_flux'].values
f2_index, f2_res = df2.loc[df1.loc[f1_index, 'nearest_df2_index'], 'Resolved'].index , df2.loc[df1.loc[f1_index, 'nearest_df2_index'], 'Resolved'].values

# Choose the right flux information (peak if unresolved, total otherwise)

f1_array = np.where(f1_res, f1_tf, f1_pf)
f2_array = np.where(f2_res, f2_tf, f2_pf)
f1_e_array = np.where(f1_res, f1_e_tf, f1_e_pf)
f2_e_array = np.where(f2_res, f2_e_tf, f2_e_pf)

# get the spectral index and errors
spec_indx, spec_indx_err = gen_spec_indx(f1_array, f2_array, f1_e_array, f2_e_array, freq1, freq2)

# save the result back into the dataframes

df1.loc[f1_index, 'SI'] = spec_indx
df2.loc[f2_index, 'SI'] = spec_indx
df1.loc[f1_index, 'E_SI'] = spec_indx_err
df2.loc[f2_index, 'E_SI'] = spec_indx_err



##########################################################################################################
#
#
# NVSS processing
#
#
##########################################################################################################


write_nvss(outputname+'_freq1.crdlst.txt', deg_to_sexagesimal(df1,' ',False), df1, nvss_comp_radius)
write_nvss(outputname+'_freq2.crdlst.txt', deg_to_sexagesimal(df2,' ',False), df2, nvss_comp_radius)
write_nvss(outputname+'_freq1.simbad.txt', deg_to_sexagesimal(df1,' ',False), df1, -1.0)
write_nvss(outputname+'_freq2.simbad.txt', deg_to_sexagesimal(df2,' ',False), df2, -1.0)

if(have_gbs):
	gbs_df = read_vizier(gbs_file, 74,['_RA','_DE','F4.5','Name'])	# read in the GBS ra and decs in degrees
	gbs_df = gbs_df.drop(gbs_df.index[0:2])	# clip units line and ----- line
	gbs_df.reset_index(inplace=True, drop=True)
	gbs_df.reindex(index=range(0,len(gbs_df)))
	
	df1['GBS'] = False
	df1['GBS_Flux'] = 0.0
	distances = scipy.spatial.distance.cdist(np.array([df1['RA'],df1['DEC']]).T , np.array([gbs_df['_RA'], gbs_df['_DE']]).T,'euclidean')
	min_dist = np.amin(distances,axis=1)
	min_dist_index = np.argmin(distances,axis=1)
	df1['GBS'] =  min_dist<gbs_radius
	df1.loc[df1['GBS'], 'GBS_Flux'] = gbs_df.loc[min_dist_index[min_dist<gbs_radius], 'F4.5'].values
	df1.loc[df1['GBS'], 'GBS_Name'] = gbs_df.loc[min_dist_index[min_dist<gbs_radius], 'Name'].values
	
	t1 = (df1.loc[df1['GBS'], 'Total_flux'].values)
	t2 = (df1.loc[df1['GBS'], 'GBS_Flux'].values).astype(np.float)
	print('GBS: Number of possible thermal sources at freq1 = '+str(np.sum(t1<t2)))

	df2['GBS'] = False
	df2['GBS_Flux'] = 0.0	
	distances = scipy.spatial.distance.cdist(np.array([df2['RA'],df2['DEC']]).T , np.array([gbs_df['_RA'], gbs_df['_DE']]).T,'euclidean')
	min_dist = np.amin(distances,axis=1)
	min_dist_index = np.argmin(distances,axis=1)
	df2['GBS'] =  min_dist<gbs_radius
	df2.loc[df2['GBS'], 'GBS_Flux'] = gbs_df.loc[min_dist_index[min_dist<gbs_radius], 'F4.5'].values
	df2.loc[df2['GBS'], 'GBS_Name'] = gbs_df.loc[min_dist_index[min_dist<gbs_radius], 'Name'].values

	
	t1 = (df2.loc[df2['GBS'], 'Total_flux'].values)
	t2 = (df2.loc[df2['GBS'], 'GBS_Flux'].values).astype(np.float)
	print('GBS: Number of possible thermal sources at freq2 = '+str(np.sum(t1<t2)))
	
	print('GBS: Detected '+str(np.sum(df1['GBS'].values))+' matches at freq1, '+str(np.sum(df1['GBS'].values))+' matches at freq2.')

if(have_nvss):
	df1['NVSS'], df1['NVSS_flux'], df1['NVSS_E_flux'], df1['NVSS_RA'], df1['NVSS_DEC'] = read_nvss(nvss_file+'_freq1.NVSS.txt', len(df1))
	df2['NVSS'], df2['NVSS_flux'], df2['NVSS_E_flux'], df2['NVSS_RA'], df2['NVSS_DEC'] = read_nvss(nvss_file+'_freq2.NVSS.txt', len(df2))

	#########################################################################################################
	# Make spectral indices
	#########################################################################################################
	
	# assign number of free parameters to fit. A= low, B = high, C =all three. N = None

	df1['NVSS_SI']= -999.0
	df2['NVSS_SI']= -999.0
	df1['NVSS_E_SI']= -999.0
	df2['NVSS_E_SI']= -999.0
	df1['NVSS_SI_FP']= 'N'
	df2['NVSS_SI_FP']= 'N'
	df1.loc[df1['NVSS'],'NVSS_SI_FP']= 'L'
	df2.loc[df2['NVSS'],'NVSS_SI_FP']= 'H'
	df1.loc[np.all([df1['Match'], df1['NVSS']],axis=0),'NVSS_SI_FP']= 'B'
#	df2.loc[np.all([df2['Match'], df2['NVSS']],axis=0),'NVSS_SI_FP']= 'B'

	# (need to make sure 'B' corresponds to an independent detection at both frequencies. Make sure both match and NVSS conditions are satisfied.
	df2.loc[df1.loc[ df1['NVSS_SI_FP']=='B', 'nearest_df2_index'], 'NVSS_SI_FP'] = 'B'
	bad_indices = df2.loc[ np.all([df2['NVSS'] == False, df2['NVSS_SI_FP']=='B'],axis=0), 'NVSS_SI_FP'].index
	for i in bad_indices:
		f1_index, f1_res = df1.loc[np.all([df1['Match'], df1['NVSS']],axis=0),'nearest_df2_index'].index.tolist(), df1.loc[np.all([df1['Match'], df1['NVSS']],axis=0),'nearest_df2_index'].values
		temp = np.argwhere(f1_res==i)
		if len(temp>0):
			f1_index = f1_index[temp]
			df1.loc[f1_index, 'NVSS_SI_FP'] = 'L'	# downgrade to a single freq match
	df2.loc[df2['NVSS'] == False, 'NVSS_SI_FP'] = 'N'
	
			

	
	if(np.sum(df1['NVSS_SI_FP']=='B') != np.sum(df2['NVSS_SI_FP']=='B')):
		print('Error matching NVSS spectral indices. Unequal number of 3 point fits across both frequencies.')
		print(str(np.sum(df1['NVSS_SI_FP']=='B'))+' != '+str(np.sum(df2['NVSS_SI_FP']=='B')))
		exit()
		

	
	# fit spectral indices and save
	
	# freq1
	f1_tf = df1.loc[df1['NVSS_SI_FP']=='L','Total_flux'].values
	f1_pf = df1.loc[df1['NVSS_SI_FP']=='L','Peak_flux'].values
	f1_e_tf = df1.loc[df1['NVSS_SI_FP']=='L','E_Total_flux'].values
	f1_e_pf = df1.loc[df1['NVSS_SI_FP']=='L','E_Peak_flux'].values
	f1_index, f1_res = df1.loc[df1['NVSS_SI_FP']=='L','Resolved'].index , df1.loc[df1['NVSS_SI_FP']=='L','Resolved'].values
	
	f1_array = np.where(f1_res, f1_tf, f1_pf)
	f1_e_array = np.where(f1_res, f1_e_tf, f1_e_pf)
	
	spec_indx, spec_indx_err = gen_spec_indx(f1_array, df1.loc[df1['NVSS_SI_FP']=='L','NVSS_flux'].values, f1_e_array, df1.loc[df1['NVSS_SI_FP']=='L','NVSS_E_flux'].values, freq1, nvss_freq)
	df1.loc[f1_index, 'NVSS_SI'] = spec_indx
	df1.loc[f1_index, 'NVSS_E_SI'] = spec_indx_err
	
	# freq2 (same again)
	f2_tf = df2.loc[df2['NVSS_SI_FP']=='H','Total_flux'].values
	f2_pf = df2.loc[df2['NVSS_SI_FP']=='H','Peak_flux'].values
	f2_e_tf = df2.loc[df2['NVSS_SI_FP']=='H','E_Total_flux'].values
	f2_e_pf = df2.loc[df2['NVSS_SI_FP']=='H','E_Peak_flux'].values
	f2_index, f2_res = df2.loc[df2['NVSS_SI_FP']=='H','Resolved'].index , df2.loc[df2['NVSS_SI_FP']=='H','Resolved'].values
	
	f2_array = np.where(f2_res, f2_tf, f2_pf)
	f2_e_array = np.where(f2_res, f2_e_tf, f2_e_pf)
	
	spec_indx, spec_indx_err = gen_spec_indx(f2_array, df2.loc[df2['NVSS_SI_FP']=='H','NVSS_flux'].values, f2_e_array, df2.loc[df2['NVSS_SI_FP']=='H','NVSS_E_flux'].values, freq2, nvss_freq)
	df2.loc[f2_index, 'NVSS_SI'] = spec_indx
	df2.loc[f2_index, 'NVSS_E_SI'] = spec_indx_err	
	
	# freq1 and freq2
	
	f1_tf = df1.loc[df1['NVSS_SI_FP']=='B','Total_flux'].values
	f1_pf = df1.loc[df1['NVSS_SI_FP']=='B','Peak_flux'].values
	f1_e_tf = df1.loc[df1['NVSS_SI_FP']=='B','E_Total_flux'].values
	f1_e_pf = df1.loc[df1['NVSS_SI_FP']=='B','E_Peak_flux'].values
	f1_index, f1_res = df1.loc[df1['NVSS_SI_FP']=='B','Resolved'].index , df1.loc[df1['NVSS_SI_FP']=='B','Resolved'].values
	
	f2_tf = df2.loc[df1.loc[f1_index, 'nearest_df2_index'], 'Total_flux'].values
	f2_pf = df2.loc[df1.loc[f1_index, 'nearest_df2_index'], 'Peak_flux'].values
	f2_e_tf = df2.loc[df1.loc[f1_index, 'nearest_df2_index'], 'E_Total_flux'].values
	f2_e_pf = df2.loc[df1.loc[f1_index, 'nearest_df2_index'], 'E_Peak_flux'].values
	f2_index, f2_res = df2.loc[df1.loc[f1_index, 'nearest_df2_index'], 'Resolved'].index , df2.loc[df1.loc[f1_index, 'nearest_df2_index'], 'Resolved'].values
	
	f3_tf = df2.loc[df1.loc[f1_index, 'nearest_df2_index'], 'NVSS_flux'].values
#	print f3_tf	# there shouldn't be a zero here!
#	print (df2['NVSS_SI_FP'].values)[df1.loc[f1_index, 'min_distance_index']]
#	print df2.loc[df2['NVSS_SI_FP']=='B','NVSS_flux'].values
	f3_e_tf = df2.loc[df1.loc[f1_index, 'nearest_df2_index'], 'NVSS_E_flux'].values
	
	f1_array = np.where(f1_res, f1_tf, f1_pf)
	f2_array = np.where(f2_res, f2_tf, f2_pf)
	f1_e_array = np.where(f1_res, f1_e_tf, f1_e_pf)
	f2_e_array = np.where(f2_res, f2_e_tf, f2_e_pf)

	if(len(f1_array)!=len(f2_array)):
		print('Error! NVSS spectral index problem')
		exit()



			
	spec_indx = np.zeros((len(f1_array),1))
	spec_indx_err = np.zeros((len(f1_array),1))
	
	x = np.log10(np.array([freq1,freq2,nvss_freq]))
	y = np.log10(np.array([f1_array,f2_array, f3_tf]))
	
	# the uncertainty should take the log into account (y = log(x), delta_y = delta_x/x)
	weights = np.array([np.divide(f1_e_array,f1_array), np.divide(f2_e_array,f2_array), np.divide(f3_e_tf,f3_tf)])
#	weights = np.square(weights)	# apparently np.polyfit should use 1/sigma as the weights (no need to square)
	weights = np.reciprocal(weights)
	
	normal_sum = np.zeros((3,1))
	for i in range(3):
		normal_sum[i] = np.sum(weights[i])
	for i in range(3):
		weights[i] = np.divide(weights[i],normal_sum[i])
		
	y = y.T
	weights = weights.T
	
	for i in range(len(f1_array)):
		p, cov = np.polyfit( x, y[i],1,w=weights[i],cov=True)	# note polyfit can be vectorized, but then only takes constant weights (=> not vectorizing it)
		spec_indx[i] = p[0]	# highest order first
		spec_indx_err[i] = cov[0,0]

# Warning: np.polyfit scales the covariance matrix by fac = resids / (len(x) - order - 2.0)	
# see https://mail.scipy.org/pipermail/numpy-discussion/2013-July/067076.html
# For a 3 point fit, we can correct by just multiplying by minus one (worse for a 4 point fit!)
	spec_indx_err = np.sqrt(-spec_indx_err)
	

	df1.loc[f1_index, 'NVSS_SI'] = spec_indx
	df2.loc[f2_index, 'NVSS_SI'] = spec_indx
	
	df1.loc[f1_index, 'NVSS_E_SI'] = spec_indx_err
	df2.loc[f2_index, 'NVSS_E_SI'] = spec_indx_err
	
	print('Sample complete 3 point NVSS S.I. = '+str(spec_indx[0])+' +/ '+str(spec_indx_err[0]))
	
	
	# calculate offsets
	
	df1['NVSS_RA_offset'] = -999.9
	df1['NVSS_DEC_offset'] = -999.9
	df2['NVSS_RA_offset'] = -999.9
	df2['NVSS_DEC_offset'] = -999.9
	
	df1.loc[df1['NVSS'],'NVSS_RA_offset'] = df1.loc[df1['NVSS'],'RA'].values - df1.loc[df1['NVSS'],'NVSS_RA'].values
	df1.loc[df1['NVSS'],'NVSS_DEC_offset'] = df1.loc[df1['NVSS'],'DEC'].values - df1.loc[df1['NVSS'],'NVSS_DEC'].values
	df2.loc[df2['NVSS'],'NVSS_RA_offset'] = df2.loc[df2['NVSS'],'RA'].values - df2.loc[df2['NVSS'],'NVSS_RA'].values
	df2.loc[df2['NVSS'],'NVSS_DEC_offset'] = df2.loc[df2['NVSS'],'DEC'].values - df2.loc[df2['NVSS'],'NVSS_DEC'].values
	
else:
	df1['NVSS'] = False
	df2['NVSS'] = False
	
	
#########################################################################################################
#
#
#	Almost done. Add cat ID and concat.
#	Make a separate concat for ML
#
#
#########################################################################################################


# Generate a cat ID

df1['Cat_ID'] = -999
df2['Cat_ID'] = -999
df1['Cat_ID'], df2['Cat_ID'] = gen_cat_id(df1, df2)

# concat the two data frames and resort the Cat IDs

dfc = (pd.concat([df1,df2], axis=0)).sort_values(by=['RA','Freq'],ascending=[1,1])	# Important to also sort frequency here (latex printout assumes low freq first)
resort_cat_id(dfc)	

# Get the normal RA and DEC coords for all sources. Also generate the name array.

dfc['S_Ra_Dec'] = deg_to_sexagesimal(dfc,'&',2)
dfc['Name'] = np.core.defchararray.replace( np.core.defchararray.add('J', deg_to_sexagesimal(dfc,'',1)) ,':','')
dfc.reset_index(inplace=True, drop=True)
dfc.reindex(index=range(0,len(dfc)))


#########################################################################################################
#
#
#	GBS counterpart check
#
#
#########################################################################################################

if have_gbs_counterparts:
	vcat_df = read_vizier(gbs_counterparts_file, 55,['_RAJ2000','_DEJ2000','Name','Other','Type'])	# read in the GBS ra and decs in degrees
	vcat_df = vcat_df.drop(vcat_df.index[0:2])	# clip units line and ----- line
	vcat_df.reset_index(inplace=True, drop=True)
	vcat_df.reindex(index=range(0,len(vcat_df)))
	
	distances = scipy.spatial.distance.cdist(np.array([dfc['RA'],dfc['DEC']]).T , np.array([vcat_df['_RAJ2000'], vcat_df['_DEJ2000']]).T,'euclidean')
	min_dist = np.amin(distances,axis=1)
	min_dist_index = np.argmin(distances,axis=1)
	dfc['GBS_C'] =  min_dist<gbs_radius
	if np.sum(dfc['GBS_C'].values)>0:
		dfc.loc[dfc['GBS_C'], 'GBS_C_Other'] = vcat_df.loc[min_dist_index[min_dist<gbs_radius], 'Other'].values
		dfc.loc[dfc['GBS_C'], 'GBS_C_Type'] = vcat_df.loc[min_dist_index[min_dist<gbs_radius], 'Type'].values
		dfc.loc[dfc['GBS_C'], 'GBS_C_Name'] = vcat_df.loc[min_dist_index[min_dist<gbs_radius], 'Name'].values
		print('\tDetected '+str(np.sum(dfc['GBS_C'].values))+' Goult belt survey matches with other counterparts at the following sources (across both gmrt freqs):')
		t1 = dfc.loc[dfc['GBS_C'], 'GBS_C_Other'].values
		t2 = dfc.loc[dfc['GBS_C'], 'GBS_C_Type'].values
		t3 = dfc.loc[dfc['GBS_C'], 'GBS_C_Name'].values
		t4 = dfc.loc[dfc['GBS_C'], 'Name'].values
		print('\t\t GMRT name \t\t\t GBS name \t\t\t GBS matches \t\t GBS ID')
		for i in range(len(t1)):
			print('\t\t'+str(t4[i])+'\t\t'+str(t3[i])+'\t '+str(t1[i])+'\t '+str(t2[i]))
	else:
		print('\tNo Goult belt survey matches with other counterparts detected')
		
if have_xmm:
	vcat_df = read_vizier(xmm_file, 0,['_RAJ2000','_DEJ2000'])
	vcat_df = vcat_df.drop(vcat_df.index[0:2])	# clip units line and ----- line
	vcat_df.reset_index(inplace=True, drop=True)
	vcat_df.reindex(index=range(0,len(vcat_df)))
	
	distances = scipy.spatial.distance.cdist(np.array([dfc['RA'],dfc['DEC']]).T , np.array([vcat_df['_RAJ2000'], vcat_df['_DEJ2000']]).T,'euclidean')
	min_dist = np.amin(distances,axis=1)
	min_dist_index = np.argmin(distances,axis=1)
	dfc['XMM'] =  min_dist<gbs_radius
	if np.sum(dfc['XMM'].values)>0:
		print('\tDetected XMM survey matches')
#		t1 = dfc.loc[dfc['XMM'], 'S_Ra_Dec'].values
#		for i in range(len(t1)):
#			print(str(t1[i]))
	else:
		print('\t No XMM survey matches detected...')
	
if have_mass:
	vcat_df = read_vizier(mass_file, 0,['_RAJ2000','_DEJ2000'])
	vcat_df = vcat_df.drop(vcat_df.index[0:2])	# clip units line and ----- line
	vcat_df.reset_index(inplace=True, drop=True)
	vcat_df.reindex(index=range(0,len(vcat_df)))
	
	distances = scipy.spatial.distance.cdist(np.array([dfc['RA'],dfc['DEC']]).T , np.array([vcat_df['_RAJ2000'], vcat_df['_DEJ2000']]).T,'euclidean')
	min_dist = np.amin(distances,axis=1)
	min_dist_index = np.argmin(distances,axis=1)
	dfc['2MASS'] =  min_dist<gbs_radius
	if np.sum(dfc['2MASS'].values)>0:
		print('\tDetected 2MASS survey matches')
#		t1 = dfc.loc[dfc['2MASS'], 'S_Ra_Dec'].values
#		for i in range(len(t1)):
#			print(str(t1[i]))
	else:
		print('\t No 2MASS survey matches detected...')
		
if have_spitzer:
	vcat_df = read_vizier(spitzer_file, 0,['_RAJ2000','_DEJ2000'])
	vcat_df = vcat_df.drop(vcat_df.index[0:2])	# clip units line and ----- line
	vcat_df.reset_index(inplace=True, drop=True)
	vcat_df.reindex(index=range(0,len(vcat_df)))
	
	if(len(vcat_df)>0):
		distances = scipy.spatial.distance.cdist(np.array([dfc['RA'],dfc['DEC']]).T , np.array([vcat_df['_RAJ2000'], vcat_df['_DEJ2000']]).T,'euclidean')
		min_dist = np.amin(distances,axis=1)
		min_dist_index = np.argmin(distances,axis=1)
		dfc['SPITZER'] =  min_dist<gbs_radius
		if np.sum(dfc['SPITZER'].values)>0:
			print('\tDetected Spitzer survey matches')
#			t1 = dfc.loc[dfc['SPITZER'], 'S_Ra_Dec'].values
#			for i in range(len(t1)):
#				print(str(t1[i]))
		else:
			print('\t No Spitzer survey matches detected...')
	else:
			print('\t No Spitzer survey matches detected...')
			have_spitzer = False
			
# print out GBS matches	
			
print('GBS matches:')
print('\t\t GMRT name \t\t\t GBS name \t\t\t Freq')
t2 = dfc.loc[dfc['GBS'], 'Freq'].values
t4 = dfc.loc[dfc['GBS'], 'Name'].values
t3 = dfc.loc[dfc['GBS'], 'GBS_Name'].values
for i in range(len(t1)):
	print('\t\t'+str(t4[i])+'\t\t'+str(t3[i])+'\t '+str(t2[i]))
print('End of GBS matches:')
#########################################################################################################
#
#
#	FR1/2 spotting
#	0 = neither
#	1 = FR1
#	2 = FR2
#
#
#########################################################################################################

dfc['MLC']=0
if have_class:
	# Read in classification list

	class_df = pd.read_csv(class_list,skiprows=0,delimiter=' ', engine='python')

	# Convert to degrees for easy maths

	class_df['RA'] = (class_df['ra_h'].values + class_df['ra_m'].values/60.0 +class_df['ra_s'].values/3600.0)*15.0
	class_df['DEC'] = (class_df['dec_h'].values + class_df['dec_m'].values/60.0 +class_df['dec_s'].values/3600.0)
	class_df['radius_deg'] = class_df['dec_h'].values/3600.0

	# Find sources that have a match in the identifier list

	distances = scipy.spatial.distance.cdist(np.array([dfc['RA'],dfc['DEC']]).T , np.array([class_df['RA'],class_df['DEC']]).T,'euclidean')
	min_dist = np.amin(distances,axis=1)
	min_dist_index = np.argmin(distances,axis=1) # Would fail if the entire array was NaN

	# Apply corresponding classification

	dfc.loc[min_dist < class_df.loc[min_dist_index, 'radius_deg'].values, 'MLC'] = class_df.loc[min_dist_index, 'type'].values


	
##############################################
#
#	Now just need to output data in a useful formats
#
#		- csv file with all columns
#		- latex table (suitable for paper)
#		- kvis annotation file
#		- list of coords for NVSS lookup (already done for each frequency)
#
#
##############################################

	
# write out full csv file, a kvis annotation file and a list of coords for NVSS lookup

dfc.to_csv(outputname+'.csv')
make_kvis_file(outputname, dfc)

# Write out fluxes for plot

df1.to_csv(outputname+'_freq1.fluxes.csv',columns=['Total_flux','E_Total_flux', 'Peak_flux', 'E_Peak_flux'],index=False)
df2.to_csv(outputname+'_freq2.fluxes.csv',columns=['Total_flux','E_Total_flux', 'Peak_flux', 'E_Peak_flux'],index=False)

# Write out island RMS as a measure of noise at each source

df1.to_csv(outputname+'_freq1.island_residuals.csv',columns=['Resid_Isl_rms'],index=False)
df2.to_csv(outputname+'_freq2.island_residuals.csv',columns=['Resid_Isl_rms'],index=False)

# Write out coords and NVSS separations for plots

dfc[['RA','DEC']].to_csv(outputname+'_positions.csv')

if(have_nvss):
	dfc.loc[np.all([dfc['NVSS'], dfc['Freq']==freq1, dfc['S_Code'] == 'S', dfc['Peak_flux'] > 20.0*dfc['Resid_Isl_rms']],axis=0),['NVSS_RA_offset','NVSS_DEC_offset']].to_csv(outputname+'_freq1.nvss_offset.csv')
	dfc.loc[np.all([dfc['NVSS'], dfc['Freq']==freq2, dfc['S_Code'] == 'S', dfc['Peak_flux'] > 20.0*dfc['Resid_Isl_rms']],axis=0),['NVSS_RA_offset','NVSS_DEC_offset']].to_csv(outputname+'_freq2.nvss_offset.csv')

# Write out spectral indicies

dfc.loc[dfc['Match']].to_csv(outputname+'.spx.csv',columns=['SI','E_SI'],index=False)
if(have_nvss):
	dfc.loc[dfc['NVSS']].to_csv(outputname+'.nvss_spx.csv',columns=['NVSS_SI','NVSS_E_SI','NVSS_SI_FP'],index=False)


##############################################
#
#	Changes made for Latex table - dfc is not suitable for anything else after this point!
#
#
##############################################

# Do some rounding to 2 decimal places. Always round up errors by adding 0.005

dfc['Total_flux'] = np.around(dfc['Total_flux'],2)
dfc['E_Total_flux'] = np.around(np.add(dfc['E_Total_flux'],0.005),2)
dfc['Peak_flux'] = np.around(dfc['Peak_flux'],2)
dfc['E_Peak_flux'] = np.around(np.add(dfc['E_Peak_flux'],0.005),2)
dfc['Maj'] = np.around(dfc['Maj'],2)
dfc['E_Maj'] = np.around(np.add(dfc['E_Maj'],0.005),2)
dfc['Min'] = np.around(dfc['Min'],2)
dfc['E_Min'] = np.around(np.add(dfc['E_Min'],0.005),2)
dfc['PA'] = np.around(dfc['PA'],2)
dfc['E_PA'] = np.around(np.add(dfc['E_PA'],0.005),2)
dfc['DC_Maj'] = np.around(dfc['Maj'],2)
dfc['DC_E_Maj'] = np.around(np.add(dfc['E_Maj'],0.005),2)
dfc['DC_Min'] = np.around(dfc['Min'],2)
dfc['DC_E_Min'] = np.around(np.add(dfc['E_Min'],0.005),2)
dfc['DC_PA'] = np.around(dfc['PA'],2)
dfc['DC_E_PA'] = np.around(np.add(dfc['E_PA'],0.005),2)
dfc['SI'] = np.around(dfc['SI'],2)
dfc['E_SI'] = np.around(np.add(dfc['E_SI'],0.005),2)
dfc['Resid_Isl_rms'] = dfc['Resid_Isl_rms']*1000000.0
dfc['Resid_Isl_rms'] = np.add(dfc['Resid_Isl_rms'],0.5).astype(int)
if(have_nvss):
	dfc['NVSS_SI'] = np.around(dfc['NVSS_SI'],2)
	dfc['NVSS_E_SI'] = np.around(np.add(dfc['NVSS_E_SI'],0.005),2)

# Write out file for as latex table

# For sources with the same cat ID report the 610 location only

# prepare the column headers and units

header_start = 'Name&Ra&Dec'
header_325 = '$S_{\\rm peak,325\,MHz}$&$S_{\\rm int,325\,MHz}$&$\\rm Resid\_Isl\_rms_{\\rm 325\,MHz}$&$\\rm S\_Code_{\\rm 325\,MHz}$'
header_610 = '$S_{\\rm peak,610\,MHz}$&$S_{\\rm int,610\,MHz}$&$\\rm Resid\_Isl\_rms_{\\rm 610\,MHz}$&$\\rm S\_Code_{\\rm 610\,MHz}$'
header_gen = '$S_{\\rm peak}$&$S_{\\rm int}$&$\\rm Resid\_Isl\_rms$&$\\rm S\_Code$'
header_end = '$\\alpha_{\\rm GMRT}$&NVSS&$\\alpha_{\\rm NVSS}$&Known'

header2 = '&&'+'&'+header_gen+'&'+header_gen+'&'+'&&&counterparts'+'\\\\\n'

units_start = 'GMRT-TAU&$(h:m:s)$&$(d:m:s)$'
units_per_freq = '(mJy beam$^{-1}$)&(mJy)&($\\mu$Jy beam$^{-1}$)&'
units_end = '&&&'

units_line = units_start+'&'+units_per_freq+'&'+units_per_freq+'&'+units_end+'\\\\\n'
cat_header = header_start+'&'+'323 MHz&&&'+'&'+'608 MHz&&&'+'&'+header_end+'\\\\\n'
cat_header_format='c'
for i in range(cat_header.count('&')):
	cat_header_format = cat_header_format + 'c'

n = len(dfc)
done = np.zeros((n,1),dtype=np.bool)
tx = 0
with open(outputname+'.table.tex', 'w') as f:
	i = 0	# points at current source, or low freq in case of match
	f.write('\\begin{tabular}{'+cat_header_format+'}\n')
	f.write('\\hline\n')
	f.write(cat_header)
	f.write(header2)
	f.write(units_line)
	f.write('\\hline\n')
	while( i < n):
		if(~done[i]):
			done[i] = True
			match = False
			cat_match = False
			cat_list=[]
			j = i	# points at current source, or low freq in case of match
			k = i	# points at current source, or high freq in case of match
#			f.write(str(dfc['Cat_ID'].iloc[i]+1)+'&')	# write out cat ID
			
			if(dfc['Match'].iloc[i]):
				match = False
				if(i +1 != n):
					for k in range(i+1,n):
						if(dfc['Cat_ID'].iloc[j] == dfc['Cat_ID'].iloc[k]):	# if there is a match, report the 610 data as k
#							print(str(dfc['S_Ra_Dec'].iloc[j])+' --- matched with '+str(dfc['S_Ra_Dec'].iloc[k]))
							done[k] = True
							tx = tx+1
							match = True
							if(dfc['Freq'].iloc[j] > dfc['Freq'].iloc[k]):
								temp = j
								j = k
								k = temp
							break
						else:
							match = False
			# generate name
			f.write(str(dfc['Name'].iloc[k]+'&'))
			f.write(str(dfc['S_Ra_Dec'].iloc[k])+'&')	# Write out Ra+Dec (for high freq if both are available for source)

			# Write out 325 data if 325 detection. Do not include Gaussian size for unresolved or M class sources.
			if(dfc['Freq'].iloc[j] == freq1):
				print_entry(dfc.iloc[[j]],True)
			else:
				print_entry(dfc.iloc[[j]],False)
			
			
			# Write out 625 fluxes if a combined detection, otherwise 625 if a 625 only detection
		
			if(match):
				if(dfc['Freq'].iloc[k] == freq2):
					print_entry(dfc.iloc[[k]],True)
				else:
					print_entry(dfc.iloc[[k]],False)
			else:
				if(dfc['Freq'].iloc[j] == freq2):
					print_entry(dfc.iloc[[j]],True)
				else:
					print_entry(dfc.iloc[[j]],False)
				
			# Spectral index time( if matched)
		
			if(match):
#				if((dfc['S_Code'].iloc[j]!='M') &(dfc['S_Code'].iloc[k]!='M')):
				if((dfc['S_Code'].iloc[j]=='S') &(dfc['S_Code'].iloc[k]=='S')):
					f.write('$%.2f'%dfc['SI'].iloc[j]+'\pm%.2f'%dfc['E_SI'].iloc[j]+'$&')
				else:
					f.write('-&')
			else:
				f.write('-&')
				
			if(dfc['NVSS'].iloc[k]&(dfc['S_Code'].iloc[k]=='S')):
				f.write(str(dfc['NVSS_SI_FP'].iloc[k])+'&$%.2f'%dfc['NVSS_SI'].iloc[k]+'\pm%.2f'%dfc['NVSS_E_SI'].iloc[k]+'$&')
				cat_match = True
				cat_list.append('N ')
			else:
				if(dfc['NVSS'].iloc[j]&(dfc['S_Code'].iloc[j]=='S')):
					f.write(str(dfc['NVSS_SI_FP'].iloc[j])+'&$%.2f'%dfc['NVSS_SI'].iloc[j]+'\pm%.2f'%dfc['NVSS_E_SI'].iloc[j]+'$&')
					cat_match = True
					cat_list.append('N')
				else:
					f.write('-&-&')
			
			if have_gbs:		
				if(dfc['GBS'].iloc[j] or dfc['GBS'].iloc[k]):
					cat_match = True
					cat_list.append('G')
				
				'''
				if(i==k):	# no match, check freq and enter label
					if(dfc['Freq'].iloc[j] == freq1):
						f.write('L&'+str(dfc['GBS_Flux'].iloc[j]))
					else:
						f.write('H&'+str(dfc['GBS_Flux'].iloc[j]))
				else:	# have match, at least one of which has a GBS match
					if(dfc['GBS'].iloc[j] and dfc['GBS'].iloc[k]):
						f.write('B&'+str(dfc['GBS_Flux'].iloc[j]))
					else:
						if(dfc['GBS'].iloc[j]):
							f.write('L&'+str(dfc['GBS_Flux'].iloc[j]))
						else:
							f.write('H&'+str(dfc['GBS_Flux'].iloc[j]))
				'''
			if(have_xmm):
				if(dfc['XMM'].iloc[j] or dfc['XMM'].iloc[k]):
					cat_match = True
					cat_list.append('X')
			if(have_mass):
				if(dfc['2MASS'].iloc[j] or dfc['2MASS'].iloc[k]):
					cat_match = True
					cat_list.append('2M')
			if(have_spitzer):
				if(dfc['SPITZER'].iloc[j] or dfc['SPITZER'].iloc[k]):
					cat_match = True
					cat_list.append('S')
				
			if(cat_match):
				f.write(' '.join(cat_list))
			else:
				f.write('-')

			f.write('\\\\\n')

		i = i + 1

	f.write('\\hline\n')
	f.write('\\end{tabular}')

if have_gbs:
	tn = (dfc.loc[dfc['GBS'], 'S_Ra_Dec'].values)
	t1 = (dfc.loc[dfc['GBS'], 'Total_flux'].values)
	t2 = (dfc.loc[dfc['GBS'], 'GBS_Flux'].values).astype(np.float)
	t1 = t1 < t2

	print('Possible thermal emission compared to GBS detected from the following sources:')
	for i in range(len(t1)):
		if t1[i]:
			print('\t'+tn[i].replace('&',' '))

print('Catalog complete.')
print(str(nmatches)+' matches detected.')
print(str(len(df1))+' sources in file 1, freq = '+str(freq1/1.0E6)+' MHz')

print(str(len(df2))+' sources in file 2, freq = '+str(freq2/1.0E6)+' MHz')
print('This implies '+str(len(df1) + len(df2) - nmatches)+' unique detections.')
