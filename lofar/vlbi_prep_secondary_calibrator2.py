#	LOFAR calibration pipeline


import os
from subprocess import Popen
from time import sleep
from numpy import ceil

# Read in environment variables from script

wd=str(os.environ['working_dir'])+'/'
sd=str(os.environ['script_dir'])+'/'
sb=int(os.environ['subband'])
nsb=int(os.environ['num_subbands'])
lstart=int(os.environ['l_start'])
lend=int(os.environ['l_end'])

casa_exec=str(os.environ['casa_exec'])

maxthread=int(os.environ['nthreads'])

shift_coords=str(os.environ['shift_coords'])

scanlist = range(lstart , lend+1, 3)

mscorpol_exec='/mnt/home_cr/coughlan/lofar/scripts/mscorpol/mscorpol.py'

sbs = range(sb,sb+nsb,1)
print(str(len(sbs)))
bad_sb=[185]	# remove bad sb for dgtau data
sb_per_if = 16
nif = int(ceil(float(nsb)/float(sb_per_if)))


# Function that will wait until all processes have finished

def waitfunc(ps):
	while True:
  		ps_status = [p.poll() for p in ps]
		if all([x is not None for x in ps_status]):
  			break
		sleep(1)	# Only check every second
	return(0)

# Function to write casa parset to concat and export files for AIPS

def casa_concat_parset(vis, concatvis):
	try:
		f = open(concatvis+'.py','w')
	except:
		print "\t Error opening ", concatvis
		return(-1)


	f.write('concat(vis='+str(vis)+', concatvis = \''+concatvis+'.ms\', async = false)\n')

	f.close
	return(0)


############################################################################################################
#
#
#
#	Start here
#
#
######################################################################################################

print('Shifting phase centre to '+shift_coords+' and preparing AIPS UVDATA file.')


ending='_uv.dppp.MS'

log_dir = 'ppc_logs/'
ps=[]
if not os.path.exists(log_dir):
	cmdstr = 'mkdir '+log_dir
	p = Popen(cmdstr, shell=True)
	ps.append(p)
	if(len(ps)==1):
		waitfunc(ps)
		ps=[]



shifted_data_parent_dir = 'ppc_shifted_data/'
shifted_data_dir = shifted_data_parent_dir+'shifted/'
shifted_phased_data_dir = shifted_data_parent_dir+'shifted_phased/'
ps=[]
if not os.path.exists(shifted_data_parent_dir):
	cmdstr = 'mkdir '+shifted_data_parent_dir
	p = Popen(cmdstr, shell=True)
	ps.append(p)
	if(len(ps)==1):
		waitfunc(ps)
		ps=[]
	cmdstr = 'mkdir '+shifted_data_dir
	p = Popen(cmdstr, shell=True)
	ps.append(p)
	if(len(ps)==1):
		waitfunc(ps)
		ps=[]
	cmdstr = 'mkdir '+shifted_phased_data_dir
	p = Popen(cmdstr, shell=True)
	ps.append(p)
	if(len(ps)==1):
		waitfunc(ps)
		ps=[]

od = 'data/'

if not os.path.exists(od):
	print(od+' folder not detected. Assuming current directory is a suitable working directory and contains the data.')
	odir=''
else:
	print(od+' folder detected. Assuming data is stored within and it is a suitable working directory.')



# Shift phase centre to desired coords
ps=[]
print 'Shifting phase centre to '+shift_coords+' on all scans/sub-bands.'
for csb in sbs:
	if(not(csb in bad_sb)):
		for i in scanlist:
			vis = wd + od +  'L'+str(i)+'_SB%03d'%csb+ending
			visout = wd +  shifted_data_dir + 'L'+str(i)+'_SB%03d'%csb+ending
			cmdstr='NDPPP '+sd+'vlbi_shift_phase_centre.dppp phase_shift.phasecenter=\''+shift_coords+'\' msin='+ vis +' msin.datacolumn=CORRECTED_DATA msout='+visout+' > '+wd+log_dir+'L'+str(i)+'_SB%03d'%csb+'_shift_phase_centre.ndppp.log'
			print(cmdstr)
			p = Popen(cmdstr, shell=True)
			ps.append(p)
			if(len(ps)==maxthread):
				waitfunc(ps)
				ps=[]
waitfunc(ps)
ps=[]


print('Making superstation and removing CS from result in single step...')
ps=[]
ifnum=1
for csb in sbs:
	if(csb-sbs[0]+1 > ifnum*sb_per_if):
		ifnum = ifnum+1
	if(not(csb in bad_sb)):
		for i in scanlist:
			cmdstr = 'NDPPP '+sd+'vlbi_make_superstation_multi.ndppp msin='+wd+shifted_data_dir+'L'+str(i)+'_SB%03d'%csb+ending+' msin.datacolumn=DATA msout='+wd+shifted_phased_data_dir+'L'+str(i)+'_SB%03d_IF'%csb+str(ifnum)+ending+' > '+wd+log_dir+'L'+str(i)+'_SB%03d_IF'%csb+str(ifnum)+'_make_superstation_multi.ndppp.log'
			print(cmdstr)
			p = Popen(cmdstr, shell=True)
			ps.append(p)
			if(len(ps)==maxthread):
				waitfunc(ps)
				ps=[]

waitfunc(ps)
ps=[]


# Merge the data by sub-band and then concat by time. Note some combination of phasing up and merging causes bad UVW points that need to be flaggged

print 'Concatentating data, first by SB, then by time.'

ps=[]
for aips_if in range(nif):
	for i in scanlist:
		ms_list=['[']
		for csb in range(sb+aips_if * sb_per_if,sb + (aips_if+1) * sb_per_if):
			ms_list.append(wd+shifted_phased_data_dir+'L'+str(i)+'_SB%03d_IF'%csb+str(aips_if+1)+ending)
			if csb != sb + (aips_if+1) * sb_per_if - 1:
				ms_list.append(',')
		ms_list.append(']')
		cmdstr = 'NDPPP '+sd+'vlbi_merge.ndppp msin='+''.join(ms_list)+' msout='+wd+shifted_data_parent_dir+'L'+str(i)+'_IF'+str(aips_if+1)+'.ms'+' > '+wd+log_dir+'L'+str(i)+'_IF'+str(aips_if+1)+'.ndppp.log'
		print(cmdstr)
		p = Popen(cmdstr, shell=True)
		ps.append(p)
		if(len(ps)==maxthread):
			waitfunc(ps)
			ps=[]
waitfunc(ps)
ps=[]

for aips_if in range(nif):
	ps=[]
	concat_list=[]
	for i in scanlist:
		concat_list.append(wd+shifted_data_parent_dir+'L'+str(i)+'_IF'+str(aips_if+1)+'.ms')
	fn='ppc_IF'+str(aips_if+1)
	casa_concat_parset(concat_list, fn)
	cmdstr=casa_exec+' -c '+fn+'.py'
	print(cmdstr)
	os.system(cmdstr)

print('\n\nFlagging bad new UVW points(if necessary)\n\n')
ps=[]
for aips_if in range(nif):
	fn='ppc_IF'+str(aips_if+1)+'.ms'
	cmdstr = 'taql'+' \'update '+fn+' set FLAG=T where isnan(UVW[0])\' > '+wd+log_dir+'ppc_IF'+str(aips_if+1)+'_flag_bad_uv_coords.taql.log'
	print(cmdstr)
	p = Popen(cmdstr, shell=True)
	ps.append(p)
	if(len(ps)==maxthread):
		waitfunc(ps)
		ps=[]

waitfunc(ps)
ps=[]
	


# Now convert from linear to circular polarisation and export

print '\n\nConverting from linear to circular polarisation and applying beam correction. Exporting result as FITS file.\n\n'

ps=[]
for aips_if in range(nif):
	fn='ppc_IF'+str(aips_if+1)+'.ms'
	cmdstr = mscorpol_exec+' -f '+fn+' > '+wd+log_dir+'ppc_IF'+str(aips_if+1)+'.mscorpol.log'
	print(cmdstr)
	p = Popen(cmdstr, shell=True)
	ps.append(p)
	if(len(ps)==maxthread):
		waitfunc(ps)
		ps=[]

waitfunc(ps)
ps=[]

for aips_if in range(nif):
	fn='ppc_IF'+str(aips_if+1)
	cmdstr='ms2uvfits ms='+fn+'.ms fitsfile='+fn+'.fits column=DATA writesyscal=FALSE > '+wd+log_dir+'ppc_IF'+str(aips_if+1)+'.ms2uvfits.log'
	print(cmdstr)
	os.system(cmdstr)



print(str(nif)+' IFs calibrated and prepared for fringe fitting in AIPS.')

