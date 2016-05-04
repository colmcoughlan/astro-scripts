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

# Function that will wait until all processes have finished

def waitfunc(ps):
	while True:
  		ps_status = [p.poll() for p in ps]
		if all([x is not None for x in ps_status]):
  			break
		sleep(1)	# Only check every second
	return(0)

# Function to write casa parset to concat and export files for AIPS

def casa_parset(vis, concatvis):
	try:
		f = open(concatvis+'.py','w')
	except:
		print "\t Error opening ", outputname
		return(-1)


	f.write('concat(vis='+str(vis)+', concatvis = \''+concatvis+'.ms\', async = false)\n')
	f.write('exportuvfits(vis = \''+str(concatvis)+'.ms\', fitsfile = \''+concatvis+'.fits\', async = false)\n')

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

log_dir = 'logs/'
ps=[]
if not os.path.exists(log_dir):
	cmdstr = 'mkdir '+log_dir
	p = Popen(cmdstr, shell=True)
	ps.append(p)
	if(len(ps)==1):
		waitfunc(ps)
		ps=[]


ps=[]
data_dir = 'data_final/'

if not os.path.exists(data_dir):
	print('Error! '+data_dir+' not detected. This script should be run after vlbi_amp_cal.py.')
	exit()

shifted_data_parent_dir = 'shifted_data3/'
shifted_data_dir = shifted_data_parent_dir+'unmerged/'
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



# Shift phase centre to desired coords

for s in range(nsb):
	csb = sb + s
	for i in scanlist:
		vis = wd + data_dir+'L'+str(i)+'_SB%03d'%csb+ending
		visout = wd + shifted_data_dir+'L'+str(i)+'_SB%03d'%csb+ending
		cmdstr='NDPPP '+sd+'vlbi_shift_phase_centre.dppp phase_shift.phasecenter=\''+shift_coords+'\' msin='+ vis +' msout='+visout+' > '+wd+log_dir+'L'+str(i)+'_SB%03d'%csb+'_shift_phase_centre.ndppp.log'
		print(cmdstr)
		p = Popen(cmdstr, shell=True)
		ps.append(p)
		if(len(ps)==maxthread):
			waitfunc(ps)
			ps=[]



ps=[]
cmdstr = 'mergeSB.py --obsDir='+wd+shifted_data_dir+' --outputDir='+wd+shifted_data_parent_dir+' --column2Merge=DATA'
print(cmdstr)
p = Popen(cmdstr, shell=True)
ps.append(p)
if(len(ps)==1):
	waitfunc(ps)
	ps=[]

# Now convert to linear pol

print 'Converting from linear to circular polarisation'
print 'Updating polarisation table'

ps=[]
for i in scanlist:
	fn=shifted_data_parent_dir+'MergedDATA/L'+str(i)+'_SB%03d'%sb+'_%03d'%(sb+nsb-1)
	cmdstr = 'python '+sd+'lin2circ.py -i '+fn+' -p True -c DATA -o CORRECTED_DATA'
	print(cmdstr)
	p = Popen(cmdstr, shell=True)
	ps.append(p)
	if(len(ps)==maxthread):
		waitfunc(ps)
		ps=[]

# Finally concat and export for AIPS

print 'Concatenating files and converting to FITSUV format.'
ps=[]
concat_list=[]
for i in scanlist:
	concat_list.append(shifted_data_parent_dir+'MergedDATA/L'+str(i)+'_SB%03d'%sb+'_%03d'%(sb+nsb-1))
fn='concat_data_shifted3'
casa_parset(concat_list, fn)
cmdstr=casa_exec+' -c '+fn+'.py'
print(cmdstr)
p = Popen(cmdstr, shell=True)
ps.append(p)
if(len(ps)==1):
	waitfunc(ps)
	ps=[]


print('Calibration run complete')
print('Final circular polarisation data in CORRECTED_DATA column')



