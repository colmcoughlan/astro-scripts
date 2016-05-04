#	LOFAR calibration pipeline


import os
from subprocess import Popen
from time import sleep
from numpy import ceil

# Read in environment variables from script

lofar_home=str(os.environ['LOFARROOT'])+'/'
wd=str(os.environ['working_dir'])+'/'
sd=str(os.environ['script_dir'])+'/'
sb=int(os.environ['subband'])
nsb=int(os.environ['num_subbands'])
lstart=int(os.environ['l_start'])
lend=int(os.environ['l_end'])
lampcal=int(os.environ['l_ampcal'])
od = str(os.environ['data_origin'])+'/'
data_dir = str(os.environ['data_dest'])+'/'
skymodel = str(os.environ['skymodel'])

parmpath=str(os.environ['parmpath'])

maxthread=int(os.environ['nthreads'])

scanlist = range(lstart , lend+1, 3)

# Function that will wait until all processes have finished

def waitfunc(ps):
	while True:
  		ps_status = [p.poll() for p in ps]
		if all([x is not None for x in ps_status]):
  			break
		sleep(1)	# Only check every second
	return(0)

# set the default ending for LOFAR files

ending='_uv.dppp.MS'
#gsmskymodel='/mnt/home_cr/coughlan/lofar/skymodels/GSMSkymodel_TTAU.bbs'

# Create some folders to keep things tidy
# If the calibration folder already exists - no need to do calibration

log_dir = 'logs/'
ps=[]
if not os.path.exists(log_dir):
	cmdstr = 'mkdir '+log_dir
	p = Popen(cmdstr, shell=True)
	ps.append(p)
	if(len(ps)==1):
		waitfunc(ps)
		ps=[]

table_dir = 'calibration_tables/'

if not os.path.exists(table_dir):
	do_cal = True
	cmdstr = 'mkdir '+table_dir
	p = Popen(cmdstr, shell=True)
	ps.append(p)
	if(len(ps)==1):
		waitfunc(ps)
		ps=[]
else:
	print(table_dir+' folder detected. No fresh amplitude calibration will be performed.')
	do_cal = False


if not os.path.exists(od):
	print(od+' folder not detected. Assuming current directory is a suitable working directory and contains the data.')
	od=''
else:
	print(od+' folder detected. Assuming data is stored within and it is a suitable working directory.')

# Remove Internation baselines from file

print 'Removing long baselines'
for s in range(nsb):
	csb = sb + s
	for i in scanlist:
		cmdstr = 'NDPPP '+sd+'remove_long_baselines.ndppp msin='+wd+od+'L'+str(i)+'_SB%03d'%csb+ending+' msout='+wd+data_dir+'L'+str(i)+'_SB%03d'%csb+ending+' > '+wd+log_dir+'L'+str(i)+'_SB%03d'%csb+'_remove_internation_stations.ndppp.log'
		print(cmdstr)
		p = Popen(cmdstr, shell=True)
		ps.append(p)
		if(len(ps)==maxthread):
			waitfunc(ps)
			ps=[]
	cmdstr = 'NDPPP '+sd+'remove_long_baselines.ndppp msin='+wd+od+'L'+str(lampcal)+'_SB%03d'%csb+ending+' msout='+wd+data_dir+'L'+str(lampcal)+'_SB%03d'%csb+ending+' > '+wd+log_dir+'L'+str(lampcal)+'_SB%03d'%csb+'_remove_internation_stations.ndppp.log'
	print(cmdstr)
	p = Popen(cmdstr, shell=True)
	ps.append(p)
	if(len(ps)==maxthread):
		waitfunc(ps)
		ps=[]
waitfunc(ps)
ps=[]



# Amplitude calibration (maxthread/nsb threads each)

ps=[]

if do_cal:

	for i in range(nsb):
		csb = sb + i
		fn = wd+data_dir+'L'+str(lampcal)+'_SB%03d'%csb+ending
		cmdstr = 'calibrate-stand-alone -f '+fn+' '+sd+'amp_cal.bbs ' +lofar_home+'share/pipeline/skymodels/3C147.skymodel > '+wd+log_dir+'SB%03d_amp_cal.log'%csb
		print('Calibrating gain amplitudes on subband '+str(csb))
		print(cmdstr)
		p = Popen(cmdstr,shell=True)
		ps.append(p)
		if(len(ps)==maxthread):
			waitfunc(ps)
			ps=[]
	waitfunc(ps)
	ps=[]

	print('Now exporting gain amplitudes:')

	for i in range(nsb):
		csb = sb + i
		fn = wd+data_dir+'L'+str(lampcal)+'_SB%03d'%csb+ending
		cmdstr=parmpath+' in='+fn+'/instrument out='+wd+table_dir+'L'+str(lampcal)+'_SB%03d_gain_solutions.table'%csb
		print(cmdstr)
		p = Popen(cmdstr, shell=True)
		ps.append(p)
		if(len(ps)==maxthread):
			waitfunc(ps)
			ps=[]
	waitfunc(ps)
	ps=[]
else:
	print('Calibration tables detected - no need to recreate them.')


# Apply amplitude calibration
# Use 1 thread per each of the 8 scans = 8 threads
print('Applying amplitude calibration to all scans/subbands')

ps=[]
for s in range(nsb):
	csb = sb + s
	for i in scanlist:
		cmdstr='calibrate-stand-alone -t 1 --parmdb '+wd+table_dir+'L'+str(lampcal)+'_SB%03d_gain_solutions.table '%csb+wd+data_dir+'L'+str(i)+'_SB%03d'%csb+ending+' '+sd+'correct_only.bbs > '+wd+log_dir+'L'+str(i)+'_SB%03d'%csb+'_apply_gain_solutions.log'
		print(cmdstr)
		p = Popen(cmdstr, shell=True)
		ps.append(p)
		if(len(ps)==maxthread):
			waitfunc(ps)
			ps=[]

waitfunc(ps)
ps=[]


# BBS based direction independent phase calibration using provided skymodel

ps=[]
for s in range(nsb):
	csb = sb + s
	for i in scanlist:
		cmdstr = 'calibrate-stand-alone -t 1 -f '+wd+data_dir+'L'+str(i)+'_SB%03d'%csb+ending+' '+sd+'phase_cal.bbs '+skymodel+' > '+wd+log_dir+'L'+str(i)+'_SB%03d_phase_cal.log'%csb
		print(cmdstr)
		p = Popen(cmdstr, shell=True)
		ps.append(p)
		if(len(ps)==maxthread):
			waitfunc(ps)
			ps=[]
waitfunc(ps)




print('Calibration run complete')


