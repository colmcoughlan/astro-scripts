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
		cmdstr = 'NDPPP '+sd+'remove_long_baselines_noUK.ndppp msin='+wd+od+'L'+str(i)+'_SB%03d'%csb+ending+'  msout='+wd+data_dir+'L'+str(i)+'_SB%03d'%csb+ending+' > '+wd+log_dir+'L'+str(i)+'_SB%03d'%csb+'_remove_internation_stations.ndppp.log'
		print(cmdstr)
		p = Popen(cmdstr, shell=True)
		ps.append(p)
		if(len(ps)==maxthread):
			waitfunc(ps)
			ps=[]
	cmdstr = 'NDPPP '+sd+'remove_long_baselines_noUK.ndppp msin='+wd+od+'L'+str(lampcal)+'_SB%03d'%csb+ending+' msout='+wd+data_dir+'L'+str(lampcal)+'_SB%03d'%csb+ending+' > '+wd+log_dir+'L'+str(lampcal)+'_SB%03d'%csb+'_remove_internation_stations.ndppp.log'
	print(cmdstr)
	p = Popen(cmdstr, shell=True)
	ps.append(p)
	if(len(ps)==maxthread):
		waitfunc(ps)
		ps=[]
waitfunc(ps)
ps=[]



# Flag data

print 'Flagging data with AOFlagger'
for s in range(nsb):
	csb = sb + s
	for i in scanlist:
		fn = wd+data_dir+'L'+str(i)+'_SB%03d'%csb+ending
		cmdstr='NDPPP '+sd+'aoflagger.ndppp msin='+fn+' msin.datacolumn=DATA msout='+fn+' > '+wd+log_dir+'L'+str(i)+'_SB%03d.aoflagger_log.ndppp'%csb
		print(cmdstr)
		os.system(cmdstr)
	fn = wd+data_dir+'L'+str(lampcal)+'_SB%03d'%csb+ending
	cmdstr = 'NDPPP '+sd+'aoflagger.ndppp msin='+fn+' msin.datacolumn=DATA msout='+fn+' > '+wd+log_dir+'L'+str(lampcal)+'_SB%03d.aoflagger_log.ndppp'%csb
	print(cmdstr)
	os.system(cmdstr)





print('Flagging and removal of Long baselines complete.')


