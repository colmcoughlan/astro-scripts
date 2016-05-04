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
data_dir = str(os.environ['data_dest'])+'/'
skymodel = str(os.environ['skymodel'])

casa_exec=str(os.environ['casa_exec'])

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


ending='_uv.dppp.MS'

# BBS based direction independent phase calibration using skymodel

ps=[]
for s in range(nsb):
	csb = sb + s
	for i in scanlist:
		cmdstr = 'calibrate-stand-alone -t 1 -f '+wd+data_dir+'L'+str(i)+'_SB%03d'%csb+ending+' '+sd+'phase_cal.bbs '+gsmskymodel+' > '+wd+log_dir+'L'+str(i)+'_SB%03d_phase_cal.log'%csb
		print(cmdstr)
		p = Popen(cmdstr, shell=True)
		ps.append(p)
		if(len(ps)==maxthread):
			waitfunc(ps)
			ps=[]
waitfunc(ps)
