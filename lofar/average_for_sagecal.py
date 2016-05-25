#	LOFAR calibration pipeline


from os import environ
from subprocess import Popen
from time import sleep
from numpy import ceil

# Read in environment variables from script

lofar_home=str(environ['LOFARROOT'])+'/'
od=str(environ['origin_dir'])+'/'
wd=str(environ['working_dir'])+'/'
sd=str(environ['script_dir'])+'/'
log_dir=str(environ['log_dir'])+'/'
sb=int(environ['subband'])
nsb=int(environ['num_subbands'])
lstart=int(environ['l_start'])
lend=int(environ['l_end'])

scanint = 8


ending='-10_uv.dppp.pre-cal.ms'

maxthread=int(environ['nthreads'])

# Function that will wait until all processes have finished

def waitfunc(ps):
	while True:
  		ps_status = [p.poll() for p in ps]
		if all([x is not None for x in ps_status]):
  			break
		sleep(1)	# Only check every second
	return(0)



# Sagecal-based phase and amplitude calibration. Average to 10s, 2 channels per subband for sagecal.
print('Averaging to 10s and 1 channel for distributed sagecal.')
ps=[]
for s in range(nsb):
	csb = sb + s
	for i in range(lstart,lend+1,scanint):
		fn = od+'L'+str(i)+'_SBgr%03d'%csb+ending
		cmdstr = 'NDPPP '+sd+'average.dppp averager.timestep=2 averager.freqstep=1 msin='+fn+' msin.datacolumn=CORRECTED_DATA msout='+wd+'L'+str(i)+'_SBgr%03d'%csb+ending+' msout.datacolumn=DATA > '+log_dir+'L'+str(i)+'_SB%03d'%csb+'_avg.dppp.log'
		print(cmdstr)
		p = Popen(cmdstr, shell=True)
		ps.append(p)
		if(len(ps) == maxthread):
			waitfunc(ps)
			ps=[]
waitfunc(ps)
ps=[]
print('Adding imaging columns.')
for s in range(nsb):
	csb = sb + s
	for i in range(lstart,lend+1,scanint):
		cmdstr = 'addImagingColumns.py '+wd+'L'+str(i)+'_SBgr%03d'%csb+ending+' > '+log_dir+'L'+str(i)+'_SB%03d'%csb+'_addImagingColumn.log'
		print(cmdstr)
		p = Popen(cmdstr, shell=True)
		ps.append(p)
		if(len(ps) == maxthread):
			waitfunc(ps)
			ps=[]
waitfunc(ps)

print('Process Complete.')
