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
sb=int(environ['subband'])
nsb=int(environ['num_subbands'])
lstart=int(environ['l_start'])
lend=int(environ['l_end'])

skymodel=str(environ['skymodel'])
cluster_file=str(environ['cluster_file'])

ending='_uv.dppp.MS'

maxthread=int(environ['nthreads'])

# Function that will wait until all processes have finished

def waitfunc(ps):
	while True:
  		ps_status = [p.poll() for p in ps]
		if all([x is not None for x in ps_status]):
  			break
		sleep(1)	# Only check every second
	return(0)



# Sagecal-based phase and amplitude calibration. Average to 10s for sagecal.
print('Doing robust sagecal calibration. Copying old corrected data to the DATA column.')
ps=[]
print('Averaging to 10s.')
for s in range(nsb):
	csb = sb + s
	for i in range(lstart,lend+1,3):
		fn = od+'L'+str(i)+'_SB%03d'%csb+ending
		cmdstr = 'NDPPP '+sd+'average.dppp averager.timestep=2 averager.freqstep=1 msin='+fn+' msin.datacolumn=CORRECTED_DATA msout='+'L'+str(i)+'_SB%03d'%csb+ending+' msout.datacolumn=DATA > '+wd+'L'+str(i)+'_SB%03d'%csb+'_avg.dppp.log'
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
	for i in range(lstart,lend+1,3):
		cmdstr = 'addImagingColumns.py '+'L'+str(i)+'_SB%03d'%csb+ending+' > '+wd+'L'+str(i)+'_SB%03d'%csb+'_addImagingColumn.log'
		print(cmdstr)
		p = Popen(cmdstr, shell=True)
		ps.append(p)
		if(len(ps) == maxthread):
			waitfunc(ps)
			ps=[]
waitfunc(ps)
ps=[]
print('Running sagecal.')
for s in range(nsb):
	csb = sb + s
	for i in range(lstart,lend+1,3):
		cmdstr='sagecal -d '+wd+'L'+str(i)+'_SB%03d'%csb+ending+' -s '+skymodel+' -c '+cluster_file+' -p '+wd+'L'+str(i)+'_SB%03d_'%csb+'sagecal_solutions.txt -I DATA -O CORRECTED_DATA -t 60 -j 5 -b 1 -F 0'
		print(cmdstr)
		p = Popen(cmdstr, shell=True)
		ps.append(p)
		if(len(ps) == maxthread):
			waitfunc(ps)
			ps=[]

waitfunc(ps)
print('Process Complete.')
