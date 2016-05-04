#	LOFAR calibration pipeline


from os import environ
from subprocess import Popen
from time import sleep
from numpy import ceil

# Read in environment variables from script

lofar_home=str(environ['LOFARROOT'])+'/'
od=str(environ['o_dir'])+'/'
wd=str(environ['working_dir'])+'/'
sd=str(environ['script_dir'])+'/'
sb=int(environ['subband'])
nsb=int(environ['num_subbands'])
lstart=int(environ['l_start'])
lend=int(environ['l_end'])
maxthread=int(environ['nthreads'])

ending='_uv.dppp.CRS.MS'

# Function that will wait until all processes have finished

def waitfunc(ps):
	while True:
  		ps_status = [p.poll() for p in ps]
		if all([x is not None for x in ps_status]):
  			break
		sleep(1)	# Only check every second
	return(0)


#print('Flagging LB with NDPPP.')
# Flag and filter long baselines on target and calibrator fields

ps=[]
for s in range(nsb):
	csb = sb + s
	for i in range(lstart,lend+1,3):
		cmdstr='NDPPP '+sd+'average_10s.dppp msin='+od+'L'+str(i)+'_SB%03d'%csb+ending+' msout='+wd+'L'+str(i)+'_SB%03d_uv.dppp.CRS.avg10.MS'%csb
		print(cmdstr)
		p = Popen(cmdstr, shell=True)
		ps.append(p)
		if(len(ps)==maxthread):
			waitfunc(ps)
			ps=[]

waitfunc(ps)
ps=[]

ps=[]
for s in range(nsb):
	csb = sb + s
	for i in range(lstart,lend+1,3):
		cmdstr='NDPPP '+sd+'make_corrected_column.dppp msin='+wd+'L'+str(i)+'_SB%03d_uv.dppp.CRS.avg10.MS'%csb+' msout='+wd+'L'+str(i)+'_SB%03d_uv.dppp.CRS.avg10.MS'%csb
		print(cmdstr)
		p = Popen(cmdstr, shell=True)
		ps.append(p)
		if(len(ps)==maxthread):
			waitfunc(ps)
			ps=[]

waitfunc(ps)
ps=[]

