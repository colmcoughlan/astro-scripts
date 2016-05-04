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
'''
ps=[]
#badsb=[330,331,332,333,334,340,347,361,362,367]
badsb=[361,362,367]
for s in badsb:
	for i in range(lstart,lend+1,3):
		fn = od+'L'+str(i)+'_SB%03d'%s+ending
		cmdstr='NDPPP '+sd+'flag_bad_channels.ndppp msin='+fn+' msout='+fn+' > '+od+'L'+str(i)+'_SB%03d'%sb+ending+'.avg10.flagchan.log'
		print(cmdstr)
		p = Popen(cmdstr, shell=True)
		ps.append(p)
		if(len(ps)==maxthread):
			waitfunc(ps)
			ps=[]
'''
ps=[]
for s in range(nsb):
	csb = sb + s
	for i in range(lstart,lend+1,3):
		cmdstr='NDPPP '+sd+'average_10s.dppp msin='+od+'L'+str(i)+'_SB%03d'%csb+ending+' msout='+wd+'L'+str(i)+'_SB%03d'%csb+ending+' > '+od+'L'+str(i)+'_SB%03d'%csb+ending+'.avg10.ndppp.log'
		print(cmdstr)
		p = Popen(cmdstr, shell=True)
		ps.append(p)
		if(len(ps)==maxthread):
			waitfunc(ps)
			ps=[]

waitfunc(ps)
ps=[]

