# Function that will wait until all processes have finished
from subprocess import Popen
from time import sleep
maxthread = 24

nsb=10
sb=120
scanlist=['260545','260553','260561','260569','260577','260585','260593','260601','260609','260617','260625','260633','260641','260649']
lampcal='260591'
fluxcaldir='/mnt/data/coughlan/exo/wd/Pre-Facet-Cal-exo1/'
targetdir='/mnt/data/coughlan/exo/initial_imaging/'
calending='_uv.dppp.ndppp_prep_cal'
ending='_uv.dppp.MS'
parmpath='parmexportcal'
sd='/mnt/home_cr/coughlan/lofar/scripts/bbs_dp3/'

def waitfunc(ps):
	while True:
  		ps_status = [p.poll() for p in ps]
		if all([x is not None for x in ps_status]):
  			break
		sleep(1)	# Only check every second
	return(0)

ps=[]
for i in range(nsb):
	csb = sb + i
	fn = fluxcaldir+'L'+str(lampcal)+'_SB%03d'%csb+calending
	cmdstr=parmpath+' in='+fn+'/instrument out='+targetdir+'L'+str(lampcal)+'_SB%03d_gain_solutions.table'%csb
	print(cmdstr)
	p = Popen(cmdstr, shell=True)
	ps.append(p)
	if(len(ps)==maxthread):
		waitfunc(ps)
		ps=[]
waitfunc(ps)


# Apply amplitude calibration
# Use 1 thread per each of the 8 scans = 8 threads
print('Applying amplitude calibration to all scans/subbands')

ps=[]
for s in range(nsb):
	csb = sb + s
	for i in scanlist:
		cmdstr='calibrate-stand-alone -t 1 --parmdb '+targetdir+'L'+str(lampcal)+'_SB%03d_gain_solutions.table '%csb+targetdir+'L'+str(i)+'_SB%03d'%csb+ending+' '+sd+'correct_only.bbs'
		print(cmdstr)
		p = Popen(cmdstr, shell=True)
		ps.append(p)
		if(len(ps)==maxthread):
			waitfunc(ps)
			ps=[]

waitfunc(ps)
ps=[]
