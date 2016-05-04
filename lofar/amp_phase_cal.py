#	LOFAR calibration pipeline


from os import environ
from subprocess import Popen
from time import sleep
from numpy import ceil

# Read in environment variables from script

lofar_home=str(environ['LOFARROOT'])+'/'
wd=str(environ['working_dir'])+'/'
sd=str(environ['script_dir'])+'/'
sb=int(environ['subband'])
nsb=int(environ['num_subbands'])
lstart=int(environ['l_start'])
lend=int(environ['l_end'])
lampcal=int(environ['l_ampcal'])

gsmskymodel=str(environ['gsmskymodel'])

parmpath=str(environ['parmpath'])

maxthread=int(environ['nthreads'])

# Function that will wait until all processes have finished

def waitfunc(ps):
	while True:
  		ps_status = [p.poll() for p in ps]
		if all([x is not None for x in ps_status]):
  			break
		sleep(1)	# Only check every second
	return(0)


#ending='_uv.dppp.CRS.MS'
ending='_uv.dppp.MS'


# Amplitude calibration (maxthread/nsb threads each)

ps=[]

use_threads = ceil(maxthread/nsb)
for i in range(nsb):
	csb = sb + i
	fn = wd+'L'+str(lampcal)+'_SB%03d'%csb+ending
	cmdstr = 'calibrate-stand-alone -f '+fn+' '+sd+'amp_cal.bbs ' +lofar_home+'share/pipeline/skymodels/3C147.skymodel > '+wd+'SB%03d_amp_cal.log'%csb
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

#ps=[]
#for i in range(nsb):
#	csb = sb + i
#	fn = wd+'L'+str(lampcal)+'_SB%03d'%csb+ending
#	file = open(wd+'SB%03d'%csb+'.parmdbm','w')
#	file.write('open tablename=\''+fn+'/instrument\'\n')
#	file.write('export Gain* tablename=\''+wd+'SB%03d_gain_clock_solutions.table'%csb+'\'\n')
#	file.write('export Clock* tablename=\''+wd+'SB%03d_gain_clock_solutions.table'%csb+'\',append=1\n')
#	file.write('close\n')
#	file.write('exit\n')
#	file.close()
#	cmdstr = 'parmdbm < '+wd+'SB%03d'%csb+'.parmdbm'
#	p = Popen(cmdstr,shell=True)
#	ps.append(p)
#waitfunc(ps)

# Note parmexportcal can be used to export the median of the gain solutions, but doesn't work for clock at the moment!
# Remeber to change solution interval from 0 to 1 in amp_cal to take advantage of parmexportcal's median if you use it
#
#
#
# Now using custom version of parmexportcal that does export the clock solutions
#

for i in range(nsb):
	csb = sb + i
	fn = wd+'L'+str(lampcal)+'_SB%03d'%csb+ending
	cmdstr=parmpath+' in='+fn+'/instrument out='+wd+'SB%03d_gain_clock_solutions.table'%csb
	print(cmdstr)
	p = Popen(cmdstr, shell=True)
	ps.append(p)
	if(len(ps)==maxthread):
		waitfunc(ps)
		ps=[]
waitfunc(ps)
ps=[]


# Apply amplitude calibration
# Use 1 thread per each of the 8 scans = 8 threads
print('Amplitude calibration complete on all scans for subband '+str(sb))


for s in range(nsb):
	csb = sb + s
	for i in range(lstart,lend+1,3):
		cmdstr='calibrate-stand-alone -t 1 --parmdb '+wd+'SB%03d_gain_clock_solutions.table '%csb+wd+'L'+str(i)+'_SB%03d'%csb+ending+' '+sd+'correct_only.bbs > '+wd+'L'+str(i)+'_SB%03d'%csb+'_apply_gain_clock_solutions.log'
		print(cmdstr)
		p = Popen(cmdstr, shell=True)
		ps.append(p)
		if(len(ps)==maxthread):
			waitfunc(ps)
			ps=[]

waitfunc(ps)
ps=[]

# BBS based direction independent phase calibration using GSM


for s in range(nsb):
	csb = sb + s
	for i in range(lstart,lend+1,3):
		cmdstr = 'calibrate-stand-alone -t 1 -f '+wd+'L'+str(i)+'_SB%03d'%csb+ending+' '+sd+'phase_cal.bbs '+gsmskymodel+' > '+wd+'L'+str(i)+'_SB%03d_phase_cal.log'%csb
		print(cmdstr)
		p = Popen(cmdstr, shell=True)
		ps.append(p)
		if(len(ps)==maxthread):
			waitfunc(ps)
			ps=[]

waitfunc(ps)




print('Calibration run complete')
