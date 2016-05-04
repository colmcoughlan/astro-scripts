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

gsmskymodel=str(os.environ['gsmskymodel'])
casa_exec=str(os.environ['casa_exec'])

parmpath=str(os.environ['parmpath'])

maxthread=int(os.environ['nthreads'])

scanlist = [192760]

# Function that will wait until all processes have finished

def waitfunc(ps):
	while True:
  		ps_status = [p.poll() for p in ps]
		if all([x is not None for x in ps_status]):
  			break
		sleep(1)	# Only check every second
	return(0)

# Function to write casa parset to concat and export files for AIPS

def casa_concat_parset(vis, concatvis):
	try:
		f = open(concatvis+'.py','w')
	except:
		print "\t Error opening ", concatvis
		return(-1)


	f.write('concat(vis='+str(vis)+', concatvis = \''+concatvis+'.ms\', async = false)\n')

	f.close
	return(0)

def casa_exportuvfits_parset(vis, fitsfile):
	try:
		f = open(fitsfile+'.py','w')
	except:
		print "\t Error opening ", vis
		return(-1)


	f.write('exportuvfits(vis = \''+str(vis)+'.ms\', fitsfile = \''+fitsfile+'.fits\', async = false)\n')

	f.close
	return(0)

# set the default ending for LOFAR files

ending='_uv.dppp.MS'

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
do_cal = False

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

od = 'data/'

if not os.path.exists(od):
	print(od+' folder not detected. Assuming current directory is a suitable working directory and contains the data.')
	odir=''
else:
	print(od+' folder detected. Assuming data is stored within and it is a suitable working directory.')

data_dir = 'data_final/'


print('Making superstation and removing CS from result in single step...')
ps=[]
for s in range(nsb):
	csb = sb + s
	for i in scanlist:
		cmdstr=cmdstr = 'NDPPP '+sd+'vlbi_make_superstation_multi.ndppp msin='+wd+od+'L'+str(i)+'_SB%03d'%csb+ending+' msout='+wd+data_dir+'L'+str(i)+'_SB%03d'%csb+ending+' > '+wd+log_dir+'L'+str(i)+'_SB%03d'%csb+'_make_superstation_multi.ndppp.log'
		print(cmdstr)
		p = Popen(cmdstr, shell=True)
		ps.append(p)
		if(len(ps)==maxthread):
			waitfunc(ps)
			ps=[]

waitfunc(ps)
ps=[]



# Merge the data by sub-band

ps=[]
for i in scanlist:
	cmdstr = 'NDPPP '+sd+'vlbi_merge.ndppp msin='+wd+data_dir+'\'*'+str(i)+'*\''+' msout='+wd+'L'+str(i)+'_merged.ms'+' > '+wd+log_dir+'L'+str(i)+'_merged.ndppp.log'
	print(cmdstr)
	p = Popen(cmdstr, shell=True)
	ps.append(p)
	if(len(ps)==maxthread):
		waitfunc(ps)
		ps=[]

ps=[]
concat_list=[]
for i in scanlist:
	concat_list.append(wd+'L'+str(i)+'_merged.ms')
fn='3C147_concat'
casa_concat_parset(concat_list, fn)
cmdstr=casa_exec+' -c '+fn+'.py'
print(cmdstr)
p = Popen(cmdstr, shell=True)
ps.append(p)
if(len(ps)==1):
	waitfunc(ps)
	ps=[]


# Now convert from linear to circular polarisation

print 'Converting from linear to circular polarisation'
print 'Updating polarisation table'

fn='3C147_concat.ms'
cmdstr = 'python '+sd+'lin2circ.py -i '+fn+' -p True -c DATA -o CORRECTED_DATA'
print(cmdstr)
p = Popen(cmdstr, shell=True)
ps.append(p)
if(len(ps)==1):
	waitfunc(ps)
	ps=[]

# Now export for AIPS


print 'Concatenating by time and exporting for AIPS'
ps=[]
fn='3C147_concat'
casa_exportuvfits_parset(fn, fn)
cmdstr=casa_exec+' -c '+fn+'.py'
print(cmdstr)
p = Popen(cmdstr, shell=True)
ps.append(p)
if(len(ps)==1):
	waitfunc(ps)
	ps=[]


print('Calibration run complete')
print('Final circular polarisation data in CORRECTED_DATA column')



