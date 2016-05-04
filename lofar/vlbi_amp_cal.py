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

casa_exec=str(os.environ['casa_exec'])
mscorpol_exec='/mnt/home_cr/coughlan/lofar/scripts/mscorpol/mscorpol.py'

parmpath=str(os.environ['parmpath'])

io_thread=1
maxthread=int(os.environ['nthreads'])

scanlist = range(lstart , lend+1, 3)

sbs = range(sb,sb+nsb,1)
print(str(len(sbs)))
bad_sb=[112,185]	# remove bad sb for dgtau data
sb_per_if = 16
nif = int(ceil(float(nsb)/float(sb_per_if)))

print('\n\n'+str(nif)+' IFs will be created''\n\n')

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



# Flag bad data
'''
print 'Flagging bad data in target field (AOFlagger)'
for s in range(nsb):
	csb = sb + s
	for i in scanlist:
		cmdstr = 'aoflagger -strategy /mnt/home_cr/coughlan/lofar/scripts/flagging_strategy_ttau_new.rfis -column DATA '+wd+od+'L'+str(i)+'_SB%03d'%csb+ending+' > '+wd+log_dir+'L'+str(i)+'_SB%03d'%csb+'.AOflagger.log'
		print(cmdstr)
		p = Popen(cmdstr, shell=True)
		ps.append(p)
		waitfunc(ps)
		ps=[]
'''
# Flag bad stations in target field
'''
print 'Flagging bad stations in target field'
for s in range(nsb):
	csb = sb + s
	for i in scanlist:
		cmdstr = 'NDPPP '+sd+'flag_bad_stations.ndppp msin='+wd+od+'L'+str(i)+'_SB%03d'%csb+ending+' msout='+wd+od+'L'+str(i)+'_SB%03d'%csb+ending+' > '+wd+log_dir+'L'+str(i)+'_SB%03d'%csb+'_flag_bad_stations.ndppp.log'
		print(cmdstr)
		p = Popen(cmdstr, shell=True)
		ps.append(p)
		if(len(ps)==maxthread):
			waitfunc(ps)
			ps=[]
waitfunc(ps)
ps=[]
'''

# Amplitude calibration (maxthread/nsb threads each)

ps=[]

if do_cal:
# Flag bad stations in amp calibrator

#	ps=[]
#	print 'Flagging bad data in target field (AOFlagger)'
#	for s in range(nsb):
#		csb = sb + s
#		fn = wd+od+'L'+str(lampcal)+'_SB%03d'%csb+ending
#		cmdstr = 'aoflagger -strategy /mnt/home_cr/coughlan/lofar/scripts/flagging_strategy_ttau_new.rfis -column DATA '+fn+' > '+wd+log_dir+'L'+str(lampcal)+'_SB%03d'%csb+'.AOflagger.log'
#		print(cmdstr)
#		p = Popen(cmdstr, shell=True)
#		ps.append(p)
#		waitfunc(ps)
#		ps=[]

	for csb in sbs:
		if(not(csb in bad_sb)):
			fn = wd+od+'L'+str(lampcal)+'_SB%03d'%csb+ending
			cmdstr = 'calibrate-stand-alone -f '+fn+' '+sd+'vlbi_amp_cal.bbs ' +lofar_home+'share/pipeline/skymodels/3C147.skymodel > '+wd+log_dir+'SB%03d.amp_cal.ndppp.log'%csb
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

	for csb in sbs:
		if(not(csb in bad_sb)):
			fn = wd+od+'L'+str(lampcal)+'_SB%03d'%csb+ending
			cmdstr=parmpath+' in='+fn+'/instrument out='+wd+table_dir+'L'+str(lampcal)+'_SB%03d_gain_solutions.table'%csb+' > '+wd+log_dir+'SB%03d.parmexportal.log'%csb
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
'''
ps=[]
for csb in sbs:
	if(not(csb in bad_sb)):
		for i in scanlist:
			cmdstr='calibrate-stand-alone -t 1 --parmdb '+wd+table_dir+'L'+str(lampcal)+'_SB%03d_gain_solutions.table '%csb+wd+od+'L'+str(i)+'_SB%03d'%csb+ending+' '+sd+'vlbi_correct_only.bbs > '+wd+log_dir+'L'+str(i)+'_SB%03d'%csb+'_apply_gain_solutions.log'
			print(cmdstr)
			p = Popen(cmdstr, shell=True)
			ps.append(p)
			if(len(ps)==maxthread):
				waitfunc(ps)
				ps=[]

waitfunc(ps)

'''
ps=[]

if not os.path.exists(data_dir):
	cmdstr = 'mkdir '+data_dir
	p = Popen(cmdstr, shell=True)
	ps.append(p)
	if(len(ps)==1):
		waitfunc(ps)
		ps=[]
else:
	print('Warning! '+data_dir+' detected. Please remove/rename old runs before starting a new one.')
	exit()

print('\n\nMaking superstation and removing CS from result in single step\n\n')
ps=[]
ifnum=1
for csb in sbs:
	if(csb-sbs[0]+1 > ifnum*sb_per_if):
		ifnum = ifnum+1
	if(not(csb in bad_sb)):
		for i in scanlist:
			cmdstr = 'NDPPP '+sd+'vlbi_make_superstation_multi.ndppp msin='+wd+od+'L'+str(i)+'_SB%03d'%csb+ending+' msin.datacolumn=CORRECTED_DATA msout='+wd+data_dir+'L'+str(i)+'_SB%03d_IF'%csb+str(ifnum)+ending+' > '+wd+log_dir+'L'+str(i)+'_SB%03d_IF'%csb+str(ifnum)+'_make_superstation_multi.ndppp.log'
			print(cmdstr)
			p = Popen(cmdstr, shell=True)
			ps.append(p)
			if(len(ps)==io_thread):
				waitfunc(ps)
				ps=[]

waitfunc(ps)
ps=[]


# Merge the data by sub-band, taking sb_per_if into account. Construct a list to include missing data filenames, if any.


ps=[]
for aips_if in range(nif):
	for i in scanlist:
		ms_list=['[']
		for csb in range(sb+aips_if * sb_per_if,sb + (aips_if+1) * sb_per_if):
			ms_list.append(wd+data_dir+'L'+str(i)+'_SB%03d_IF'%csb+str(aips_if+1)+ending)
			if csb != sb + (aips_if+1) * sb_per_if - 1:
				ms_list.append(',')
		ms_list.append(']')
		cmdstr = 'NDPPP '+sd+'vlbi_merge.ndppp msin='+''.join(ms_list)+' msout='+wd+'L'+str(i)+'_IF'+str(aips_if+1)+'.ms'+' > '+wd+log_dir+'L'+str(i)+'_IF'+str(aips_if+1)+'.ndppp.log'
		print(cmdstr)
		p = Popen(cmdstr, shell=True)
		ps.append(p)
		if(len(ps)==io_thread):
			waitfunc(ps)
			ps=[]
waitfunc(ps)
ps=[]

print '\n\nConcatenating by time\n\n'



for aips_if in range(nif):
	ps=[]
	concat_list=[]
	for i in scanlist:
		concat_list.append(wd+'L'+str(i)+'_IF'+str(aips_if+1)+'.ms')
	fn='target_IF'+str(aips_if+1)
	casa_concat_parset(concat_list, fn)
	cmdstr=casa_exec+' -c '+fn+'.py'
	print(cmdstr)
	os.system(cmdstr)


print('\n\nFlagging bad new UVW points(if necessary)\n\n')
ps=[]
for aips_if in range(nif):
	fn='target_IF'+str(aips_if+1)+'.ms'
	cmdstr = 'taql'+' \'update '+fn+' set FLAG=T where isnan(UVW[0])\' > '+wd+log_dir+'target_IF'+str(aips_if+1)+'_flag_bad_uv_coords.taql.log'
	print(cmdstr)
	p = Popen(cmdstr, shell=True)
	ps.append(p)
	if(len(ps)==maxthread):
		waitfunc(ps)
		ps=[]

waitfunc(ps)
ps=[]
	
	
# Now convert from linear to circular polarisation and export

print '\n\nConverting from linear to circular polarisation and applying beam correction. Exporting result as FITS file.\n\n'

ps=[]
for aips_if in range(nif):
	fn='target_IF'+str(aips_if+1)+'.ms'
	cmdstr = mscorpol_exec+' -f '+fn+' > '+wd+log_dir+'target_IF'+str(aips_if+1)+'.mscorpol.log'
	print(cmdstr)
	p = Popen(cmdstr, shell=True)
	ps.append(p)
	if(len(ps)==maxthread):
		waitfunc(ps)
		ps=[]

waitfunc(ps)
ps=[]

for aips_if in range(nif):
	fn='target_IF'+str(aips_if+1)
	cmdstr='ms2uvfits ms='+fn+'.ms fitsfile='+fn+'.fits column=DATA writesyscal=FALSE > '+wd+log_dir+'target_IF'+str(aips_if+1)+'.ms2uvfits.log'
	print(cmdstr)
	os.system(cmdstr)



print(str(nif)+' IFs calibrated and prepared for AIPS.')



