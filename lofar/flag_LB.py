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
lampcal=int(environ['l_ampcal'])
maxthread=int(environ['nthreads'])
aoflagger=str(environ['aoflagger_cmd'])
aostrategy=str(environ['aostrategy'])

ending='_uv.dppp.MS'

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
		cmdstr='NDPPP '+sd+'flag_LB.ndppp msin='+od+'L'+str(i)+'_SB%03d'%csb+ending+' msout='+wd+'L'+str(i)+'_SB%03d'%csb+ending+' > '+wd+'L'+str(i)+'_SB%03d'%csb+ending+'.flag_LB.ndppp.log'
		print(cmdstr)
		p = Popen(cmdstr, shell=True)
		ps.append(p)
		if(len(ps)==maxthread):
			waitfunc(ps)
			ps=[]

waitfunc(ps)
ps=[]


for i in range(nsb):
	csb = sb + i
	fn = 'L'+str(lampcal)+'_SB%03d'%csb+ending
	cmdstr='NDPPP '+sd+'flag_LB.ndppp msin='+od+fn+' msout='+wd+fn+' > '+wd+fn+'.flag_LB.log'
	print(cmdstr)
	p = Popen(cmdstr,shell=True)
	ps.append(p)
	if(len(ps)==maxthread):
		waitfunc(ps)
		ps=[]
waitfunc(ps)
ps=[]


# Flag target field

ps=[]
for s in range(nsb):
	csb = sb + s
	for i in range(lstart,lend+1,3):
		cmdstr=aoflagger+' -strategy '+aostrategy+' -j '+str(maxthread)+' '+wd+'L'+str(i)+'_SB%03d'%csb+ending+' > '+wd+'L'+str(i)+'_SB%03d'%csb+'_aoflagger.log'
		print(cmdstr)
		p = Popen(cmdstr, shell=True)
		ps.append(p)
		if(len(ps)!=0):
			waitfunc(ps)
			ps=[]

waitfunc(ps)

ps=[]
for s in range(nsb):
	csb = sb + s
	cmdstr=aoflagger+' -strategy '+aostrategy+' -j '+str(maxthread)+' '+wd+'L'+str(lampcal)+'_SB%03d'%csb+ending+' > '+wd+'L'+str(lampcal)+'_SB%03d'%csb+'_aoflagger.log'
	print(cmdstr)
	p = Popen(cmdstr, shell=True)
	ps.append(p)
	if(len(ps)!=0):
		waitfunc(ps)
		ps=[]

waitfunc(ps)

ps=[]
for s in range(nsb):
	csb = sb + s
	for i in range(lstart,lend+1,3):
		cmdstr='NDPPP '+sd+'flag_LB.ndppp msin='+wd+'L'+str(i)+'_SB%03d'%csb+ending+' msout='+wd+'L'+str(i)+'_SB%03d'%csb+ending+' > '+wd+'L'+str(i)+'_SB%03d'%csb+ending+'.flag_LBr2.ndppp.log'
		print(cmdstr)
		p = Popen(cmdstr, shell=True)
		ps.append(p)
		if(len(ps)==maxthread):
			waitfunc(ps)
			ps=[]

waitfunc(ps)
ps=[]

ps=[]
for i in range(nsb):
	csb = sb + i
	fn = 'L'+str(lampcal)+'_SB%03d'%csb+ending
	cmdstr='NDPPP '+sd+'flag_LB.ndppp msin='+wd+fn+' msout='+wd+fn+' > '+wd+fn+'.flag_LBv2.log'
	print(cmdstr)
	p = Popen(cmdstr,shell=True)
	ps.append(p)
	if(len(ps)==maxthread):
		waitfunc(ps)
		ps=[]
waitfunc(ps)
ps=[]

