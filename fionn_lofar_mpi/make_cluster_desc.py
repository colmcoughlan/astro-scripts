# Colm Coughlan, DIAS
# This file prepares a cluster description for LOFAR pipelines using information from Fionn

import os

cd_name =str(os.environ['cluster_desc_file'])
localdisk = '/ichec/work/dsast016b'

try:
	pbs_fn = str(os.environ['PBS_NODEFILE'])	# read cluster information
except:
	print('Error: PBS_NODEFILE not defined')
	exit()

with open(pbs_fn) as f:
	node_info = f.readlines()

head_node = node_info[0].rstrip('\n')	# set head node to first computer node (I don't think the head node is used much!)

nodes = set(node_info)	# get set of nodes
nodes = (list(nodes))
nnodes = len(nodes)
for node in range(nnodes):
	nodes[node] = nodes[node].rstrip('\n')


with open(cd_name, 'w') as outfile:
	outfile.write('ClusterName = fionn\n')
	outfile.write('\n')
	outfile.write('# Compute nodes (varies with number of cores used in qsub)\n')
	outfile.write('Compute.Nodes = [')
	for node in range(nnodes):
		outfile.write(nodes[node])
		if(node + 1 != nnodes):
			outfile.write(',')
	outfile.write(']\n')
	outfile.write('Compute.LocalDisks = [ '+localdisk+' ] \n')
	outfile.write('\n')
	outfile.write('# Head nodes (using compute node as head node)\n')
	outfile.write('Head.Nodes = [ '+head_node+' ] \n')
	outfile.write('Head.LocalDisks = [ '+localdisk+' ] \n')

print('Cluster description written to '+cd_name)
