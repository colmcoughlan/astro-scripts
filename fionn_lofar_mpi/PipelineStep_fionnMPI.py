#!/usr/bin/env python
import os
from lofarpipe.support.data_map import DataMap
from lofarpipe.support.clusterdesc import ClusterDesc, get_compute_nodes, get_head_node
import itertools as it


# Colm Coughlan, DIAS, May 2016
# Change hostnames in mapfiles to allow for MPI execution on Fionn (and probably other clusters with shared filesystems like Fionn)


def plugin_main(args, **kwargs):
    """
    Takes in mapfiles and change host names to allow for efficient MPI reduction

    Parameters
    ----------
    mapfiles : list of strs
        List of the names of the input mapfiles. WILL BE MODIFIED!
    mapfile_dir : str
        Name of the directory containing the mapfile
    head_node_only : str
        String: Either True or False. Describes whether to use just the head node or not.

    Returns
    -------
    result : empty dictionary

    """

    result = {}
    mapfiles = (kwargs['mapfiles'][1:-1]).split(',')    # read in list of mapfiles from string (separated by commas)
    mapfile_dir = kwargs['mapfile_dir']
    head_only = (kwargs['head_node_only'] in ['True','true','T','t','1'])
    fn_list=[]
    for mf in mapfiles:
        fn_list.append( os.path.join(mapfile_dir,mf) )

    # caution: remember to reload the compute node iterable for every mapfile to ensure corresponding entries have the same node set as host

    for fn in fn_list:
        if(head_node_only):
            cn_cycle = it.cycle( get_head_node( ClusterDesc(str(os.environ['cluster_desc_file'])) ) )   # Read in head node. Set up iterator (unnessary with just one node, but better to have less code!)
        else:
            cn_cycle = it.cycle( get_compute_nodes( ClusterDesc(str(os.environ['cluster_desc_file'])) ) ) # Read in list of compute nodes. Set up iterator to cyclically iterate over them.

        data = DataMap.load(fn)                 # read in current data map file (probably with all host values set to "localhost")
        iterator = DataMap.SkipIterator(data)   # set up iterator for all values in mapfile

        for value in iterator:
            value.host = cn_cycle.next()   # iterate through map file, assigning each entry a host from the available compute nodes in a cyclical fashion
        data.save(fn)   # overwrite original file
    return result

