#!/usr/bin/env python
import os
from lofarpipe.support.data_map import DataMap
from lofarpipe.support.clusterdesc import ClusterDesc, get_compute_nodes
import itertools as it


# Colm Coughlan, DIAS, May 2016
# Change hostnames in mapfiles to allow for MPI execution on Fionn (and probably other clusters with shared filesystems like Fionn)


def plugin_main(args, **kwargs):
    """
    Takes in mapfile and change host names to allow for efficient MPI reduction

    Parameters
    ----------
    mapfile : str
        Name of the input mapfile. WILL BE MODIFIED!
    mapfile_dir : str
        Name of the directory containing the mapfile

    Returns
    -------
    result : empty dictionary

    """

    result = {}
    mapfile = kwargs['mapfile']
    mapfile_dir = kwargs['mapfile_dir']
    fn = os.path.join(mapfile_dir,mapfile)
    cn_cycle = it.cycle( get_compute_nodes( ClusterDesc(str(os.environ['cluster_desc_file'])) ) ) # Read in list of compute nodes. Set up iterator to cyclically iterate over them.

    data = DataMap.load(fn)                 # read in current data map file (probably with all host values set to "localhost")
    iterator = DataMap.SkipIterator(data)   # set up iterator for all values in mapfile

    for value in iterator:
        value.host = cn_cycle.next()   # iterate through map file, assigning each entry a host from the available compute nodes in a cyclical fashion
    data.save(fn)   # overwrite original file
    return result

