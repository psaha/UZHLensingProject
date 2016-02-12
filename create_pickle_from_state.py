#!/usr/bin/env python 
# -*- coding: utf-8 -*-
"""
Converts a glass state file (*.state) to a pickles file used by LO and her
parameter fitting program


Use either on command line:
    ./process_state model [ext]
    
    loads the statefile model.ext (where ext defaults to ".state")
    and results in an pickle file "model.pickle"
 

Or use as module:
    import process_state as procstate
    
then either create pickle:
    procstate.mk_pkl(model, ext=".state")
    
    model:  string that points to the state file, without extension
            ("ASW0001234/HASHOFMODEL")
    ext:    extension of the state file (default='.state')
    
or only use the data loader:
    data = procstate.get_data_from_state(statefn)
    
    statefn: string that points to the state file, with extension
             ("ASW0001234/HASHOFMODEL.state")
    data:    returnvalue, dict of data


!!! an install of glass is required! make sure it's in the python path!
!!! best to use interactive glass and module:

$ /path/to/glass/interactive_glass
[1] import process_state as ps
[2] ps


rk 2016.02.05: Heavily modified to either work on the command line or as an
               module

modified p sahas ensem to work from command line with args

Created on Wed Jul  8 12:57:52 2015
Modified on Fri Feb  5 14:29:35 2016 rafik
@author: psaha, rafik
"""

import sys
import os
import glob
import pickle



def mk_pkl(model, ext = ".state"):
    """ Main file, loads a state file, specified by a path without extension
    (model), optionally give [ext]enstion as well.
    Saves a resulting pickles file with filename [model].pkl in the same place
    """

    statefn = model + ext
    chk_filename(statefn)
    
    pklfn = model + '.pkl'
    d = get_data_from_state(statefn)
    save_as_pickle(d, pklfn)


def mk_pkl2(fpath):
    """ accept full path with extension """
    model,ext = os.path.splitext(fpath)
    mk_pkl(model, ext)


def get_data_from_state(statefn):
    """ loads the state file, creates and returns the data dictionary"""
    
    #check if the file really exists
    if not os.path.isfile(statefn):
        print "no valid file! (%s)" % statefn
        sys.exit('no valid file')
    
    state = loadstate(statefn)
    
    grids = []
    for m in state.models:
        obj,data = m['obj,data'][0]
        g = obj.basis._to_grid(data['kappa DM'])
        grids.append(g)
        
    imgradii = [abs(img.pos) for img in state.objects[0].sources[0].images]
    
    # create dict to be pickled afterwards
    d = {
        'grids': grids,
        'maprad': obj.basis.maprad,
        'pixrad': obj.basis.pixrad,
        'maximgpos': max(imgradii),
    }

    return d



def save_as_pickle(data, filename):
    # save as pickle
    with open(filename,'w') as fil:
        pickle.dump(data,fil)
        
    print "done (saved as '%s')" % filename
    

def chk_filename(statefn):

    if not os.path.isfile(statefn):
        print "no valid file! (%s)" % statefn
        sys.exit('no valid file')






def create_conv_for_paper():
    """ shorthand funcion for creating the pickles for the paper """

    state_files =  glob.glob('for_paper/*/*/*.state')
    
    for fpath in state_files:
        path, statefn = os.path.split(fpath)
        mid, ext = os.path.splitext(statefn)
        _, tjpe, asw = path.split(os.sep)
        print "working on:", asw, tjpe, mid
        
        mk_pkl2(fpath)
        



# main routine. is only run if executed from command line
# otherwise use mk_pkl()
if __name__ == "__main__":

    # read first parameter from command line (remember: 0 is the script file name)
    try: 
        model = sys.argv[1]
    except IndexError:
        print "no file given as first argument. give a path to a state file WITHOUT extension"
        sys.exit('no argument')
    
    try: 
        model = sys.argv[1]
    except IndexError:
        ext = '.state'
    
    mk_pkl(model, ext)






