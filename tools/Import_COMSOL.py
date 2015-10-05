# -*- coding: utf-8 -*-
"""
Created on Mon Sep 14 14:06:41 2015
Import_COMSOL.py: Generate lookup tables from COMSOL output.
@author: jkk3
"""

import numpy as np
import sys

fname = sys.argv[1]

data = np.loadtxt(fname, skiprows = 8, delimiter = ',')
#data = np.loadtxt(fname, skiprows = 7, delimiter = ',', dtype = 'string')

np.save(fname.split('.')[0] + '_x.npy', data[:,0]*1E6)
np.save(fname.split('.')[0] + '_c.npy', data[:,1::2])
np.save(fname.split('.')[0] + '_dc.npy', data[:,2::2])