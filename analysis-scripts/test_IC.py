#!/usr/bin/env python
"""
Loop over all fields in a given initial condition file and make 
sure none of the values are nan or weird.
"""

from netCDF4 import Dataset
import numpy as np
import pylab as pl

f=Dataset('/global/cfs/cdirs/e3sm/bhillma/scream/data/init/screami_ne4np4L128_20220512.nc')
#'/global/cfs/cdirs/e3sm/bhillma/scream/data/init/screami_ne120np4L72_20220503.nc')
#'/global/cfs/cdirs/e3sm/bhillma/scream/data/init/screami_ne30np4L72_20220503.nc')
#'/global/cfs/cdirs/e3sm/bhillma/scream/data/init/screami_ne512np4L128_20220505.nc')

vars=f.variables.keys()
for var in vars:
    x=f.variables[var]
    xmin=np.min(x)
    xmax=np.max(x)
    sz=x.shape
    print(var+str(sz)+': min=%f, max=%f'%(xmin,xmax))

