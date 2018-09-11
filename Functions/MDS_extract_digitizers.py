# -*- coding: utf-8 -*-
"""
MDS_extract.py

Author:  Michael Ross
Description:  This code extracts data from MDSplus and saves that data as a CSV
file for later use in MATLAB.  This enables MATLAB users to avoid the memory
leakage issues associated with reading MDSplus data directly into MATLAB.  
These memory issues typically cause MATLAB to crash.
"""

import MDSplus as MDS
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import *
import sys
import numpy as np

tree = 'zaphd'
server = 'zappa.zap'
prefix = '\\'

# Note that the inputs from the command line must be typecast from strings to numerical datatypes:
shotnum = np.float64(sys.argv[1]) # using float 64-bit necessary to extracting shotnums less than 010; if future problems arise, try type-casting to integer type
signame = sys.argv[2]
#shotnum = 160524007
#signame = 'y_p0'
#print "change python back"

# Using the MDSPlus interface to extract time series data:
mds = MDS.Connection(server)
mds.openTree(tree,shotnum)
data = mds.get(prefix + signame)
#time = mds.get('DIM_OF('+ prefix + signame +')')
mds.closeTree(tree,shotnum)

#print np.shape(data)
#print type(data)
#print np.shape(time)
#print type(time)

# pre-cast arrays into which data and time will be placed:
b = np.zeros((np.size(data),1))
#b[:,0] = time
b = data


# Saving the extracted data into a temporary CSV file:
np.savetxt("M:\Users\Michael_Ross\Coding\Functions\MDS_temp_digitizers.csv", b, delimiter=",")

#plt.figure(1)
#plt.plot(b[:,0],b[:,1])


"""
np.shape(time)
np.shape(data)
np.shape(data[:])
np.shape(time[:])
np.shape(a)
a
b = np.nan(np.zeros(2,np.size(time)))
AAAA = np.zeros(2,2)
b = np.zeros(2,np.size(time))
b = np.zeros((2,np.size(time)))
b
b[0,] = time
b[1,] = data
b"""
