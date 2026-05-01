#%% Import modules:
import os

import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits import mplot3d

import readphilips as rp

# Working directory:
location = os.getcwd()
    
#%% 2D spiral data:
    
# Read in the raw-list files:
filename = location + '\\data\\2DSpiralVentilationCCHMC\\' + 'raw_004.data' 

# Read in the data using the ReadPhilips.py script:
inputfile = rp.PhilipsData(filename)
inputfile.compute()

# Extract data:
data = abs(inputfile.data)

# Remove singleton dimensions
data = data.squeeze()

# Plot an FID to check data:
fid = data[0,1,1,:] # First dynamic, first slice, first spiral
plt.figure()
plt.plot(fid)
plt.ylabel('k-space intensity')
plt.xlabel('Sample number')
plt.title('Free induction decay')
plt.show()

# Read in the raw-lab-sin files:
filename = location + '\\data\\2DSpiralVentilationCCHMC\\' + '20211013_113717_CPIR_Vent_HANNING_2DSOS_WIP.sin' 

# Read in the data using the ReadPhilips.py script:
inputfile =  rp.PhilipsData(filename)
inputfile.compute()

# Extract spiral coordinates:
coords = inputfile.spparams.get('COORDS')

# Plot the coordinates:
plt.figure()
plt.plot(coords[0,:], coords[1,:])
plt.ylabel('$k_y$')
plt.xlabel('$k_x$')
plt.title('Spiral coordinates')
plt.show()

# Extract expanded spiral coordinates:
coords_expanded = inputfile.spparams.get('COORDS_EXPANDED')

# Plot the coordinates:
plt.figure()
plt.plot(coords_expanded[0,:,:], coords_expanded[1,:,:])
plt.ylabel('$k_y$')
plt.xlabel('$k_x$')
plt.title('Spiral coordinates')
plt.show()

# Extract data:
data = abs(inputfile.data)

# Remove singleton dimensions
data = data.squeeze()

# Plot an FID to check data:
fid = data[0,4,0,:] # First dynamic, 5th slice, first spiral
plt.figure()
plt.plot(fid)
plt.ylabel('k-space intensity')
plt.xlabel('Sample number')
plt.title('Free induction decay')
plt.show()

#%% 3D radial data:
    
# Read in the raw-list files:
filename = location + '\\data\\3DRadialGas-exchangeCCHMC\\' + 'raw_007.list' 
# Note: the Duke data is in the correct radial format, whereas the CCHMC data 
# has had some editing so the data comes out strange.

# For future reference, CCHMC data looks like this:
    # "Approximate translation" from complicated Philips settings to MRI:
    # Dim 1 = mixes; used to separate spectroscopy data and imaging data.
    # Dim 2 = dynamics; used to separate gas, dissolved, and off-res acquisitions.
    # Dim 3 = kz; Archimedean spiral interleaves.
    # Dim 4 = ky; projections per interleave.
    # Dim 5 = kx; read out.
    # Useful paper: https://onlinelibrary.wiley.com/doi/full/10.1002/mrm.22898

# Read in the data using the ReadPhilips.py script:
inputfile =  rp.PhilipsData(filename)
inputfile.compute()

# Extract data:
data = np.squeeze(inputfile.data)

# Plot an FID to check data:
fid = abs(data[0,0,0,0,0:58])

plt.figure()
plt.plot(fid)
plt.ylabel('k-space intensity')
plt.xlabel('Sample number')
plt.title('Free induction decay')
plt.show()

# Read in the raw-lab-sin files:
location = os.getcwd()
filename = location + '\\data\\3DRadialGas-exchangeDuke\\' + '20220112_113653_DukeIPF_Gas_Exchange.sin'

# Read in the data using the ReadPhilips.py script:
inputfile =  rp.PhilipsData(filename)
inputfile.readParamOnly = False # Read raw as well
inputfile.trajtype = 2
inputfile.delay = 1.25
inputfile.compute()

# Extract radial coordinates:
coords = inputfile.radparams.get('COORDS')

# Plot the coordinates:
plt.figure()
ax = plt.axes(projection='3d')
ax.scatter(coords[:,:,-1,0], coords[:,:,-1,1], coords[:,:,-1,2], c = 'r', s = 1)
ax.set_zlabel('$k_z$')
ax.set_ylabel('$k_y$')
ax.set_xlabel('$k_x$')
plt.title('Radial coordinates')
plt.show()

# Extract data:
data = abs(inputfile.data)

# Remove singleton dimensions
data = data.squeeze()

# Plot an FID to check data:
fid = data[0,0,:] # First dyn, first proj
plt.figure()
plt.plot(fid)
plt.ylabel('k-space intensity')
plt.xlabel('Sample number')
plt.title('Free induction decay')
plt.show()

#%% 3D FLORET data:
    
# Read in the raw-lab-sin files:
filename = location + '\\data\\3DFloretVentilationCCHMC\\' + '20240321_114557_CPIR_VENT_FLORET.sin' 

# Read in the data using the ReadPhilips.py script:
inputfile =  rp.PhilipsData(filename)
inputfile.compute()

# Extract data:
data = abs(inputfile.data)

# Remove singleton dimensions
data = data.squeeze()

# Plot an FID to check data:
fid = data[0,0,:] # First interleave, first proj
plt.figure()
plt.plot(fid)
plt.ylabel('k-space intensity')
plt.xlabel('Sample number')
plt.title('Free induction decay')
plt.show()

# Extract spiral coordinates:
coords = inputfile.spparams.get('COORDS_EXPANDED')

# Plot some coordinates to check trajectory:
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
for proj in (0,26):
    for intlv in (0,1,2):
        plt.plot(coords[0,:,proj,intlv],coords[1,:,proj,intlv],coords[2,:,proj,intlv])
plt.ylabel('ky')
ax.view_init(elev=10,azim=10)
plt.xlabel('kx')
ax.set_zlabel('kz')
plt.title('Base FLORET Spiral')
plt.show()


# %%
