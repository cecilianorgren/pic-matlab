Particle-in-cell simulation of magnetotail reconnection with cold plasma protons in inflow regions. Used for manuscript titled 'On the presence and thermalization of cold ions in the exhaust of antiparallel symmetric reconnection' by Norgren et al., 2021

Contact information:
Cecilia Norgren
Space Plasma Physics Group, University of Bergen, 5007 Bergen, Norway.
E-mail: cecilia.norgren@uib.no
ORCID ID of author: https://orcid.org/0000-0002-6561-2337

Data file description:
This data set contains 5 files of binary data describing the output of the particle-in-cell simulation of magnetic reconnection. The data includes electric and magnetic fields, current density and density. The documentation of the variables and arrays are given below. The dat files are made using Fortran 90. They are named fields-*.dat, for which the * has the usual meaning in a linux environment. The name "fields" refers to the all quantities in Maxwells equations. The number behind the fields, e.g. 00200, refer to the time in units of the inverse of the electron cyclotron frequency. The variables are structured identically in each file, only the time of evaluation is of difference.

Because the dat files contain binary data, they can not be viewed in a text editor. Instead, the variables and arrays can be read into e.g. Python, an open source programming language. The file 01_ReadData.py is the code for reading the fields-*.dat files into Python. 01_ReadData.py can be opened in any text editor.

The simulation setup and units are described in "On the presence and thermalization of cold ions in the exhaust of antiparallel symmetric reconnection' by Norgren et al., 2021.


########################################################################
########################################################################
########################################################################
###################### Documentation of Variables ######################
########################################################################
########################################################################
########################################################################

header # Binary file information
time # Current timestep
dt # Distance from one timestep to another
teti # Ratio of the temperatures for the electron vs protons
xmax # Max value of the x-position of the simulation domain
zmax # Max value of the z-position of the simulation domain
nnx, nnz # nnx and nnz denote the number of grid cells in the x and z direction, respectively.
nns # number of species = 6


########################################################################
###################### Arrays of size (nnx * nnz) ######################
########################################################################

bx # Magnetic field in the x-direction
by # Magnetic field in the y-direction
bz # Magnetic field in the z-direction

ex # Electric field in the x-direction
ey # Electric field in the y-direction
ez # Electric field in the z-direction


########################################################################
############## Arrays of size nnx (xe) and nnz (ze) ####################
########################################################################

xe # Grid layout in the x-direction
ze # Grid layout in the z-direction

wpewce # Ratio of plasma frequency to cyclotron frequency for the electrons


########################################################################
######################## Arrays of size nns ############################
########################################################################

# For python: Index position [0,1,2,3,4,5] corresponds to different species
# 0: Ion species of hot Harris sheet
# 1: Electron of hot Harris sheet
# 2: Cold ions originating from the north
# 3: Cold electrons originating from the north
# 4: Cold ions originating from the south
# 5: Cold electrons originating from the south

mass # Mass for the different species
q # Charge


########################################################################
################# Arrays of size (nss * nnx * nnz) #####################
########################################################################

dns # Information about the density for each species.
dfac # Density factor
pxx # Pressure tensor components
pyy # Pressure tensor components
pzz # Pressure tensor components
pxy # Pressure tensor components
pxz # Pressure tensor components
pyz # Pressure tensor components

vxs # Velocity for the different species in the x-direction
vys # Velocity for the different species in the y-direction
vzs # Velocity for the different species in the z-direction
