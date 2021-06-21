Particle-in-cell simulation of magnetotail reconnection with cold plasma protons in inflow regions. Used for manuscript titled 'On the presence and thermalization of cold ions in the exhaust of antiparallel symmetric reconnection' by Norgren et al., 2021

Contact information:
Cecilia Norgren
Space Plasma Physics Group, University of Bergen, 5007 Bergen, Norway.
E-mail: cecilia.norgren@uib.no
ORCID ID of author: https://orcid.org/0000-0002-6561-2337

Data file description:
This data set contains 5 files of binary data describing the output of the particle-in-cell simulation of magnetic reconnection. The data includes electric and magnetic fields, current density and density. The documentation of the variables and arrays are given below. The dat files are made using Fortran 90. They are named fields-*.dat, for which the * has the usual meaning in a linux environment. The name "fields" refers to the all quantities in Maxwells equations. The number behind the fields, e.g. 24000, refer to the time in units of the inverse of the electron cyclotron frequency. The variables are structured identically in each file, only the time of evaluation is of difference.

Because the dat files contain binary data, they can not be viewed in a text editor. Instead, the variables and arrays can be read into e.g. Matlab, or an open source programming language like Python. The example loading files we provide are written in the Matlab programming language, but the files can be opened in any text editor, and the code can easily be translated to another language. 

The file example_script.m loads the data into a data structure (described below) using the function function_load_data.m, and then plots some example quantities. The data is accessed using dot-indexing (see example_script.m).

The simulation setup and units are described in "On the presence and thermalization of cold ions in the exhaust of antiparallel symmetric reconnection' by Norgren et al., 2021.


########################################################################
########################################################################
########################################################################
###################### Documentation of Variables ######################
########################################################################
########################################################################
########################################################################

species  # description of plasma species (6 in total)
time_wpe # timestep in units of inverse electron plasma frequencies
time_wci # timestep in units of inverse ion cyclotron frequencies
nx, nz   # number of grid cells in the x and z directions, respectively
mass     # mass of plasma species
q        # charge of plasma species
wpewce   # Ratio of plasma frequency (wpe) to cyclotron frequency (wce) for the electrons

########################################################################
###################### Arrays of size nx or nz #########################
########################################################################

x_de # x-grid in units of electron inertial lengths
z_de # z-grid in units of electron inertial lengths
x_di # x-grid in units of ion inertial lengths
z_di # z-grid in units of ion inertial lengths

########################################################################
###################### Arrays of size (nx * nz) ########################
########################################################################

Bx # Magnetic field in the x-direction
By # Magnetic field in the y-direction
Bz # Magnetic field in the z-direction

Ex # Electric field in the x-direction
Ey # Electric field in the y-direction
Ez # Electric field in the z-direction

# Plasma moments 
# The plasma moments are provided for all individual species in sub structures:
#   hot_ion      - ion species of hot Harris sheet
#   hot_ele.     - electron of hot Harris sheet
#   cold_ion_top - cold ions originating from the north
#   cold_ele_top - cold electrons originating from the north
#   cold_ion_bot - cold ions originating from the south
#   cold_ele_bot - cold electrons originating from the south
# and for the cold species combined:
#   cold_ion     - cold ions originating from the north and south
#   cold_ele     - cold electrons originating from the north and south

n   # density
jx  # particle flux in x direction
jy  # particle flux in y direction
jz  # particle flux in z direction
vx  # particle speed in x direction, vx = jx/n
vy  # particle speed in y direction, vy = jy/n
vz  # particle speed in z direction, vz = jz/n
pxx # pressure tensor components
pyy # pressure tensor components
pzz # pressure tensor components
pxy # pressure tensor components
pxz # pressure tensor components
pyz # pressure tensor components
p   # scalar pressure p = (pxx+pyy+pzz)/3 
t   # scalar temperature t = p/n

