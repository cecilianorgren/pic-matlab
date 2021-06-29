Particle-in-cell simulation of magnetotail reconnection with cold plasma protons in inflow regions. Used for manuscript titled 'On the presence and thermalization of cold ions in the exhaust of antiparallel symmetric reconnection' by Norgren et al., 2021

Contact information:
Cecilia Norgren
Space Plasma Physics Group, University of Bergen, 5007 Bergen, Norway.
E-mail: cecilia.norgren@uib.no
ORCID ID of author: https://orcid.org/0000-0002-6561-2337

Data file description:
This data set contains 5 files of binary data describing the output of the particle-in-cell simulation of magnetic reconnection. The data includes electric and magnetic fields, current density and density. The documentation of the variables and arrays are given below. The dat files are made using Fortran 90. They are named fields-*.dat, for which the * has the usual meaning in a linux environment. The name "fields" refers to the all quantities in Maxwells equations. The number behind the fields, e.g. 24000, refer to the time in units of the inverse of the electron cyclotron frequency. The variables are structured identically in each file, only the time of evaluation is of difference.

This dataset contains data files for five timesteps:
fields-18000.dat
fields-20000.dat
fields-22000.dat
fields-24000.dat
fields-25000.dat

Because the dat files contain binary data, they can not be viewed in a text editor. Instead, the variables and arrays can be read into e.g. Matlab (the example files were provided were written and executed in Matlab version R2017b), or an open source programming language like Python. The example loading files we provide are written in the Matlab programming language, but the files can be opened in any text editor, and the code can easily be translated to another language. 

The file example_script.m loads the data into a data structure (described below) using the function function_load_data.m, and then plots some example quantities. The data is accessed using dot-indexing (see example_script.m). Example usage in Matlab:

>> data = function_load_data('/Path/to/your/file/fields-24000.dat');
>> n_hot = data.hot_ion.n;            % density for hot Harris sheet population
>> n_cold = data.cold_ion.n;          % density for cold inflow populations 
>> n_cold_top = data.cold_ion_top.n;  % density for cold inflow population, originating from the north
>> n_cold_bot = data.cold_ion_bot.n;  % density for cold inflow population, originating from the south

###################### Simulation set-up ###############################

The simulation is initalized with a Harris current sheet:

Bx = B0*tanh(z/L) - magnetic field prfile
n  = n0/(cosh(z/L)^2) - density prfile

where L = 2c/ωpi is the thickness of the current sheet, B0 is the symptotic magnetic field value, and n0 is the Harris current sheet peak density.
The ion-to-electron mass ratio is mi/me=100, the ion-to-electron temperature ratio of the hot Harris sheet species is Ti/Te=5, with n0*(Ti+Te)=0.5B0^2, and the electron plasma-to-cyclotron frequency ratio is ωpe/Ωce=2. In addition to the hot Harris sheet populations, we add cold ions and electrons to the northern and southern inflow regions, with a density of nc = 0.2n0. To initiate magnetic reconnection, a magnetic pertubation is added to the center of the box at the start of the simulation.

The simulation setup and units are further described in detail in “On the presence and thermalization of cold ions in the exhaust of antiparallel symmetric reconnection' by Norgren et al., 2021, and adetailed description of the code and post-processing normalization are described below.

###################### Normalization ###################################

Inside the code, we employ a Vlasov scaling, referenced to the electrons.

ωpe - electron plasma frequency 
Ωce - electron cyclotron frequency
vA = B0^2/(mu0*mi*n0)^0.5 - Alfven velocity
mu0 - magnetic permeability in a vacuum
e - electric charge
c - speed of light

Quantity		code		after postprocessing
————————————————————————————————————————————————————————————————————————
Density			n0		n0
Length			c/ωpe		c/ωpi
Time			1/ωpe		1/Ωci
B - magnetic field	ωpe/Ωce*B0	B0
Velocity		c		vA
E - electric field	c*ωpe/Ωce*B0	vA*B0
Fluxes			e*n0*c		e*n0*vA
Charge density		e*n0		e*n0

The transformation between code units and postprocessing units is done by function_load_data.m.


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

