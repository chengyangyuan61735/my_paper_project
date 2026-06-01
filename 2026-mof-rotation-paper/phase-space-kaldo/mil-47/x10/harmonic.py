# Import necessary packages

from ase.io import read
from kaldo.controllers import plotter
from kaldo.conductivity import Conductivity
from kaldo.forceconstants import ForceConstants
from kaldo.helpers.storage import get_folder_from_label
from kaldo.phonons import Phonons
import matplotlib.pyplot as plt
import numpy as np
import os
import seaborn
plt.style.use('seaborn-v0_8-poster')

from ase.dft.kpoints import get_bandpath



### Set up force constant objects via interface to LAMMPS ####

# Replicate the unit cell 'nrep'=3 times
#nrep = 3
supercell = np.array([2, 2, 4])

# Load in computed 2nd, 3rd IFCs from LAMMPS outputs
forceconstants = ForceConstants.from_folder(folder='fc',supercell=supercell,format='lammps')


# Configure phonon object
# 'k_points': number of k-points
# 'is_classic': specify if the system is classic, True for classical and False for quantum
# 'temperature: temperature (Kelvin) at which simulation is performed
# 'folder': name of folder containing phonon property and thermal conductivity calculations
# 'storage': Format to storage phonon properties ('formatted' for ASCII format data, 'numpy' 
#            for python numpy array and 'memory' for quick calculations, no data stored)

# Define the k-point mesh using 'kpts' parameter
#k_points = 3#'k_points'=3 k points in each direction
phonons_config = {'kpts': [2, 2, 4],
                  'is_classic': True, 
                  'temperature': 300, #'temperature'=300K
                  'folder': 'ALD',
		   'storage': 'formatted'}

# Set up phonon object by passing in configuration details and the forceconstants object computed above
phonons = Phonons(forceconstants=forceconstants, **phonons_config)

# Visualize phonon dispersion, group velocity and density of states with 
# the build-in plotter.

# 'with_velocity': specify whether to plot both group velocity and dispersion relation
# 'is_showing':specify if figure window pops up during simulation

atoms=phonons.atoms
# BandPath
bandpath = get_bandpath(
    path='GX', 
    cell=atoms.cell,
    npoints=20
)

#manually_defined_path = 'G_X'

plotter.plot_dispersion(phonons,with_velocity =False,is_showing=False, manually_defined_path=bandpath)
plotter.plot_dos(phonons,is_showing=False)

# Visualize heat capacity vs frequency and 
# 'order': Index order to reshape array, 
# 'order'='C' for C-like index order; 'F' for Fortran-like index order

# Define the base folder to contain plots
# 'base_folder':name of the base folder
folder = get_folder_from_label(phonons, base_folder='plots')
if not os.path.exists(folder):
        os.makedirs(folder)
# Define a boolean flag to specify if figure window pops during sumuatlion
is_show_fig = False

frequency = phonons.frequency.flatten(order='C')

ps = phonons.phase_space.flatten(order='C')

# 
data = np.column_stack((frequency, ps))

# 
np.savetxt('frequency_phase_space.txt', data, 
           fmt='%.6e',  # 
           delimiter='\t',  # 
           header='frequency\tphase_space')  # 


plt.figure()
plt.scatter(frequency[3:], ps[3:], s=5)
plt.xlabel("$\\nu$ (THz)", fontsize=16)
plt.ylabel("$phase_space_{v} \ (10^{23} \ J/K)$", fontsize=16)
plt.savefig(folder + '/ps_vs_freq.png', dpi=300)
if not is_show_fig:
  plt.close()
else:
  plt.show()
