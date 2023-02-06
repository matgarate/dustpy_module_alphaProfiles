import numpy as np
import sys
import os


import dustpy
from dustpy import Simulation
from dustpy import constants as c



sim = Simulation()

# Star and disk parameters
sim.ini.star.M = 1.0 * c.M_sun
sim.ini.gas.Mdisk = 0.1 * sim.ini.star.M
sim.ini.gas.SigmaRc = 60 * c.au
sim.ini.gas.alpha = 1.e-3


################################
# GRID SETUP
################################

# Mass Grid Parameters
sim.ini.grid.Nmbpd = 7
sim.ini.grid.mmin = 1.e-12
sim.ini.grid.mmax = 1.e5


# Radial Grid Parameters
sim.ini.grid.Nr = 150
sim.ini.grid.rmin = 2.5 * c.au
sim.ini.grid.rmax = 250 * c.au

sim.initialize()

################################
# ALPHA PROFILE SETUP
################################

from setup_alphaProfiles import setup_profile_bumps, setup_alphaProfiles


# NOTE: Pick only one of the available profiles. Mixing them requires further work.

# Basic setup of a single gap in the gas surface density profile
# By default the bump profile initializes the perturbation in alpha and density profiles
setup_profile_bumps(sim, Location = 40 * c.au, Amplitude = 4., Width = 1., GasBumpType = 'GAP')


# Basic setup of a dead zone profile
# By default it is applied only in the alpha
#setup_profile_deadzone(sim, alpha_active = sim.ini.gas.alpha, alpha_dead = 1.e-4, r_dz_outer = 10*c.au, width_dz_outer = 1 * c.au)




################################
# RUN SIMULATION
################################
print("Running Simulation")

sim.writer.datadir = "Simulation/"
sim.t.snapshots = np.linspace(0.5, 5.0, 10) * 1.e5 * c.year
sim.writer.overwrite = True

sim.run()
