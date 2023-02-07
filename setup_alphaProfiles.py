import dustpy
from dustpy import constants as c
import numpy as np

from functions_alphaProfiles import Alpha_Bump
from functions_alphaProfiles import Alpha_DeadZone


################################################################################################
# Helper routine to setup a bumpy profile to create gaps (or bumps) in the gas surface density
################################################################################################

def setup_profile_bumps(sim, Location = 40 * c.au, Amplitude = 4., Width = 1., GasBumpType = 'GAP',
                        apply_to_sigma = True, correct_mass = False, copy_alpha_to_delta = False):
    '''
    Add one or multiple gaussian bumps to the alpha profile to create a GAP or a BUMP in the gas surface density.
    The local pressure maximums in the perturbed gas surface density act as dust traps.
    The bump profile follows the Eq. 5 of Stadler et al.(2022).

    Call the setup function after the initialization and then run, as follows:

    sim.initialize()
    setup_profile_bumps(sim)
    sim.run()
    ----------------------------------------------

    Location [cm]:              Location of the bumps   (int, float, or array)
    Amplitude []:               Amplitude of the bumps  (int, float, or array)
    Width [Hp]:                 FWHM of the bumps in gas scale heights  (int, float, or array)
    BumpType [string]:          'GAP' or 'BUMP' profile outcome in the gas surface density

    apply_to_sigma [Bool]:      if True, add the perturbation into the surface density profile from the beginning of the simulation.
    correct_mass [Bool]:        if True, correct the disk mass to match the sim.ini.gas.Mdisk value (apply_to_sigma must be True)
    copy_alpha_to_delta [Bool]: if True, copy the value of the gas (alpha) turbulence, to the dust (delta) turbulence values

    ----------------------------------------------
    '''

    # Create the description of the gaussian bumps and specify their type
    sim.gas.addgroup("GaussianBumps", description="Gaussian bump parameters, as in Stadler et al. (2022)")
    sim.gas.GaussianBumps.addfield("Location", Location, description = "Location of the bump center (cm)")
    sim.gas.GaussianBumps.addfield("Amplitude", Amplitude, description = "Amplitude of the gaussian bump")
    sim.gas.GaussianBumps.addfield("Width", Width, description = "FWHM of the gaussiam bump (in gas scale heights)")
    sim.gas.GaussianBumps.addfield("Type", GasBumpType, description = "Perturbation type in the gas surface density: GAP or BUMP")

    # The bumps can be updated if the Location, Amplitude, or Width field updaters are asssigned
    # for example with migration, a time-dependent amplitude, or even with bump number
    sim.gas.GaussianBumps.updater = ['Location', 'Amplitude', 'Width']

    # The gas updater is modified.
    # The gaussian bumps and alpha profile are updated after the scale height, since the bump standard deviation is scale height dependant.
    sim.gas.updater = ['gamma', 'mu', 'T', 'cs', 'Hp', 'GaussianBumps', 'alpha', 'nu', 'rho', 'n', 'mfp', 'P', 'eta', 'S']


    # Assign the updated alpha profile
    sim.gas.alpha.updater = Alpha_Bump
    sim.gas.GaussianBumps.update()
    sim.gas.alpha.update()



    if apply_to_sigma:
        rescale_Sigma_with_Alpha(sim, sim.ini.gas.alpha, correct_mass)

    if copy_alpha_to_delta:
        assign_delta_updaters(sim)

    # Update
    sim.update()



################################################################################################
# Helper routine to implement a parametric dead zone model in the gas turbulent alpha parameter
################################################################################################


def setup_profile_deadzone(sim, alpha_active = 1.e-3, alpha_dead = 1.e-4, r_dz_outer = 10 * c.au, width_dz_outer = 1 * c.au,
                        apply_to_sigma = False, correct_mass = True, copy_alpha_to_delta = True):
    '''
    Add a dead zone profile following the parameterization of Garate et al.(2019, 2021).
    This implements a smooth exponential transition between the active and dead regions at the dead zone outer edge.

    Call the setup function after the initialization and then run, as follows:

    sim.initialize()
    setup_profile_deadzone(sim)
    sim.run()
    ----------------------------------------------

    alpha_active:               Alpha turbulence parameter in the MRI active region
    alpha_dead:                 Alpha turbulence parameter in the dead zone region
    r_dz_outer [cm]:            Location of the dead-to-active transition.
    width_dz_outer [cm]:        Width of the transition between dead and active regions

    apply_to_sigma [Bool]:      if True, add the perturbation into the surface density profile from the beginning of the simulation.
    correct_mass [Bool]:        if True, correct the disk mass to match the sim.ini.gas.Mdisk value (apply_to_sigma must be True)
    copy_alpha_to_delta [Bool]: if True, copy the value of the gas (alpha) turbulence, to the dust (delta) turbulence values

    ----------------------------------------------
    '''

    # Create the description of the dead zone parameters
    sim.gas.addgroup("DeadZone", description="Dead Zone parameters, as in Garate et al.(2019, 2021)")
    sim.gas.DeadZone.addfield("alpha_active", alpha_active, description = "Alpha turbulence in the MRI active region")
    sim.gas.DeadZone.addfield("alpha_dead", alpha_dead, description = "Alpha turbulence in the dead zone active region")
    sim.gas.DeadZone.addfield("outer_radii", r_dz_outer, description = "Outer boundary of the dead zone (cm)")
    sim.gas.DeadZone.addfield("transition_width", width_dz_outer, description = "Transition width of the dead zone (cm)")



    # The dead zone parameters can be evolved in time with their corresponding updaters
    sim.gas.DeadZone.updater = ['alpha_active', 'alpha_dead', 'outer_radii', 'transition_width']

    # The gas updater is modified.
    # The dead zone and alpha profile are updated after the scale height, to be consistent with other alpha profile models.
    sim.gas.updater = ['gamma', 'mu', 'T', 'cs', 'Hp', 'DeadZone', 'alpha', 'nu', 'rho', 'n', 'mfp', 'P', 'eta', 'S']


    # Assign the updated alpha profile
    sim.gas.alpha.updater = Alpha_DeadZone
    sim.gas.DeadZone.update()
    sim.gas.alpha.update()


    if apply_to_sigma:
        rescale_Sigma_with_Alpha(sim, alpha_active, correct_mass)


    if copy_alpha_to_delta:
        assign_Delta_updaters(sim)

    # Update
    sim.update()




################################################################################################
# Other helper routines
################################################################################################
def get_DiskMass(sim):
    return np.sum(sim.grid.A * sim.gas.Sigma)

def rescale_Sigma_with_Alpha(sim, alpha_0, correct_mass):
    '''
    When called, the simulation starts with the the surface density scaled in inverse proportion to the alpha profile.
    This causes the simulation to be in quasi-steady state.

    ----------------------------------------------
    alpha_0:        Base turbulence value. The surface density is scaled by alpha_0 / alpha
    correct_mass:   Rescaling the surface density causes the total disk mass to change, the surface density profile is rescaled accordingly.

    ----------------------------------------------
    '''
    original_disk_mass = get_DiskMass(sim)

    sim.gas.Sigma *= alpha_0/sim.gas.alpha
    sim.dust.Sigma *= alpha_0/sim.gas.alpha[:, None]

    scaled_disk_mass = get_DiskMass(sim)

    # Warning - some mass loss/gain occurs depending on the mass removed/added from the bumps.
    mass_ratio = scaled_disk_mass / original_disk_mass
    print('Alpha profile (inversily) applied to the surface density.')
    if correct_mass:
        sim.gas.Sigma /= mass_ratio
        sim.dust.Sigma /= mass_ratio
        print('The surface density profiles were corrected by a scale factor of  {:.3f} to match the initial disk mass.'.format(1./mass_ratio))
    else:
        print('The total disk mass is modified by a factor of  {:.3f}.'.format(mass_ratio))



def get_Alpha(sim):
    return sim.gas.alpha

def assign_Delta_updaters(sim):
    '''
    The dust delta turbulence parameters are assigned the same value as the gas turbulence profile
    '''
    sim.dust.delta.rad.updater = get_Alpha
    sim.dust.delta.turb.updater = get_Alpha
    sim.dust.delta.vert.updater = get_Alpha
