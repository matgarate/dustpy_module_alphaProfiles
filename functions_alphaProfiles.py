import numpy as np
from dustpy import constants as c
from scipy.interpolate import interp1d


##################################################################################
#
#   Alpha profile to create bump or gap in the gas surface density profile
#
##################################################################################

def Gaussian(r, r0, A, sigma):
    '''
    Standard gaussian shape:
    r:      Radial grid [array, cm]
    r0:     Center of the gaussian [cm]
    A:      Amplitude of the gaussian []
    sigma:  Standard deviation of the gaussian [cm]
    '''

    return A * np.exp(-0.5 * ((r - r0)/sigma)**2)


def get_BumpProfile(sim):
    '''
    Returns a bump profile that follows the following equation:
    BumpProfile = 1 + Sum_i(Bump_i), where Bump_i is a gaussian.
    BumpProfile is an array of size Nr
    '''



    r_bumps = sim.gas.GaussianBumps.Location
    A_bumps = sim.gas.GaussianBumps.Amplitude

    # The interpolation is done in terms of the pressure scale height.
    # The standard deviation is computed assuming that the Hp * width = Gaussian FWHM
    Hp_interpolator = interp1d(sim.grid.r, sim.gas.Hp)
    sigma_bumps = sim.gas.GaussianBumps.Width * Hp_interpolator(r_bumps) / (2. * np.sqrt(2 * np.log(2)))


    bump_profile = np.ones_like(sim.grid.r)


    if type(r_bumps) == int or type(r_bumps) == float:
        bump_profile += Gaussian(sim.grid.r, r_bumps, A_bumps, sigma_bumps)
    else:
        for i in range(len(r_bumps)):
            bump_profile += Gaussian(sim.grid.r, r_bumps[i], A_bumps[i], sigma_bumps[i])


    if sim.gas.GaussianBumps.Type == 'GAP':
        return bump_profile
    elif sim.gas.GaussianBumps.Type == 'BUMP':
        return 1./bump_profile
    else:
        print("PROFILE TYPE NOT DEFINED")
        exit(0)



def Alpha_Bump(sim):
    '''
    Alpha profiled scaled up/down by gaussian bumps
    Feel free to modify the bump profile to the user convenience
    '''

    alpha_0 = sim.ini.gas.alpha
    bump_profile = get_BumpProfile(sim)
    alpha = alpha_0 * bump_profile

    return alpha



##################################################################################
#
#   Alpha profile of a parametric dead zone model
#
##################################################################################

def Alpha_DeadZone(sim):
    '''
    Dead zone profile.
    The outer edge of the dead zone follows a smooth exponential transition as in Garate et al.(2019, 2021)
    '''
    alpha_active = sim.gas.DeadZone.alpha_active
    alpha_dead = sim.gas.DeadZone.alpha_dead
    r_out = sim.gas.DeadZone.outer_radii
    w_out = sim.gas.DeadZone.transition_width


    r = sim.grid.r
    alpha_shape = np.ones_like(r)
    x = r - r_out

    #Define the shape of the dead zone outer boundary
    alpha_shape[x<0] = 0.5 * np.exp(x[x < 0] / w_out)
    alpha_shape[x>=0] = 1.0 - 0.5 * np.exp(-x[x >= 0] / w_out)

    #Rescale the Parameter
    alpha = alpha_dead + (alpha_active - alpha_dead) * alpha_shape

    return alpha
