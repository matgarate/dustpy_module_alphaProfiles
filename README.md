# Alpha Profiles Module for Dustpy

Includes various parametric radial profiles for the alpha turbulence parameter into [DustPy](https://github.com/stammler/dustpy) (Stammler and Birnstiel, [2022](https://ui.adsabs.harvard.edu/abs/2022ApJ...935...35S/abstract)).

The module includes:
* Bump/Gap profile following the implementation of Stadler et al.[(2022)](https://ui.adsabs.harvard.edu/abs/2022A%26A...668A.104S/abstract).
* Dead Zone profile following the implementation of Garate et al.([2019](https://ui.adsabs.harvard.edu/abs/2019ApJ...871...53G/abstract), [2021](https://ui.adsabs.harvard.edu/abs/2021A%26A...655A..18G/abstract))


To setup the alpha profile module add the following lines after initialization

`from setup_alphaProfiles import setup_profile_bumps`

`setup_profile_bumps(sim, Location = 40 * c.au, Amplitude = 4., Width = 1., GasBumpType = 'GAP')`


or



`from setup_alphaProfiles import setup_profile_deadzone`

`setup_profile_deadzone(sim, alpha_active = 1.e-3, alpha_dead = 1.e-4)`


See the `run_alphaProfiles_XX.py` files for an example.

If you use this module, please cite the corresponding papers mentioned above.
