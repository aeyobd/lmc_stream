This repository contains data for the simulations of Large Magellanic Cloud (LMC)'s
interaction with the Milky Way (MW).
Some of the data and all scripts require the Agama framework for stellar dynamics.
Units are: length = 1 kpc, velocity = 1 km/s, time = 0.978 Gyr (for simplicity,
referred to as Gyr throughout the text), mass = 232500 Msun, G = 1.

There are six simulations with different combinations of initial LMC and MW mass
(see Section 2 in Vasiliev 2023 for a definition of initial profiles of both galaxies):

name  M_LMC[Msun]  M_MW[Msun]  time[Gyr]  comment
L2M10       2e11       10e11         11   second passage, Tperi ~ 8.8 Gyr
L2M11       2e11       11e11         10   --"--, Tperi ~ 6.5 Gyr
L3M10       3e11       10e11         11   --"--, Tperi ~ 8.2 Gyr
L3M11       3e11       11e11         10   --"--, Tperi ~ 6.3 Gyr
L2M10first  2e11       10e11          4   first passage starting from apocentre
L3M10rad    3e11       10e11         11   radially anisotropic MW halo (beta=0.5)

Each simulation's folder contains the following data:
- potential.ini:  time-dependent potential of the entire system represented by
  two series of Multipole expansions for each galaxy (the one for the LMC is moving),
  plus the fixed MW stellar potential, plus non-inertial acceleration of the MW-centric
  reference frame. This potential can be used to integrate orbits (resimulation).
- potential_bse.ini:  same, but with potentials represented by BasisSet expansions
  (not recommended, because it is generally slower to evaluate, but this format can
  be more readily used in other galaxy modelling packages).
- potential_lmc_bound.ini:  time-dependent moving potential of particles still bound
  to the Magellanic system (defined so that their total energy in this potential is
  negative).
- potential_lmc_init.ini:  LMC potential at the start of the simulation.
- potential_mw_now.ini:  Milky Way potential at present time (without the non-inertial
  acceleration term).
- potential_mw_init.ini:  Milky Way potential at the start of the simulation.
- trajlmc.txt:  Galactocentric trajectory of the LMC (also used in the potentials).
- boundmass.txt:  evolution of bound mass of the LMC (in units of Msun).
- snapshot.npz:  numpy zip archive with the following arrays:
  * pos:  Nx3 array of particle positions at present time (t=0) in the standard
  Galactocentric coordinate system from Astropy (where the Sun is at -8.122,0,0.02);
  * vel:  Nx3 array of particle velocities;
  * mass:  array of particle masses of length N;
  * tstrip:  time at which the particle became unbound from the LMC (length Nlmc),
  if greater than zero, it is still bound at present time.
  The first Nlmc = 2e6 particles represent the Magellanic system, the next 7e6
  represent the MW halo, and the remaining 1e6 represent the stellar disc+bulge;
  all particles from the same category have identical masses.

The folder "scripts" contains the following Python scripts or Jupyter notebooks:
- find_lmc_orbit.ipynb:  notebook illustrating the method for finding orbital initial
  conditions by running a bunch of simulations with slightly different initial points.
  Instead of running an actual simulation, it creates a fake trajectory by integrating
  equations of motion of LMC and MW as two extended rigid bodies, taking into account
  dynamical friction (but note that this approximation is not particularly realistic).
  One may adapt the code to actual simulations by replacing one specific routine.
- resimulate_orbits.py:  script for performing orbit integration in the recorded
  approximated potential of both galaxies; the actual particle trajecties in the
  simulation may differ individually, but on average this gives a realistic result.
- membership_analysis.py:  script for computing the probability of Magellanic
  association for nearby dwarf galaxies, using data from the text file satellites.txt
  (positions and distances taken from the catalogue of McConnachie,
  https://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/en/community/nearby/
  and proper motions from Battaglia et al.2022).
- plot_snapshot.py:  script for plotting the spatial distribution and velocity map of
  Magellanic debris, together with the observed positions and velocities of dwarfs.
- adaptive_histogram.py:  auxiliary module for variable-bandwidth KDE construction.
