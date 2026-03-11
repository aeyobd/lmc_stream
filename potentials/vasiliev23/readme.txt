This repository contains final snapshots of N-body simulations of the Milky Way and LMC interaction,
accompanying the review paper "The effect of the LMC on the Milky Way system" (Vasiliev 2023).
There are two models for the Milky Way - in both cases the disk and bulge are the same, while
the halo is either spherical (and slightly heavier) or triaxial with a radially changing shape;
the first one is an ad hoc but quite realistic model, and the second one is taken from the "Tango"
paper (Vasiliev et al. 2021, MNRAS, 501, 2279) and provides a good fit for the Sagittarius stream,
but in this repository it is rerun for 5 Gyr into the past instead of 3 Gyr in the original paper.
The LMC is a single-component spherical truncated NFW model with mass 5e10 or 15e10
(the second variant provides a better fit for many features in the Milky Way).
The initial conditions for the LMC orbit are adjusted so that it arrives to the right place
at the right time (i.e., its present-day position and velocity match observations to within 1 kpc
and a couple km/s), but the choice of its current phase-space coordinates differs between the
older Tango simulation and the more recent simplified spherical halo simulation.
The files "snapshot.npz" in each folder contain the present-day snapshots in Galactocentric
coordinates: the Sun is located at -8.2,0,0 and the LMC sits around -0.5,-41,-27 kpc.
"posvel" is the Nx6 array of positions and velocities of particles and "mass" is the array of
particle masses, in the units of 1 kpc, 1 km/s, 1 Msun; the time unit is close to 1 Gyr.
The first 1e6 particles are the LMC, the next 1e6 particles are the stellar disk and bulge of
the Milky Way, and the remaining 4e6 particles are its dark halo.
Other snapshots in the simulation are not provided, but instead the time-dependent potential
of the entire system is represented in the format of the Agama stellar-dynamical framework
(https://agama.software). The total potential is given in the non-inertial reference frame
centered on the Milky Way center, and consists of three components: the Milky Way itself,
the moving LMC (its trajectory in the Galactocentric coordinates is stored in trajlmc.txt),
and the spatially uniform but time-dependent acceleration associated with this non-inertial frame.
There are two versions of the potential in each directory: "frozen" retains the initial potential
of both galaxies (i.e. they are non-deforming, though of course the LMC is moving in space),
and "evolving" contains potentials extracted from the actual N-body snapshots in the last 2 Gyr,
represented by multipole expansions (one for the entire LMC, the other is for the Milky Way halo,
while its disk+bulge are assumed to be fixed). The latter variant is more accurate and tracks
the deformations of both galaxies, but is more expensive when used for orbit integrations.
The example Python script illustrates the usage of these potentials and reproduces one of the
figures from the paper: the kinematic perturbations in the Milky Way halo at present day,
which disappear when the orbits of stars are rewound back in time in the provided time-dependent
potentials, starting from the current phase-space coordinates of particles.
