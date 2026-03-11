#!/usr/bin/python
"""
An example of using the time-dependent potentials of the Milky Way--LMC interaction.
We take the present-day snapshot, select a subset of MW halo particles, and
integrate their orbits back in time for 0.5 Gyr, showing the kinematic maps of
velocity components as observed from the current Solar position.
At the present time, there are prominent asymmetries in these maps,
but they disappear when the perturbation from the LMC is gradually undone.
To run the script, select one of the simulations in this repository
(provide the directory name as the command-line argument).
The plots for spherical_mw_lmc15e10 should be similar to Figure 5 in the paper.
"""
import sys, agama, numpy, matplotlib, matplotlib.pyplot as plt, scipy.special

if len(sys.argv)<=1:
    print('Provide the name of the directory with the simulation snapshot and the potential')
    exit()

def plotMaps(posvel):
    from numpy import sin, cos, arcsin, pi
    lon, lat, dist, pml, pmb, vlos = agama.getGalacticFromGalactocentric(*posvel.T, galcen_v_sun=[0,0,0])
    data = [numpy.ones(len(posvel)), pml * dist, pmb * dist, vlos]
    # create maps for three quantities (vl, vb, vv=vlos) in three radial bins in the Mollweide projection
    sizeX = 120
    sizeY = 60
    gridX = 2 * numpy.linspace(-1 + 1./sizeX, 1 - 1./sizeX, sizeX)
    gridY = 1 * numpy.linspace(-1 + 1./sizeY, 1 - 1./sizeY, sizeY)
    gridXY = numpy.repeat(gridX, sizeY), numpy.tile(gridY, sizeX)
    # inverse Mollweide projection - determine l,b from X,Y
    theta = arcsin(gridXY[1])
    gridb = arcsin( (2*theta + sin(2*theta)) / pi )
    gridl = pi/2 * gridXY[0] / cos(theta)
    usepix= abs(gridl) < pi   # remove pixels outside the ellipse with semiaxes 2x1
    radii = [0, 30, 60, 100]  # radial bins to plot
    # the smooth maps of averaged quantities (three velocity components) are created from the input points
    # by first creating a spherical-harmonic representation of these quantities,
    # and then plotting it on a regular grid in X,Y (which are transformed back to sky coordinates l,b)
    maps = numpy.zeros((len(data), len(radii)-1, len(usepix)), dtype=complex)
    lmax = 6   # order of spherical-harmonic expansion, which determines the angular resolution of maps
    for l in range(lmax+1):
        for m in range(0,l+1):
            Ylm_data = scipy.special.sph_harm(m, l, lon, pi/2-lat)
            Ylm_plot = scipy.special.sph_harm(m, l, -gridl[usepix], pi/2-gridb[usepix])
            for ir in range(len(radii)-1):
                filt = (dist >= radii[ir]) * (dist < radii[ir+1])  # particles in the given radial range
                for iq in range(len(data)):
                    # forward sph-harm transform - compute coefficients from input points for four quantities:
                    # density of points in l,b (unweighted), and density times each velocity component
                    coef = numpy.sum((Ylm_data * data[iq])[filt])
                    # inverse sph-harm transform - evaluate the maps on the grid of points for plotting
                    maps[iq, ir, usepix] += Ylm_plot * coef + (Ylm_plot.conj() * coef.conj() if m!=0 else 0)
    labels = [r'$D\times \mu_l \;\sf[km/s]$', r'$D\times \mu_b \;\sf[km/s]$', r'$v_\mathsf{los} \;\sf[km/s]$']
    for ir in range(len(radii)-1):
        for iq in range(len(data)-1):
            ax = plt.axes([ir*0.30+0.005, 0.65-iq*0.32, 0.28, 0.28])
            ax.set_axis_off()
            # the mean velocity is given by <v*rho> / <rho>
            ax.imshow((maps[iq+1,ir] / maps[0,ir]).real.reshape(sizeX, sizeY).T,
                 cmap='bluered', vmin=-80, vmax=80, extent=[-2,2,-1,1], origin='lower', interpolation='nearest')
            ax.add_artist(matplotlib.patches.Ellipse((0,0), 4, 2, color='k', fill=False, lw=0.5, clip_on=False))
            ax.set_xlim(2,-2)
            ax.set_ylim(-1,1)
            if iq==0:
                ax.text(0, 1.2, '%.0f < D [kpc] < %.0f' % (radii[ir], radii[ir+1]), ha='center', va='center', clip_on=False)
    # add colorbars
    for iq in range(3):
        axc = plt.axes([0.90, 0.65-0.32*iq, 0.02, 0.28])
        axc.imshow(numpy.linspace(0,1,256).reshape(-1,1), extent=[0,1,-80,+80], aspect='auto', interpolation='nearest',
            origin='lower', cmap='bluered', vmin=0, vmax=1)
        axc.set_xticks([])
        axc.yaxis.tick_right()
        axc.set_ylabel(labels[iq], labelpad=-32, fontsize=6)

agama.setUnits(length=1, velocity=1, mass=1)
dirname = sys.argv[1]
# two choices of time-dependent potentials:
# "frozen" - the LMC is moving, but both galaxies retain the initial potentials;
# "evolving" - both galaxies are deforming (in the Milky Way, only the halo, whereas the disk retains the initial profile).
# the latter one is more accurate but also more expensive;
# the main features in the halo kinematics are recovered already in the frozen potential.
if len(sys.argv)>2 and 'evolving' in sys.argv[2].lower():
    pot = agama.Potential(dirname + '/potential_evolving.ini')
else:
    pot = agama.Potential(dirname + '/potential_frozen.ini')
snap = numpy.load(dirname + '/snapshot.npz')
xv   = snap['posvel']
# select a subsample of halo particles
xv   = xv[2000000::10]

plt.rc('axes', linewidth=0.5)
plt.rc('font', size=6)
plt.rc('ytick.major', size=1)
plt.figure(figsize=(4.8,2.4), dpi=250)
plt.ion()
# integrate their orbits backward in time in the time-dependent potentials of the moving LMC and accelerating MW
times= numpy.linspace(0, -0.5, 9)
for k in range(len(times)):
    if k>0:
        xv = numpy.vstack(agama.orbit(potential=pot, ic=xv, time=times[k]-times[k-1], timestart=times[k-1], trajsize=1)[:,1])
    print('plotting kinematic maps at time %.3f Gyr' % times[k])
    plt.clf()
    plotMaps(xv)
    plt.text(1.0, 0.98, 'time: %.2f Gyr ' % times[k], ha='right', va='top', transform=plt.gcf().transFigure)
    plt.draw()
    plt.pause(.1)

plt.ioff()
plt.show()
