"""
Estimate the probability of Magellanic association for satellite galaxies in the Milky Way system,
using one of the simulations of the Milky Way--LMC encounter
(the choice of simulation is provided as a command-line argument).
"""
import sys, numpy, agama, matplotlib.pyplot as plt

if len(sys.argv) <= 1:
    print('Provide the path to the simulation to examine (e.g., "../L3M11")')
    exit(1)

path = sys.argv[1]

# read in the stored present-day snapshot
data = numpy.load(path + '/snapshot.npz')
# the first nblmc particles originate from the LMC, the next nbhalo represent the Milky Way dark halo,
# and the remaining ones represent its stellar component (not used in the analysis)
nblmc = 2000000
nbhalo = 7000000
thinfactor = 1000  # use only every one in N particles to resimulate
posvel = numpy.hstack([data['pos'], data['vel']])[:nblmc+nbhalo:thinfactor]
tstrip = data['tstrip'][::thinfactor]  # only defined for the first nblmc particles

# time-dependent potential of both galaxies centred on the Milky Way at all times
# (including the acceleration caused by this reference frame being non-inertial)
pot = agama.Potential(path + '/potential.ini')

# integrate the orbits of selected particles back in time in the pre-recorded evolving potential
times = numpy.linspace(-9.0, 0.0, 16*9+1)
orbits = numpy.vstack(
    agama.orbit(ic=posvel, potential=pot, time=times[0], timestart=times[-1], trajsize=len(times), accuracy=1e-4)[:,1]
    ).reshape(len(posvel), len(times), 6)[:,::-1]
plt.ion()
plt.figure(figsize=(12, 8))
plt.axes([0.08, 0.08, 0.9, 0.9])
order = numpy.random.choice(len(posvel), len(posvel), replace=False)   # show particles in random order to avoid bias
for i in range(len(times)):
    plt.cla()
    # colour particles as follows: Milky Way halo by gray, bound LMC by blue, stripped by red
    bound = numpy.zeros(len(posvel), bool)
    bound[:len(tstrip)] = times[i] < tstrip
    color = numpy.zeros((len(posvel), 3))
    color[ bound] = [0, 0, 1]
    color[~bound] = [1, 0, 0]
    color[len(tstrip):] = [0.6, 0.6, 0.6]  # MW particles
    plt.scatter(orbits[order, i, 1], orbits[order, i, 2], c=color[order], s=5, linewidths=0, edgecolors='none')
    plt.xlim(400, -200)
    plt.ylim(-150, 250)
    plt.text(0.98, 0.98, 'time=%4.2f' % times[i], ha='right', va='top', transform=plt.gca().transAxes, fontsize=12)
    plt.draw()
    plt.pause(0.1)

plt.ioff()
plt.show()