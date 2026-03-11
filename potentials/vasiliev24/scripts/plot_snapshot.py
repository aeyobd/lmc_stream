import sys, numpy, agama, matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt, adaptive_histogram
plt.rc('axes', linewidth=0.5)
plt.rc('font', size=8)
plt.rc('xtick.major', size=1.5)
plt.rc('ytick.major', size=1.5)
numpy.set_printoptions(linewidth=999, precision=4, suppress=True)

if len(sys.argv) <= 1:
    print('Provide the path to the simulation to examine (e.g., "../L3M11")')
    exit(1)

path = sys.argv[1]
outfile = path.replace('/','').replace('.','')
try:
    probs = numpy.loadtxt('satellites_%s.txt' % outfile, dtype=str)
except:
    print('Results file satellites_%s.txt not found, you need to run the script "get_satellites.py" to create it' % outfile)
    exit(1)

# read the stored LMC trajectory
trajlmc = numpy.loadtxt(path + '/trajlmc.txt')
trajtime = trajlmc[:, 0]
poslmc = trajlmc[:, 1:4]
vellmc = trajlmc[:, 4:7]

# read in the stored present-day snapshot
data = numpy.load(path + '/snapshot.npz')
# the first nblmc particles originate from the LMC
nblmc = 2000000
posvel = numpy.hstack([data['pos'], data['vel']])[:nblmc]
mass = data['mass'][:nblmc]
tstrip = data['tstrip']  # only defined for the first nblmc particles
bound = tstrip > 0.0
tstrip[bound] = 0.0

# read in the list of Galactic satellites
satellites = numpy.loadtxt('satellites.txt', dtype=str)
names = satellites[:,1]  # use shortened names
(obs_ra, obs_dec, _, obs_dist, obs_dist_err, obs_vlos, obs_vlos_err,
    obs_pmra, obs_pmra_err, obs_pmdec, obs_pmdec_err, _) = satellites.T[2:].astype(float)
# replace the observed values with the posterior means computed by the classification script
obs_pmra, obs_pmra_err, obs_pmdec, obs_pmdec_err, obs_vlos, obs_vlos_err = probs.T[6:12].astype(float)

# colour the satellites by probability of association with the LMC (gray=low, red=high),
# and switch to blue for the current satellites
problmc = probs[:,4].astype(float)
pcol  = numpy.zeros((len(satellites),3))
pcol[:,0] = numpy.minimum(problmc*2, 1)
for i in range(len(satellites)):
    if names[i] in ['CarII', 'CarIII', 'Del2', 'EriIII', 'HorI', 'HorII', 'HyiI', 'PheII', 'PicII', 'RetII', 'SMC']:
        pcol[i] = [0,0,1]

d2r = numpy.pi/180   # conversion factor from degrees to radians
masyr2kms = 4.74     # conversion factor from mas/yr to km/s (after multiplying by distance in kpc)
obs_l, obs_b, obs_pml, obs_pmb =  agama.transformCelestialCoords(agama.fromICRStoGalactic,
    obs_ra*d2r, obs_dec*d2r, obs_pmra, obs_pmdec)
obs_posvel = numpy.column_stack(agama.getGalactocentricFromGalactic(
    obs_l, obs_b, obs_dist, obs_pml * masyr2kms, obs_pmb * masyr2kms, obs_vlos))

# current position and velocity of the LMC
LMCposvel = numpy.array([-0.61, -41.02, -26.83, -69.84, -221.66, 214.12])
# rotation matrix that aligns the LMC orbital plane with the image plane
ez = numpy.array([
    LMCposvel[1] * LMCposvel[5] - LMCposvel[2] * LMCposvel[4],
    LMCposvel[2] * LMCposvel[3] - LMCposvel[0] * LMCposvel[5],
    LMCposvel[0] * LMCposvel[4] - LMCposvel[1] * LMCposvel[3]
])
LMCangmom = sum(ez**2)**0.5
ez/= LMCangmom
ex = numpy.array([0., -1., 0.])
ex-= ez * ex.dot(ez)
ex/= sum(ex**2)**0.5
ey = numpy.array([ ez[1]*ex[2]-ez[2]*ex[1], ez[2]*ex[0]-ez[0]*ex[2], ez[0]*ex[1]-ez[1]*ex[0] ])

# rotation matrix for the face-on view of Magellanic plane, and corresponding viewing angles
rotmat_f = numpy.vstack([ex, ey, ez])
beta_f   = numpy.arccos (rotmat_f[2,2])
alpha_f  = numpy.arctan2(rotmat_f[2,0],-rotmat_f[2,1])
gamma_f  = numpy.arctan2(rotmat_f[0,2], rotmat_f[1,2])

# another matrix for the edge-on view of Magellanic plane, and corresponding angles
rotmat_e = numpy.vstack([ex, ez, -ey])
beta_e   = numpy.arccos (rotmat_e[2,2])
alpha_e  = numpy.arctan2(rotmat_e[2,0],-rotmat_e[2,1])
gamma_e  = numpy.arctan2(rotmat_e[0,2], rotmat_e[1,2])

pot_lmc_all  = agama.Potential(type='multipole', symmetry='n', lmax=16, rmin=0.5, rmax=1000, gridsizer=50,
    particles=(posvel[:, 0:3] - poslmc[-1], mass), center=poslmc[-1] )
pot_lmc_bound= agama.Potential(type='multipole', symmetry='n', lmax=8,  rmin=0.5, rmax=150,  gridsizer=40,
    particles=(posvel[:, 0:3][bound] - poslmc[-1], mass[bound]), center=poslmc[-1] )

# show surface density of LMC debris in two projections (face-on and edge-on),
# coloured by mean stripping time, together with the velocity vectors of particles and satellites
cmap = plt.get_cmap('circle')
gridx = numpy.linspace(-250, 250, 251)
gridy = numpy.linspace(-150, 250, 201)
grid2d = numpy.column_stack([numpy.repeat(gridx, len(gridy)), numpy.tile(gridy, len(gridx))])

dens_all   = pot_lmc_all  .projectedDensity(grid2d, alpha=alpha_f, beta=beta_f, gamma=gamma_f
   ).reshape(len(gridx), len(gridy))
dens_bound = pot_lmc_bound.projectedDensity(grid2d, alpha=alpha_f, beta=beta_f, gamma=gamma_f
   ).reshape(len(gridx), len(gridy))
bound_frac = numpy.nan_to_num(dens_bound / dens_all)

mod_pos_f = posvel[:,0:3].dot(rotmat_f.T)
mod_vel_f = posvel[:,3:6].dot(rotmat_f.T)
poslmc_f = poslmc.dot(rotmat_f.T)
hist = adaptive_histogram.adaptiveHistogram(mod_pos_f[:, 0:2],
    numpy.column_stack([mass, mass*tstrip, mass*mod_vel_f[:,0], mass*mod_vel_f[:,1]]), gridx, gridy)
dens_max = 2000
dens_min = dens_max * 1e-4
tstrip_min = -8.0
mean_tstrip = hist[:,:,1] / hist[:,:,0]
hist[:,:,0][dens_all < dens_min] = numpy.nan
mean_vx = hist[:,:,2] / hist[:,:,0]
mean_vy = hist[:,:,3] / hist[:,:,0]

plt.figure(figsize=(7.0, 8.4))
ax = plt.axes([0.06, 0.37, 0.93, 0.62])
ax.contour(gridx, gridy, numpy.log10(dens_all.T / dens_max), cmap='Reds',
    levels=numpy.linspace(-4, 0, 17), vmin=-4, vmax=0, linewidths=0.5)
for ll in ax.contour(gridx, gridy, bound_frac.T, levels=[0.5], colors='brown', linewidths=0.5).collections:
    ll.set_dashes([(0, (5,3))])
tstrip_color = cmap(numpy.clip(1 - mean_tstrip.T / tstrip_min, 0, 1))
tstrip_color = 0.9 - (1 - tstrip_color) * 0.4 * numpy.clip(dens_all.T / dens_min - 0.25, 0, 1)[:,:,None]
ax.imshow(tstrip_color[:,:,0:3], extent=(min(gridx), max(gridx), min(gridy), max(gridy)),
    origin='lower', interpolation='bilinear', aspect='auto')

meshx, meshy = numpy.meshgrid(gridx[2::5], gridy[2::5])
vscale = 40.0
ax.quiver(meshx, meshy, mean_vx[2::5,2::5].T, mean_vy[2::5,2::5].T, color='k',
    width=0.001, headwidth=3, headlength=4, scale=vscale, scale_units='x', alpha=0.33)

# plot rainbow-coloured past LMC orbit
for io in range(len(trajtime)-1):
    if trajtime[io+1] < tstrip_min: continue
    ax.plot(poslmc_f[io:io+2,0], poslmc_f[io:io+2,1], color=cmap(1 - trajtime[io+1] / tstrip_min), lw=1.0)
ax.set_xlim(min(gridx), max(gridx))
ax.set_ylim(min(gridy), max(gridy))
ax.set_ylabel('Y [kpc]', labelpad=-6)

# plot satellites
obs_pos_f = obs_posvel[:,0:3].dot(rotmat_f.T)
obs_vel_f = obs_posvel[:,3:6].dot(rotmat_f.T)
arrowparams = dict(edgecolor=None, linewidth=0.5, head_length=1.0, head_width=0.75, overhang=-0.2,
    length_includes_head=True, zorder=10)
fontdict = dict(family='serif', size=5)
for i in range(len(names)):
    plt.arrow(obs_pos_f[i,0], obs_pos_f[i,1], obs_vel_f[i,0]/vscale, obs_vel_f[i,1]/vscale,
        color=pcol[i], **arrowparams)
    plt.text(obs_pos_f[i,0], obs_pos_f[i,1], names[i], zorder=11, fontdict=fontdict, clip_on=True,
        ha='center', va='top' if obs_vel_f[i,1]>0 else 'bottom', color=pcol[i])

lmcp = LMCposvel[0:3].dot(rotmat_f.T)
lmcv = LMCposvel[3:6].dot(rotmat_f.T)
plt.arrow(lmcp[0], lmcp[1], lmcv[0]/vscale, lmcv[1]/vscale, color='b', **arrowparams)
plt.text(lmcp[0], lmcp[1], 'LMC', ha='center', va='top', fontdict=fontdict, zorder=11, color='b')

# MW disc projected onto the LMC orbit plane
ax.add_artist(matplotlib.patches.Ellipse((0,0), 20, 20*abs(ez[2]),
    facecolor=[0,0,0,0.5], edgecolor='k', zorder=4))

gridy = numpy.linspace(-100, 100, 101)
grid2d = numpy.column_stack([numpy.repeat(gridx, len(gridy)), numpy.tile(gridy, len(gridx))])

dens_all   = pot_lmc_all  .projectedDensity(grid2d, alpha=alpha_e, beta=beta_e, gamma=gamma_e
   ).reshape(len(gridx), len(gridy))
dens_bound = pot_lmc_bound.projectedDensity(grid2d, alpha=alpha_e, beta=beta_e, gamma=gamma_e
   ).reshape(len(gridx), len(gridy))
bound_frac = numpy.nan_to_num(dens_bound / dens_all)

mod_pos_e = posvel[:,0:3].dot(rotmat_e.T)
poslmc_e = poslmc.dot(rotmat_e.T)
hist = adaptive_histogram.adaptiveHistogram(mod_pos_e[:, 0:2],
    numpy.column_stack([mass, mass*tstrip]), gridx, gridy)
mean_tstrip = hist[:,:,1] / hist[:,:,0]
hist[:,:,0][dens_all < dens_min] = numpy.nan

ax = plt.axes([0.06, 0.04, 0.93,0.31])
ax.contour(gridx, gridy, numpy.log10(dens_all.T / dens_max), cmap='Reds',
    levels=numpy.linspace(-4, 0, 17), vmin=-4, vmax=0, linewidths=0.5)
for ll in ax.contour(gridx, gridy, bound_frac.T, levels=[0.5], colors='brown', linewidths=0.5).collections:
    ll.set_dashes([(0, (5,3))])
tstrip_color = cmap(numpy.clip(1 - mean_tstrip.T / tstrip_min, 0, 1))
tstrip_color = 0.9 - (1 - tstrip_color) * 0.4 * numpy.clip(dens_all.T / dens_min - 0.25, 0, 1)[:,:,None]
ax.imshow(tstrip_color[:,:,0:3], extent=(min(gridx), max(gridx), min(gridy), max(gridy)),
    origin='lower', interpolation='bilinear', aspect='auto')

# plot rainbow-coloured past LMC orbit
for io in range(len(trajtime)-1):
    if trajtime[io+1] < tstrip_min: continue
    ax.plot(poslmc_e[io:io+2,0], poslmc_e[io:io+2,1], color=cmap(1 - trajtime[io+1] / tstrip_min))
ax.set_xlim(min(gridx), max(gridx))
ax.set_ylim(min(gridy), max(gridy))
ax.set_xlabel('X [kpc]', labelpad=2)
ax.set_ylabel('Z [kpc]', labelpad=-6)

# plot satellites
obs_pos_e = obs_posvel[:,0:3].dot(rotmat_e.T)
obs_vel_e = obs_posvel[:,3:6].dot(rotmat_e.T)
for i in range(len(names)):
    plt.arrow(obs_pos_e[i,0], obs_pos_e[i,1], obs_vel_e[i,0]/vscale, obs_vel_e[i,1]/vscale,
        color=pcol[i], **arrowparams)
    plt.text(obs_pos_e[i,0], obs_pos_e[i,1], names[i], zorder=11, fontdict=fontdict, clip_on=True,
        ha='center', va='top' if obs_vel_e[i,1]>0 else 'bottom', color=pcol[i])

# colorbar for time axis
ax = plt.axes([0.97, 0.07, 0.02, 0.25])
ax.imshow(numpy.linspace(0.125,1,224).reshape(-1,1), extent=(0, 1, tstrip_min+1, 0),
    aspect='auto', interpolation='nearest', origin='lower', cmap=cmap, vmin=0, vmax=1, alpha=0.6)
ax.set_xticks([])
ax.set_ylabel('time [Gyr]')

plt.savefig('fig_debris_%s.pdf' % outfile)
