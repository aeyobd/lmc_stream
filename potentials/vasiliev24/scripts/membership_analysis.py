"""
Estimate the probability of Magellanic association for satellite galaxies in the Milky Way system,
using one of the simulations of the Milky Way--LMC encounter
(the choice of simulation is provided as a command-line argument).
"""
import sys, numpy, agama, matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt, matplotlib.backends.backend_pdf, sklearn.mixture
try:
    from ConfigParser import RawConfigParser  # python 2
except ImportError:
    from configparser import RawConfigParser  # python 3
plt.rc('axes', linewidth=0.5)
plt.rc('font', size=6)
plt.rc('xtick.major', size=1.5)
plt.rc('ytick.major', size=1.5)
numpy.random.seed(42)
numpy.set_printoptions(linewidth=9999, suppress=True, precision=5)

if len(sys.argv) <= 1:
    print('Provide the path to the simulation to examine (e.g., "../L3M11")')
    exit(1)

path = sys.argv[1]

# read the stored LMC trajectory
trajlmc = numpy.loadtxt(path + '/trajlmc.txt')
trajtime = trajlmc[:, 0]
poslmc = trajlmc[:, 1:4]
vellmc = trajlmc[:, 4:7]

# read in the stored present-day snapshot
data = numpy.load(path + '/snapshot.npz')
# the first nblmc particles originate from the LMC, the next nbhalo represent the Milky Way dark halo,
# and the remaining ones represent its stellar component (not used in the analysis)
nblmc = 2000000
nbhalo = 7000000
posvel = numpy.hstack([data['pos'], data['vel']])[:nblmc+nbhalo].astype(float)
mass = data['mass'][:nblmc+nbhalo].astype(float)
tstrip = data['tstrip']  # only defined for the first nblmc particles

# time-dependent potential of both galaxies centred on the Milky Way at all times
# (including the acceleration caused by this reference frame being non-inertial)
pot_all = agama.Potential(path + '/potential.ini')
# potential of the bound particles in the LMC (moving around the Milky Way centre)
pot_lmc = agama.Potential(path + '/potential_lmc_bound.ini')
# initial potential of the Milky Way (centred on itself)
pot_mw_init = agama.Potential(path + '/potential_mw_init.ini')
# initial potential of the LMC (also centred on itself)
pot_lmc_init = agama.Potential(path + '/potential_lmc_init.ini')

# create initial distribution function model of the Milky Way halo
iniFile = RawConfigParser()
iniFile.read(path + '/potential_mw_init.ini')
den_mw_halo = agama.Density(dict(iniFile.items('Potential halo')))
df_mw_halo = agama.DistributionFunction(density=den_mw_halo, potential=pot_mw_init, **dict(iniFile.items('DF halo')))
gm_mw_halo = agama.GalaxyModel(pot_mw_init, df_mw_halo)
# same for the LMC (a single-component spherical model is somewhat simpler to set up)
df_lmc = agama.DistributionFunction(type='QuasiSpherical', potential=pot_lmc_init)
gm_lmc = agama.GalaxyModel(pot_lmc_init, df_lmc)

# read in the list of Galactic satellites
satellites = numpy.loadtxt('satellites.txt', dtype=str)
names = satellites[:,0]
(obs_ra, obs_dec, _, obs_dist, obs_dist_err, obs_vlos, obs_vlos_err,
    obs_pmra, obs_pmra_err, obs_pmdec, obs_pmdec_err, _) = satellites.T[2:].astype(float)

d2r = numpy.pi/180   # conversion factor from degrees to radians
masyr2kms = 4.74     # conversion factor from mas/yr to km/s (after multiplying by distance in kpc)
obs_coord_err = 4.0  # assumed "uncertainty" on celestial coordinates [degrees]
# replace missing values in PM and Vlos by zeroes with very large uncertainty
novl = numpy.isnan(obs_vlos_err)
nopm = numpy.isnan(obs_pmra_err)
obs_vlos     [novl] = 0
obs_vlos_err [novl] = 1000
obs_pmra     [nopm] = 0
obs_pmdec    [nopm] = 0
obs_pmra_err [nopm] = 1000 / obs_dist[nopm] / masyr2kms
obs_pmdec_err[nopm] = 1000 / obs_dist[nopm] / masyr2kms

# impose a lower limit on the velocity uncertainty and relative distance uncertainty
min_vel_err  = 10.0  # km/s
min_pm_err   = 0.05  # mas/yr
min_dist_err = 0.1   # fractional
obs_vlos_err = numpy.maximum(numpy.nan_to_num(obs_vlos_err), min_vel_err)
obs_pmra_err = numpy.maximum(numpy.nan_to_num(obs_pmra_err), numpy.maximum(min_pm_err, min_vel_err/obs_dist/masyr2kms))
obs_pmdec_err= numpy.maximum(numpy.nan_to_num(obs_pmdec_err),numpy.maximum(min_pm_err, min_vel_err/obs_dist/masyr2kms))
obs_dist_err = numpy.maximum(numpy.nan_to_num(obs_dist_err), min_dist_err*obs_dist)

mod_l, mod_b, mod_dist, mod_pml, mod_pmb, mod_vlos = agama.getGalacticFromGalactocentric(*posvel.T)
mod_pml /= masyr2kms  # convert from km/s/kpc to mas/yr
mod_pmb /= masyr2kms

# distance between two points on celestial sphere
def angular_distance(ra0, dec0, ra1, dec1):
    return 2 * numpy.arcsin( (numpy.sin( (dec0-dec1)*0.5 * d2r )**2 +
        numpy.cos(dec0 * d2r) * numpy.cos(dec1 * d2r) * numpy.sin( (ra0-ra1)*0.5 * d2r )**2 )**0.5 ) / d2r


def classify(obs_ind):
    """
    # Determine the probability of being [formerly or currently] associated with the LMC for a given satellite galaxy
    """
    # rotate the celestial coordinates to align the object with (lon,lat)=(0,0):
    # this does not affect pmra,pmdec for points very close to the object coordinates,
    # but straightens out the distortions of more distant points caused by the curvilinear
    # celestial coordinate system if the object is near north/south pole.
    rotmat  = agama.makeCelestialRotationMatrix(obs_ra[obs_ind]*d2r, obs_dec[obs_ind]*d2r, 0)
    obs_lon, obs_lat, obs_pmlon, obs_pmlat, obs_pmlon_err, obs_pmlat_err = agama.transformCelestialCoords(rotmat,
        obs_ra[obs_ind]*d2r, obs_dec[obs_ind]*d2r, obs_pmra[obs_ind],
        obs_pmdec[obs_ind], obs_pmra_err[obs_ind], obs_pmdec_err[obs_ind], 0)[0:6]
    obs_val = numpy.array([obs_lon, obs_lat, obs_dist[obs_ind], obs_pmlon, obs_pmlat, obs_vlos[obs_ind]])
    obs_err = numpy.array([1e-3, 1e-3, 10*obs_dist_err[obs_ind], obs_pmlon_err, obs_pmlat_err, obs_vlos_err[obs_ind]])

    # particles in the snapshot also rotated into the same celestial frame
    mod_lon, mod_lat, mod_pmlon, mod_pmlat = agama.transformCelestialCoords(rotmat.dot(agama.fromGalactictoICRS),
        mod_l, mod_b, mod_pml, mod_pmb)
    mod_val = numpy.column_stack([mod_lon, mod_lat, mod_dist, mod_pmlon, mod_pmlat, mod_vlos])

    # probability of spatial matching (incorporating the distance uncertainty)
    angdist = angular_distance(obs_lon/d2r, obs_lat/d2r, mod_lon/d2r, mod_lat/d2r)
    difpos2 = (angdist / obs_coord_err)**2 + ((obs_dist[obs_ind] - mod_dist) / obs_dist_err[obs_ind])**2
    probpos = numpy.exp(-0.5 * difpos2)

    '''
    Method 1: find the particles in the snapshot that are close to the observed galaxy in all 6d phase-space coordinates,
    normalized by observational errors (Mahalanobis distance);
    if the number of nearby particles is too low, keep increasing the errors in velocity until get enough of them.
    '''
    errmul = 1.0
    while True:
        probkin = numpy.exp(-0.5 * numpy.sum((obs_val - mod_val)[:,3:6]**2 / (obs_err[3:6] * errmul)**2, axis=1) )
        probtot = mass * probpos * probkin
        probsum = numpy.sum(probtot)
        if numpy.max(probtot) < 0.1 * probsum: break
        errmul *= 1.5
    problmc_part = numpy.sum(probtot[:nblmc]) / probsum
    # obtain a probability-weighted sample of matching particles
    nsamples = 50
    samples_part_indices = numpy.searchsorted(numpy.cumsum(probtot) / probsum,
        numpy.linspace(0.5/nsamples, 1-0.5/nsamples, nsamples))
    # randomize the indices to avoid perceptual biases in subsequent plotting
    numpy.random.shuffle(samples_part_indices)

    '''
    Method 2: Gaussian Mixture Model in the 6d phase space,
    but with the three spatial dimensions pre-filtered by spatial matching probability
    (so effectively a GMM for the kinematic variables in the given spatial region)
    '''
    ndim  = 6  # all dimensions of phase space
    ncat  = 2  # two alternative categories in the mixture model (LMC=0, MW=1)
    probgmm = numpy.zeros(ncat)  # un-normalized probabilities of matching for each category
    posterior_moment1 = numpy.zeros(ndim)  # accumulators for computing the posterior mean and error estimate
    posterior_moment2 = numpy.zeros(ndim)  # (only keep track of diagonal elements of posterior covariance matrix)
    gmm_param = numpy.zeros(ncat, object)  # parameters of fitted GMM for each category
    for cat in range(ncat):
        filt = numpy.hstack([numpy.repeat(1-cat, nblmc), numpy.repeat(cat, nbhalo)]).astype(bool)
        csum = numpy.cumsum(probpos * mass * filt)
        catmass = csum[-1]  # total mass of all particles of this category weighted by spatial probability
        meanmass = numpy.mean(mass[filt])
        if catmass < 5 * meanmass:   # tooooo few particles of a given subset, bail out
            gmm_param[cat] = dict(norm=0, ampl=numpy.zeros(1), mean=numpy.zeros((1,ndim)), covar=numpy.ones((1,ndim,ndim)))
            continue
        ncomp = 1 if catmass < 50*meanmass else 2   # number of GMM components for this category

        # we construct a GMM for all particles in the snapshot, but with their weights already pre-multiplied
        # by the spatial matching probability. Since the GMM operates on equal-weight samples,
        # we construct such a sample by repeating each particle nmul*probpos times
        # (for the vast majority of them, probpos is so small that they will not make it into this sample even once)
        nmul = 10
        nsamples = round(catmass / meanmass * nmul)
        inds = numpy.searchsorted(csum / catmass, numpy.linspace(0.5/nsamples, 1-0.5/nsamples, nsamples))
        gmm  = sklearn.mixture.GaussianMixture(ncomp, 'full', tol=1e-6, max_iter=1000, n_init=1).fit(mod_val[inds])
        ampl = gmm.weights_ / sum(gmm.weights_)  # normalized amplitudes (summing up to unity)
        mod_mean = gmm.means_.copy()
        mod_covar = gmm.covariances_.copy()
        gmm_param[cat] = dict(norm=catmass, ampl=ampl, mean=mod_mean, covar=mod_covar)

        # to compute the probability of association for the given observed object, we need to convolve the model
        # distribution with the Gaussian measurement uncertainty, and then evaluate it at the observed values.
        # For a GMM, this is trivial - just add the covariance matrices of the model and the observed uncertanties;
        # we do not need to add the positional uncertainty, as it was already incorporated into the GMM construction.
        obs_covar = numpy.diag(obs_err**2)
        inv_covar = numpy.linalg.inv(mod_covar + obs_covar)
        inv_mod_covar = numpy.linalg.inv(mod_covar)  # weight matrices for computing the posterior mean and variance
        inv_obs_covar = numpy.linalg.inv(obs_covar)
        for comp in range(ncomp):
            probcomp = catmass * numpy.linalg.det(inv_covar[comp])**0.5 * numpy.exp(
                -0.5 * inv_covar[comp].dot(obs_val - mod_mean[comp]).dot(obs_val - mod_mean[comp]) )
            probgmm[cat] += probcomp
            posterior_covar = numpy.linalg.inv(inv_mod_covar[comp] + inv_obs_covar)
            posterior_mean  = posterior_covar.dot(inv_mod_covar[comp].dot(mod_mean[comp]) + inv_obs_covar.dot(obs_val) )
            posterior_moment1 += probcomp * posterior_mean
            posterior_moment2 += probcomp * (posterior_mean**2 + numpy.diag(posterior_covar))
    problmc_gmm = probgmm[0] / sum(probgmm)
    # posterior mean values and uncertainties of the object coordinates
    posterior_val = posterior_moment1 / sum(probgmm)
    posterior_err = numpy.maximum(0, posterior_moment2 / sum(probgmm) - posterior_val**2)**0.5

    '''
    Method 3: integrate orbits uniformly sampled from measurement uncertainties back in time,
    record the phase-space coordinates at the beginning of evolution,
    then reweigh them according to the prior probability of finding such coords in each of the two galaxies.
    '''
    # initial conditions for orbits are unweighted samples from measurement uncertainties,
    # transformed into Galactocentric coordinates
    num_ic_raw  = 100000  # initial number of random samples
    num_ic_filt = 5000    # max number of computed orbits after pre-filtering the initial sample
    rand = numpy.random.normal(size=(ndim, num_ic_raw))
    samp_l, samp_b, samp_pml, samp_pmb = agama.transformCelestialCoords(
        agama.fromICRStoGalactic.dot(rotmat.T),
        rand[0] * obs_coord_err * d2r, rand[1] * obs_coord_err * d2r,
        obs_pmlon + rand[2] * obs_pmlon_err, obs_pmlat + rand[3] * obs_pmlat_err)
    ic = numpy.column_stack(agama.getGalactocentricFromGalactic(
        samp_l, samp_b, obs_dist[obs_ind] + rand[4] * obs_dist_err[obs_ind],
        samp_pml * masyr2kms, samp_pmb * masyr2kms, obs_vlos[obs_ind] + rand[5] * obs_vlos_err[obs_ind]))
    # pre-filter the sample, removing values of velocity that fall outside of the range spanned by actual particles
    # in the simulation, as these points are unlikely to be bound to either galaxy, and limit the remaining sample
    mod_range_min, mod_range_max = numpy.array(numpy.percentile(posvel[probpos > 1e-4], [0,100], axis=0))
    ic = ic[numpy.all((ic >= mod_range_min) * (ic <= mod_range_max), axis=1)][:num_ic_filt]

    # integrate the trajectories back in time in the pre-recorded evolving potential
    orbittime = numpy.linspace(trajtime[0], trajtime[-1], int(round((trajtime[-1] - trajtime[0]) / 0.0625)) + 1)
    orbits_rewind = numpy.vstack(
        agama.orbit(ic=ic, potential=pot_all, time=trajtime[0], timestart=trajtime[-1], trajsize=len(orbittime))[:,1]
        ).reshape(len(ic), len(orbittime), 6)[:,::-1]

    # phase-space coordinates at the beginning of evolution relative to each galaxy
    posvel_begin_mw  = orbits_rewind[:,0]
    posvel_begin_lmc = orbits_rewind[:,0] - trajlmc[0, 1:7]

    # compute [unnormalized] prior probabilities of each trajectory for two alternative origins (LMC or MW),
    # by taking the value of the corresponding DF at the actions in the respective galaxy's initial potential,
    # and replacing invalid values (when the initial point is unbound) by zeroes
    prob_mw  = numpy.nan_to_num(df_mw_halo(gm_mw_halo.af(posvel_begin_mw )))
    prob_lmc = numpy.nan_to_num(df_lmc    (gm_lmc    .af(posvel_begin_lmc)))
    prob_orb = prob_mw + prob_lmc + 1e-100  # overall posterior probability of each orbit [still unnormalized]
    problmc_rewind = numpy.sum(prob_lmc) / numpy.sum(prob_orb)  # overall normalized odds ratio for this object

    # create a subset of orbits from this sample, choosing them in proportion to their posterior probability,
    # and assign a category to each orbit based on the ratio of posterior probabilities for this orbit
    nsamples = 50
    samples_rewind_indices = numpy.searchsorted(numpy.cumsum(prob_orb) / numpy.sum(prob_orb),
        numpy.linspace(0.5/nsamples, 1-0.5/nsamples, nsamples))
    samples_rewind_cat = numpy.random.random(size=nsamples) < (prob_mw / prob_orb)[samples_rewind_indices]

    '''
    Done with computing, now create plots
    '''
    fig = plt.figure(figsize=(3.3*3/2, 4.4))

    # display only particles with non-negligible probability of spatial matching
    filt = numpy.where(probpos > 1e-3)[0]
    # randomize their order to avoid perception bias in case they start to overlap
    numpy.random.shuffle(filt)
    # assign particle colour based on its category, and transparency - based on probability of spatial matching
    colors = numpy.array([[1,0,0], [0.5,0.5,0.5]])  # red=LMC, gray=MW
    color  = numpy.vstack([numpy.tile(colors[0], nblmc).reshape(nblmc,3), numpy.tile(colors[1], nbhalo).reshape(nbhalo,3)])
    pcolor = numpy.column_stack([color, probpos])[filt]
    # four dimensions to show on the corner plot (distance and three kinematic dimensions)
    labels = ['pmra', 'pmdec', 'vlos', 'distance']
    dims   = (3,4,5,2)
    # plotting range should enclose the vast majority of particles
    limits = [ numpy.percentile(mod_val[filt, k], [1,99]) for k in dims ]
    # and also make sure to enclose the observed values (if available),
    # or replace missing values with the posterior means
    if not nopm[obs_ind]:
        limits[0][0] = min(limits[0][0], obs_pmlon - 1.2*obs_pmlon_err)
        limits[0][1] = max(limits[0][1], obs_pmlon + 1.2*obs_pmlon_err)
        limits[1][0] = min(limits[1][0], obs_pmlat - 1.2*obs_pmlat_err)
        limits[1][1] = max(limits[1][1], obs_pmlat + 1.2*obs_pmlat_err)
    else:
        obs_val[dims[0]] = posterior_val[dims[0]]
        obs_val[dims[1]] = posterior_val[dims[1]]
    if not novl[obs_ind]:
        limits[2][0] = min(limits[2][0], obs_vlos[obs_ind] - 1.2*obs_vlos_err[obs_ind])
        limits[2][1] = max(limits[2][1], obs_vlos[obs_ind] + 1.2*obs_vlos_err[obs_ind])
    else:
        obs_val[dims[2]] = posterior_val[dims[2]]
    obs_err[2] = obs_dist_err[obs_ind]  # for display purposes, put back the distance uncertainty

    # corner plot - 2d projections of particle distributions
    for i in range(1,4):
        for j in range(i):
            ax = fig.add_axes([(0.12 + j*0.296)*2/3, 0.99 - i*0.222, 0.28*2/3, 0.21])
            # all particles satisfying the spatial selection, which are considered in the GMM of Method 2
            ax.scatter(mod_val[:,dims[j]][filt], mod_val[:,dims[i]][filt],
                c=pcolor, s=1, rasterized=True, linewidths=0, edgecolors='none')
            # particles selected to match the object in both spatial and kinematic dimensions (Method 1)
            ax.scatter(mod_val[:,dims[j]][samples_part_indices], mod_val[:,dims[i]][samples_part_indices],
                c=colors[(samples_part_indices>=nblmc).astype(int),], s=2, marker='x', linewidths=0.5)
            # measured values and uncertainties of the object
            ax.errorbar(obs_val[dims[j]], obs_val[dims[i]], xerr=obs_err[dims[j]], yerr=obs_err[dims[i]],
                color='b', marker=None, linewidth=1, capsize=1.5)
            ax.set_xticklabels([])
            if j==0:
                plt.yticks(rotation=45)
                ax.tick_params(pad=2)
                ax.text(-0.35, 0.5, labels[i], rotation=90, transform=ax.transAxes, fontsize=7, ha='center', va='center')
            else:
                ax.set_yticklabels([])
            ax.set_xlim(limits[j])
            ax.set_ylim(limits[i])

    # bottom row - marginalized 1d projections
    for j in range(3):
        d = dims[j]
        ax = fig.add_axes([(0.12 + j*0.296)*2/3, 0.99 - 4*0.222, 0.28*2/3, 0.21])
        plt.xticks(rotation=45)
        ax.tick_params(pad=2)
        ax.text(0.5, -0.40, labels[j], transform=ax.transAxes, fontsize=7, ha='center', va='center')
        ax.set_yticks([])
        xgrid = numpy.linspace(limits[j][0], limits[j][1], 100)
        ymax = 0.
        for cat in range(ncat):
            yvals = xgrid*0
            for c in range(len(gmm_param[cat]['ampl'])):
                yvals +=  (gmm_param[cat]['norm'] * gmm_param[cat]['ampl'][c] * gmm_param[cat]['covar'][c,d,d]**-0.5 *
                    numpy.exp(-0.5 * (xgrid - gmm_param[cat]['mean'][c,d])**2 / gmm_param[cat]['covar'][c,d,d]) )
            ymax = max(ymax, max(yvals))
            ax.plot(xgrid, yvals, color=colors[cat])
        ax.set_ylim(0, ymax * 1.1)
        ax.set_xlim(limits[j])
        ax.set_yticks([])
        ax.fill_between([obs_val[j+3] - obs_err[j+3], obs_val[j+3] + obs_err[j+3]], [0, 0], [1, 1],
            color='#c0c0ff', zorder=-3, lw=0, transform=ax.get_xaxis_transform())

    # side panels: possible past orbits of the object according to methods 1 and 3.
    def show_orbit(orbit, tstrip):
        dist_mw  = numpy.sum( orbit[:,0:3]**2, axis=1)**0.5
        dist_lmc = numpy.sum((orbit[:,0:3] - poslmc[numpy.searchsorted(trajtime, orbittime)])**2, axis=1)**0.5
        params = dict(lw=0.5, alpha=0.5, zorder=-1)
        if tstrip is None:  # native MW orbit
            ax_lmc.plot(orbittime, dist_lmc, color='grey', **params)
            ax_mw .plot(orbittime, dist_mw , color='grey', **params)
        elif tstrip > 0:
            ax_lmc.plot(orbittime, dist_lmc, color='b', **params)
            ax_mw .plot(orbittime, dist_mw , color='b', **params)
        else:
            istrip = numpy.searchsorted(orbittime, tstrip)
            ax_lmc.plot(orbittime[:istrip ], dist_lmc[:istrip ], color='b', **params)
            ax_mw .plot(orbittime[:istrip ], dist_mw [:istrip ], color='b', **params)
            ax_lmc.plot(orbittime[ istrip:], dist_lmc[ istrip:], color='r', **params)
            ax_mw .plot(orbittime[ istrip:], dist_mw [ istrip:], color='r', **params)
            ax_lmc.plot(orbittime[ istrip-1:istrip+1], dist_lmc[istrip-1:istrip+1], color='m', **params)
            ax_mw .plot(orbittime[ istrip-1:istrip+1], dist_mw [istrip-1:istrip+1], color='m', **params)

    def setup_orbit_axes(label):
        ax_mw .plot(trajtime, numpy.sum(poslmc**2, axis=1)**0.5, color='m', dashes=[2,2], lw=1)
        ax_mw .set_xticklabels([])
        ax_lmc.set_xlabel('time', labelpad=2, fontsize=7)
        ax_mw .set_xlim(trajtime[0], trajtime[-1])
        ax_lmc.set_xlim(trajtime[0], trajtime[-1])
        ax_mw .set_ylim(0, 400)
        ax_lmc.set_ylim(0, 400)
        ax_mw .set_ylabel('distance to MW',  labelpad=2, fontsize=7)
        ax_lmc.set_ylabel('distance to LMC', labelpad=2, fontsize=7)
        ax_lmc.text(1.08, 1.02, label, ha='center', va='center', rotation=90, fontsize=7, transform=ax_lmc.transAxes)

    # Method 1: selected particles from the actual simulation, but rewound in the smooth potential approximation
    orbits_part = numpy.vstack(
        agama.orbit(ic=posvel[samples_part_indices], potential=pot_all,
        time=trajtime[0], timestart=trajtime[-1], trajsize=len(orbittime))[:,1]
        ).reshape(len(samples_part_indices), len(orbittime), 6)[:,::-1]
    ax_mw  = fig.add_axes([0.76, 0.780, 0.21, 0.21])
    ax_lmc = fig.add_axes([0.76, 0.555, 0.21, 0.21])
    setup_orbit_axes('Method 1 (particles)')
    for i in range(len(samples_part_indices)):
        show_orbit(orbits_part[i], tstrip[samples_part_indices[i]] if samples_part_indices[i] < nblmc else None)

    # Method 3: orbits sampled from the posterior distribution of measured phase-space coordinates
    orbits_rewind = orbits_rewind[samples_rewind_indices]
    ax_mw  = fig.add_axes([0.76, 0.285, 0.21, 0.21])
    ax_lmc = fig.add_axes([0.76, 0.060, 0.21, 0.21])
    setup_orbit_axes('Method 3 (rewinding)')
    for i in range(len(samples_rewind_indices)):
        if samples_rewind_cat[i]:
            tstrip_i = None  # native MW orbit
        else:
            # determine the stripping time of the orbit by finding when it first acquired positive energy w.r.t. LMC
            unbound = numpy.where(
                pot_lmc.potential(orbits_rewind[i,:,0:3], t=orbittime) +
                0.5 * numpy.sum( (orbits_rewind[i,:,3:6] - vellmc[numpy.searchsorted(trajtime, orbittime)])**2, axis=1)
                > 0)[0]
            tstrip_i = orbittime[unbound[0]] if len(unbound)>0 else numpy.inf
        show_orbit(orbits_rewind[i], tstrip_i)

    # various annotations on the entire figure
    ax = fig.add_axes([0,0,2./3,1])
    ax.set_axis_off()
    ax.set_xlim(0,1)
    ax.set_ylim(0,1)
    ax.errorbar(0.74, 0.62, xerr=0.01, yerr=0.01, marker=None, capsize=1.5, lw=1, color='b')
    ax.text(0.45, 0.98, names[obs_ind], ha='left', va='top', fontsize=10)
    ax.text(0.45, 0.91, 'probability of Magellanic association', fontsize=7)
    ax.text(0.45, 0.87, 'particles',    fontsize=7)
    ax.text(0.45, 0.83, 'velocity GMM', fontsize=7)
    ax.text(0.45, 0.79, 'orbit rewind', fontsize=7)
    ax.text(0.72, 0.87, '%.2f' % problmc_part,   fontsize=7)
    ax.text(0.72, 0.83, '%.2f' % problmc_gmm,    fontsize=7)
    ax.text(0.72, 0.79, '%.2f' % problmc_rewind, fontsize=7)
    ax.text(0.72, 0.74, 'N-body particles:', fontsize=7)
    ax.text(0.72, 0.70, 'Milky Way', color='gray', fontsize=7)
    ax.text(0.72, 0.66, 'Magellanic system', color='r', fontsize=7)
    ax.text(0.77, 0.615,'observations', color='b', fontsize=7)
    pdf.savefig()
    #fig.savefig(names[obs_ind]+'.pdf', dpi=150)
    plt.close()

    print('%-15s:  Plmc,part=%.3f,  Plmc,gmm=%.3f,  Plmc,rew=%.3f' % (names[obs_ind], problmc_part, problmc_gmm, problmc_rewind))
    return (problmc_part, problmc_gmm, problmc_rewind), posterior_val[3:6], posterior_err[3:6]

outfile = path.replace('/','').replace('.','')
with matplotlib.backends.backend_pdf.PdfPages('fig_satellites_' + outfile + '.pdf') as pdf, \
    open('satellites_' + outfile + '.txt', 'w') as out:

    # output file contains the estimates of probability of Magellanic association according to three methods,
    # their mean value and scatter, and posterior values of PM and Vlos (measured values reweighted by GMM distributions)
    out.write('# name         \tparticl\tGMM\trewind\tmean\tscatter\tpmra\te_pmra\tpmdec\te_pmdec\tvlos\te_vlos\n')
    for obs_ind in range(len(satellites)):
        prob, post_val, post_err = classify(obs_ind)
        out.write(('%-15s' + '\t%7.3f'*9 + '%7.1f'*2 + '\n') %
            (names[obs_ind], prob[0], prob[1], prob[2], numpy.mean(prob), numpy.std(prob),
            post_val[0], post_err[0], post_val[1], post_err[1], post_val[2], post_err[2]))
