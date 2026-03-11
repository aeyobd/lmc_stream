"""
Construction of N-dimensional kernel density estimators with variable bandwidth
proportional to the mean local interparticle distance.
Run the script without parameters to get an illustration.
"""
import numpy, scipy.spatial, matplotlib.pyplot as plt

def simpleHistogram(point_coords, quantities, gridx, gridy):
    """
    Compute a simple histogram of quantities on a 2d grid from points
    Arguments:
      point_coords: coordinates and masses of points (shape: Nx2);
      quantities:   array of NxQ quantities to sum up (e.g., point mass, mass times velocity, etc.);
      gridx, gridy: coords of the 2d grid at which to compute the sums (centers of bins/pixels);
    Return:         3d array of shape (len(gridx), len(gridy), Q), where each grid pixel contains
                    the sum of quantities of points that lie within this pixel divided by the pixel area
    """
    quantities = numpy.atleast_2d(quantities.T).T  # if the input was a 1d array, make it a Nx1 2d array
    # gridx, gridy are the coordinates of grid cell centers, need to convert them to cell boundaries
    binx = numpy.hstack([ gridx[0]*1.5 - gridx[1]*0.5, (gridx[1:] + gridx[:-1])*0.5, gridx[-1]*1.5 - gridx[-2]*0.5 ])
    biny = numpy.hstack([ gridy[0]*1.5 - gridy[1]*0.5, (gridy[1:] + gridy[:-1])*0.5, gridy[-1]*1.5 - gridy[-2]*0.5 ])
    return numpy.dstack([
        numpy.histogram2d(point_coords[:,0], point_coords[:,1], bins=(binx, biny), weights=quantities[:,i])[0]
        for i in range(quantities.shape[1])
    ]) / (binx[1:] - binx[:-1])[:,None,None] / (biny[1:] - biny[:-1])[None,:,None]

def adaptiveHistogram(point_coords, quantities, gridx, gridy, K=200, maxdist=numpy.inf, nthreads=-1):
    """
    Compute a histogram-like statistics for some quantities on a 2d grid from points,
    using a kernel density estimate from K nearest neighbours (i.e. with a variable kernel size).
    Arguments:
      point_coords: coordinates and masses of points (shape: Nx2);
      quantities:   array of NxQ quantities to sum up (e.g., point mass, mass times velocity, etc.);
      gridx, gridy: coords of the 2d grid at which to compute the sums (centers of bins/pixels);
      K:            number of neighbour points used in density estimate for each grid node;
      maxdist:      if not infinite, suppress reporting the results for grid nodes (pixels)
                    that are more than maxdist away from the nearest point (replace by NaN);
      nthreads:     number of CPU threads (-1 means use all available cores), available on scipy>=0.16.
    Return:         3d array of shape (len(gridx), len(gridy), Q), where each grid pixel contains
                    a distance-weighted sum of quantities of K nearest points to this pixel,
                    using the Epanichnikov kernel with adaptive bandwidth equal to the K'th nearest neighbor;
                    the normalization is the same as for simpleHistogram (i.e. sum divided by the area)
    """
    quantities = numpy.atleast_2d(quantities.T).T  # if the input was a 1d array, make it a Nx1 2d array
    # 2d grid with shape (len(gridx)*len(gridy), 2)
    grid2d = numpy.column_stack((numpy.repeat(gridx, len(gridy)), numpy.tile(gridy, len(gridx))))
    # find out K closest particles for each pixel of the 2d grid
    KDTree = scipy.spatial.cKDTree(point_coords)
    try:      Dngbr, Ingbr = KDTree.query(grid2d, K, workers=nthreads)  # new API
    except:
      try:    Dngbr, Ingbr = KDTree.query(grid2d, K,  n_jobs=nthreads)  # old API
      except: Dngbr, Ingbr = KDTree.query(grid2d, K)  # very old API without multithreading
    bandwidth = Dngbr[:,-1]     # size of the smoothing kernel for each pixel in the 2d grid
    # G is the total number of grid nodes len(grid2d) = len(gridx) * len(gridy)
    # N is the number of input points
    # Q is the number of quantities to average
    # K is the number of points contributing to the average quantities at each grid node
    # quantities: shape (N, Q)  array of quantities for each input point
    # Dngbr:      shape (G, K)  distances to all K points contributing to each grid node
    # Ingbr:      shape (G, K)  indices (range: 0..N-1) of K points contributing to each grid node
    # bandwidth:  shape (G)     kernel bandwidth for each grid node = max(Dngbr, axis=1)
    # result:     shape (G, Q)  average quantities for each grid node
    # quantities[Dngbr]: shape (G, K, Q)  array of quantities for K points contributing to each grid node
    result = 2/numpy.pi * bandwidth[:,None]**-2 * numpy.sum(quantities[Ingbr] *
        (1 - (Dngbr[:,:,None] / bandwidth[:,None,None])**2),  # 2d Epanechnikov kernel
        axis=1)
    if numpy.isfinite(maxdist):
        # discard results for pixels that are more than maxdist away from from any input point [optional]
        result[Dngbr[:,0] > maxdist] = numpy.nan
    return result.reshape(len(gridx), len(gridy), quantities.shape[1])

if __name__ == '__main__':    # demo
    # hint: change npoints to 1e5 and see how the simple histogram improves dramatically
    # while the adaptive KDE was already good at 1e4 points
    npoints = 10000
    # make two gaussian blobs
    nfirst  = int(npoints*0.7)
    x, y  = numpy.random.normal(size=(npoints,2)).T
    x[:nfirst] += 1.5
    y[:nfirst] += 1.0
    x[nfirst:] -= 1.5
    y[nfirst:] -= 1.0
    xy    = numpy.column_stack([x, y])
    # assign unequal weights to points
    mass  = numpy.random.random(size=npoints)
    # assign some other quantities (say, velocity) to each point, differently for the two blobs
    vel1  = numpy.random.normal(size=nfirst)         * 0.5 + (x[:nfirst] + 2*y[:nfirst]) - 1.0
    vel2  = numpy.random.normal(size=npoints-nfirst) * 1.0 + (x[nfirst:]*2 - y[nfirst:])
    vel   = numpy.hstack([vel1, vel2])

    binx  = numpy.linspace(-5, 5, 101)
    biny  = numpy.linspace(-4, 4, 81)
    gridx = (binx[1:]+binx[:-1])/2  # use centres of 2d grid cells as reference points for KDE
    gridy = (biny[1:]+biny[:-1])/2

    quantities = numpy.column_stack([mass, mass*vel, mass*vel**2])
    results    = [
        simpleHistogram  (xy, quantities, gridx, gridy),
        adaptiveHistogram(xy, quantities, gridx, gridy, K=100, maxdist=binx[1]-binx[0])
    ]
    clim_logdensity = [numpy.log10(npoints)-3, numpy.log10(npoints)-1]
    clim_meanvel    = [-10, 10]
    clim_sigmavel   = [0, 2]
    plt.figure(figsize=(16.67,10))
    for p in range(2):
        logdensity = numpy.log10(results[p][:,:,0])
        meanvel    = results[p][:,:,1] / results[p][:,:,0]
        sigmavel   =(results[p][:,:,2] / results[p][:,:,0] - meanvel**2)**0.5
        plt.axes([0.03, 0.05 + (1-p)*0.44, 0.3, 0.4])
        plt.imshow(logdensity.T, extent=[min(binx), max(binx), min(biny), max(biny)],
            aspect='auto', interpolation='nearest', origin='lower', cmap='gray_r',
            vmin=clim_logdensity[0], vmax=clim_logdensity[1])
        plt.text(0.5, 0.95, ['Simple', 'Adaptive'][p]+' histogram',
            ha='center', va='center', transform=plt.gca().transAxes)
        plt.axes([0.36, 0.05 + (1-p)*0.44, 0.3, 0.4])
        plt.imshow(meanvel.T, extent=[min(binx), max(binx), min(biny), max(biny)],
            aspect='auto', interpolation='nearest', origin='lower', cmap='rainbow',
            vmin=clim_meanvel[0], vmax=clim_meanvel[1])
        plt.axes([0.69, 0.05 + (1-p)*0.44, 0.3, 0.4])
        plt.imshow(sigmavel.T, extent=[min(binx), max(binx), min(biny), max(biny)],
            aspect='auto', interpolation='nearest', origin='lower', cmap='rainbow',
            vmin=clim_sigmavel[0], vmax=clim_sigmavel[1])

    # add colorbars
    plt.axes([0.03, 0.95, 0.3, 0.04]).set_yticks([])
    plt.imshow(numpy.linspace(0,1,256).reshape(1,-1), extent=clim_logdensity+[0,1],
        aspect='auto', interpolation='nearest', cmap='gray_r')
    plt.xlabel('log(density)')
    plt.axes([0.36, 0.95, 0.3, 0.04]).set_yticks([])
    plt.imshow(numpy.linspace(0,1,256).reshape(1,-1), extent=clim_meanvel+[0,1],
        aspect='auto', interpolation='nearest', cmap='rainbow')
    plt.xlabel('mean vel')
    plt.axes([0.69, 0.95, 0.3, 0.04]).set_yticks([])
    plt.imshow(numpy.linspace(0,1,256).reshape(1,-1), extent=clim_sigmavel+[0,1],
        aspect='auto', interpolation='nearest', cmap='rainbow')
    plt.xlabel('sigma vel')

    plt.show()
