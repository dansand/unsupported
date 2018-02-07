#
from unsupported_dan import interfaces
import numpy as np

def particlesToAdd(markerLine, A, _lowdist, _updist = False):

    #all_particle_coords = markerLine.kdtree.data
    all_particle_coords = markerLine.data

    #We want only the lower half of the matrix, including the upper half would add particles twice
    Alow = np.tril(A)

    pd = markerLine.pairDistanceMatrix()

    #Here is the distance mask
    if _updist:
        pdMask = np.logical_and(pd > _lowdist, pd < _updist)
    else:
        pdMask = pd > _lowdist

    #We only want to choose those particles that have two nearest neighbours (this hopefully excludes endpoints)
    mask = np.where(A.sum(axis=1) != 2)
    #Set those rows to zero
    pdMask[mask,:]  = 0

    #the magic is here - simply mutiply the neigbours matrix by the distance mask
    AF = Alow*pdMask

    uniques = np.transpose(np.nonzero(AF))
    #First, store a complete copy of the new particle positions (mean pair positions)
    newPoints = np.copy(0.5*(all_particle_coords[uniques[:,0]] + all_particle_coords[uniques[:,1]]))

    return newPoints


def shadowMask(markerLine):

    """
    Builds a boolean mask that filters only the local particles,
    from an markerLine particle array that includes both shadow and local
    """

    #this function relies on updated information
    #markerLine.rebuild()


    #allcs = markerLine.kdtree.data
    allcs = markerLine.data
    localcs = markerLine.swarm.particleCoordinates.data

    xmatch =np.in1d(allcs[:,0], localcs[:,0])
    ymatch =np.in1d(allcs[:,1], localcs[:,1])

    indexes = xmatch[(xmatch == ymatch)]
    return indexes


def laplaceVector(markerLine, k,  limit=0.25):
    """
    this includes my current appraoch to managing the local/shadow complexity
    We build the Laplacian using the full particle coordinate data (accessed through markerLine.data)
    Then apply a mask to the resulting update vector.


    """

    #First we want to get the distances between neighbours for each point
    A = markerLine.neighbourMatrix( k =k)
    nbIndexes = [ list(np.nonzero(t)[0]) for t in A ]
    for item in nbIndexes:
        if len(item) == 1:
            item.append(item[0])
        #Adding this for case of empty items
        if len(item) == 0:
            item.append(0)
            item.append(0)

    n1 = np.array(nbIndexes)[:,0]
    n2 = np.array(nbIndexes)[:,1]
    #distances = np.linalg.norm(markerLine.kdtree.data[n2] - markerLine.kdtree.data[n1], axis = 1)
    distances = np.linalg.norm(markerLine.data[n2] - markerLine.data[n1], axis = 1)

    #Now compute the lapacian
    L = markerLine.laplacianMatrix(k = k)
    #dl = np.dot(L,markerLine.kdtree.data)
    dl = np.dot(L,markerLine.data)
    dlNorms = np.linalg.norm(dl, axis = 1)

    #we only want to dampen the Laplacian if returns high values
    #try to mitigate the damage if bad neighbour information is passed through
    #lapalacian vectors should never be bigger than limit*distances
    distanceMask = dlNorms > limit*distances
    #x and y components
    dl[distanceMask,0] *= (limit*distances[distanceMask]/dlNorms[distanceMask])
    dl[distanceMask,1] *= (limit*distances[distanceMask]/dlNorms[distanceMask])

    #and we only want to return the vector corresponding to local points
    mask = shadowMask(markerLine)


    return dl[mask]

def neighbourDistanceQuery(markerLine, A, _lowdist=1e-10, _updist = False):

    """
    Queries the marker line and returns info about particles,
    where the inter-particle distance fulfils the distance conditions

    Returns:
    newPoints: midpoints of any pairs of particles fulfilling the conditions.
    These are generated with both the local and shadow info,
    meaning that some of the returned particles may no be on the local processor.

    localIds: the ids of local part of the markerLine where the conditions were satisfied.
    These are intended for deletion, so only the local part of the markerLine is wanted.

    """

    #all_particle_coords = markerLine.kdtree.data
    all_particle_coords = markerLine.data

    #We want only the lower half of the matrix, including the upper half would add particles twice
    Alow = np.tril(A)

    pd = markerLine.pairDistanceMatrix()

    #Here is the distance mask
    if _updist:
        pdMask = np.logical_and(pd > _lowdist, pd < _updist)
    else:
        pdMask = pd > _lowdist

    #We only want to choose those particles that have two nearest neighbours (this hopefully excludes endpoints)
    mask = np.where(A.sum(axis=1) != 2)
    #Set those rows to zero
    pdMask[mask,:]  = 0

    #to get the ids simply mutiply the neigbours matrix by the distance mask
    #this bit works from the full set of points local + shadow
    AF = Alow*pdMask
    allIds = np.transpose(np.nonzero(AF))
    newPoints = np.copy(0.5*(all_particle_coords[allIds[:,0]] + all_particle_coords[allIds[:,1]]))

    #Now we want to get the local ids - typically use these for removing points
    mask = shadowMask(markerLine)
    AF = Alow[mask]*pdMask[mask]
    localIds = np.transpose(np.nonzero(AF))

    return newPoints, localIds.flatten()


def repair_markerLines(markerLine, ds, smoothCycles=1, k=4, _lambda = 0.5, laplaceLimit = 0.25 ):

    """
    smoothCycles ...
    k = max number of particles to search for neighbour information
    _lambda = 0.5         #A dampening  applied to the entire laplacian vector
    laplaceLimit = 0.25   #fraction of inter-particle distance that is the maximum laplace displacement
    """

    ###########
    #Removal
    ###########
    #dummy arrays to use in case there's no markerLine on the proc

    midPoints = np.empty((0,2))
    currentIds = np.empty((0,)).astype('bool')

    #***I'm pairing these, because the shadow mask Fn only works is line is rebuilt.
    A = markerLine.neighbourMatrix(k =k, jitter=1e-8)
    markerLine.rebuild()
    #***
    if not markerLine.empty:
        midPoints, currentIds = neighbourDistanceQuery(markerLine, A, _lowdist=0.,_updist= 0.5*ds)
        #Need to delete those points first, before the
    with markerLine.swarm.deform_swarm():
        markerLine.swarm.particleCoordinates.data[currentIds] = (9999999., 9999999.)

    markerLine.add_points(midPoints[:,0], midPoints[:,1])


    ###########
    #Smoothing
    ###########
    for cyc in range(smoothCycles):
        markerLine.rebuild()
        if not markerLine.empty:
            Dl = laplaceVector(markerLine, k = k, limit=laplaceLimit)
        else:
            Dl = 0.0

        with markerLine.swarm.deform_swarm():
                markerLine.swarm.particleCoordinates.data[:] -= _lambda *Dl

    ###########
    #Addition
    ###########
    markerLine.rebuild()
    if not markerLine.empty:
        A = markerLine.neighbourMatrix( k =k)
        newPoints = particlesToAdd(markerLine, A, _lowdist=2.*ds) #start addding particles above _lowdist
    else:
        newPoints = np.empty((0,2))
    markerLine.add_points(newPoints[:,0], newPoints[:,1])
