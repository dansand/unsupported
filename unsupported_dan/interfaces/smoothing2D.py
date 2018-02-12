#
from unsupported_dan import interfaces
import numpy as np

def particlesToAdd(markerLine, A, _lowdist, _updist = False):

    #Build all objects that need to be called collectively
    all_particle_coords = markerLine.data
    pd = markerLine.pairDistanceMatrix()
    #We want only the lower half of the matrix, including the upper half would add particles twice
    Alow = np.tril(A)

    if all_particle_coords.shape[1]:
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

    else:
        newPoints = np.empty((2,0))

    return newPoints


def shadowMask(markerLine):

    """
    Builds a boolean mask that filters only the local particles,
    from an markerLine particle array that includes both shadow and local
    """


    #allcs = markerLine.kdtree.data
    allcs = markerLine.data
    localcs = markerLine.swarm.particleCoordinates.data

    indexes = np.empty((0,)).astype('bool')

    if localcs.shape[0]:
        #indexes = np.full((localcs.shape[0]), True, dtype=bool)
        xmatch =np.in1d(allcs[:,0], localcs[:,0])
        ymatch =np.in1d(allcs[:,1], localcs[:,1])
        indexes = xmatch[(xmatch == ymatch)]

    return indexes


def laplaceVector(markerLine, k,  limit=0.25):
    """
    this includes my current approach to managing the local/shadow complexity
    We build the Laplacian using the full particle coordinate data (accessed through markerLine.data)
    Then apply a mask to the resulting update vector.


    """

    #First we want to get the distances between neighbours for each point
    A = markerLine.neighbourMatrix( k =k)
    #and we only want to return the vector corresponding to local points
    shmask = shadowMask(markerLine)
    all_particle_coords = markerLine.data
    #build the lapacian matrix
    L = markerLine.laplacianMatrix(k = k)

    dlF = np.empty((0,2))

    if len(all_particle_coords) == 1: #safeguard for little arrays
        return dlF

    if all_particle_coords.shape[1]:
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
        distances = np.linalg.norm(all_particle_coords[n2] - all_particle_coords[n1], axis = 1)


        #dl = np.dot(L,markerLine.kdtree.data)
        dl = np.dot(L,all_particle_coords)
        dlNorms = np.linalg.norm(dl, axis = 1)

        #we only want to dampen the Laplacian if returns high values
        #try to mitigate the damage if bad neighbour information is passed through
        #lapalacian vectors should never be bigger than limit*distances
        distanceMask = dlNorms > limit*distances
        #x and y components
        dl[distanceMask,0] *= (limit*distances[distanceMask]/dlNorms[distanceMask])
        dl[distanceMask,1] *= (limit*distances[distanceMask]/dlNorms[distanceMask])
        dlF = dl[shmask]

    return dlF

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

    newPoints = np.empty((0,2))
    localIds = np.empty((0,)).astype('bool')

    #these methods need to be called by all procs
    all_particle_coords = markerLine.data
    shmask = shadowMask(markerLine)
    Alow = np.tril(A) #We want only the lower half of the matrix,
    pd = markerLine.pairDistanceMatrix()

    #Now we can isolate the processors with particles
    if all_particle_coords.shape[1]:

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

        AF = Alow[shmask]*pdMask[shmask]
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
    #Smoothing
    ###########
    for cyc in range(smoothCycles):
        markerLine.rebuild()
        Dl = laplaceVector(markerLine, k = k, limit=laplaceLimit)

        #print(markerLine.swarm.particleCoordinates.data[:].shape, Dl.shape)
        #if Dl.shape[0]

        with markerLine.swarm.deform_swarm():
                markerLine.swarm.particleCoordinates.data[:] -= _lambda *Dl

    ###########
    #Addition
    ###########

    #these methods need to be called by all procs
    markerLine.rebuild()
    A = markerLine.neighbourMatrix( k =k)
    newPoints = particlesToAdd(markerLine, A, _lowdist=2.*ds) #start addding particles above _lowdist

    #Now we can isolate the processors that actually have points
    if markerLine.empty:
        newPoints = np.empty((0,2))
    markerLine.add_points(newPoints[:,0], newPoints[:,1])

    ###########
    #Removal
    ###########

    #dummy arrays to use in case there's no markerLine on the proc
    midPoints = np.empty((0,2))
    currentIds = np.empty((0,)).astype('bool')
    #these methods need to be called by all procs
    markerLine.rebuild()
    A = markerLine.neighbourMatrix( k =k)
    midPoints, currentIds = neighbourDistanceQuery(markerLine, A, _lowdist=0.,_updist= 0.5*ds)

    with markerLine.swarm.deform_swarm():
        markerLine.swarm.particleCoordinates.data[currentIds] = (9999999., 9999999.)

    markerLine.add_points(midPoints[:,0], midPoints[:,1])
