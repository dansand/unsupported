from unsupported_dan import interfaces
import numpy as np

def particlesToAdd(markerLine, A, _lowdist, _updist = False):

    all_particle_coords = markerLine.kdtree.data

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

    allcs = markerLine.kdtree.data
    localcs = markerLine.swarm.particleCoordinates.data

    #xmatch = allcs[:,0].searchsorted(localcs[:,0])
    #ymatch = allcs[:,1].searchsorted(localcs[:,1])
    xmatch =np.in1d(allcs[:,0], localcs[:,0])
    ymatch =np.in1d(allcs[:,1], localcs[:,1])


    indexes = xmatch[(xmatch == ymatch)]
    return indexes


def laplaceVector(markerLine, k,  limit=0.25):
    """
    this includes my current appraoch to managing the local/shadow complexity
    We build the Laplacian using the full particle coordinate data (accessed through kdtree.data)
    Then apply a mask to the resulting update vector.


    """

    #First we want to get the distances between neighbours for each point
    A = markerLine.neighbourMatrix( k =k)
    nbIndexes = [ list(np.nonzero(t)[0]) for t in A ]
    for item in nbIndexes:
        if len(item) == 1:
            item.append(item[0])

    n1 = np.array(nbIndexes)[:,0]
    n2 = np.array(nbIndexes)[:,1]
    distances = np.linalg.norm(markerLine.kdtree.data[n2] - markerLine.kdtree.data[n1], axis = 1)

    #Now compute the lapacian
    L = markerLine.laplacianMatrix(k = k)
    dl = np.dot(L,markerLine.kdtree.data)
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
