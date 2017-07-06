import numpy as np

def slab_top(trench, normal, gradientFn, ds, maxDepth, mesh):
    """
    Create points representing the top of a slab from trench to maxDepth
    
    Parameter
    ---------
    trench : list or list like 
            Points represnting trench location, 
    normal: list or list like
            vector in the horizontal plane normal to trench
    gradientFn: function
             function that returns the dip or the slab dz/ds 
             where s is the distance along the normal vector
    ds: float
            distance between points, in model coordinates
    
    max depth: float, or list or list like
            Maximum depth of slab
    mesh: uw 2 mesh   
    
    """
    
    #convert everything to numpy arrays
    trench = np.array(trench)
    normal = np.array(normal)/np.linalg.norm(normal)
    maxDepth = np.array(maxDepth)
    
    #test if gradientFn is a function   
    points = []
    points.append(list(trench))
    
    #set starting values
    #normal/= np.linalg.norm(normal)#unitize
    vertical = np.zeros(mesh.dim)
    vertical[-1] = -1.

    P0 = trench.copy()
    F = gradientFn(0.)
    #print(F)
    H = 0.
    V = 0.
    #print(normal, F)
    S = normal.copy()

    S[-1] = F     #assumes the last component is vertical
    S/= np.linalg.norm(S) #unitize

    
    while V < maxDepth:
        
        #Get next set of points
        P1 = P0 + S*ds
        points.append(list(P1))
        P0 = P1
        
        #update, H, V, F, S
        H +=  np.dot(S*ds, normal)
        V +=  abs(((S*ds)[-1]))
        F = gradientFn(H)        
        S = normal.copy()
        S[-1] = F     #assumes the last component is vertical
        S/= np.linalg.norm(S) #unitize

        
        
        
    return(np.array(points))##
