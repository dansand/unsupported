import numpy as np
import underworld as uw
from underworld import function as fn

def cosine_taper(increasingFn, startval, width):
    cosinusoid= 0.5*(1. - fn.math.cos((np.pi*(increasingFn - startval))/(width)) )
    taperFn = fn.branching.conditional( ((increasingFn < startval, 0. ),
                                            (increasingFn > (startval + width), 1. ), 
                                           (True,                      cosinusoid)  ))
    return taperFn

