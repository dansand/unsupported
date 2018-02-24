from easydict import EasyDict as edict
import math
import numpy as np
import pint
import unsupported_dan.geodynamics.scaling as sub_scaling;

u = sub_scaling.UnitRegistry

#####################
#A set of physical parameters, that can be used to run basic subduction simulations
#####################

#dp refers to dimensional paramters
dp = edict({})

#Main physical paramters (thermal convection parameters)
dp.refDensity = 3300.* u.kilogram / u.meter**3                    #reference density
dp.refGravity = 9.8* u.metre / u.second**2                        #surface gravity
dp.refDiffusivity = 1e-6 *u.metre**2 / u.second                   #thermal diffusivity
dp.refExpansivity = 3e-5/u.kelvin                                        #surface thermal expansivity
dp.refViscosity = 1e20* u.pascal* u.second
dp.refLength = 2900*u.km
dp.gasConstant = 8.314*u.joule/(u.mol*u.kelvin)                   #gas constant
dp.specificHeat = 1250.4*u.joule/(u.kilogram* u.kelvin)           #Specific heat (Jkg-1K-1)
dp.potentialTemp = 1573.*u.kelvin                                 #mantle potential temp (K)
dp.surfaceTemp = 273.*u.kelvin                                    #surface temp (K)
#these are the shifted dimensionless temps
dp.potentialTemp_ = dp.potentialTemp - dp.surfaceTemp
dp.surfaceTemp_ = dp.surfaceTemp - dp.surfaceTemp
#dp.deltaT = dp.potentialTemp - dp.surfaceTemp
#main rheology paramters (thermal convection parameters)
dp.cohesionMantle = 20.*u.megapascal                              #mantle cohesion in Byerlee law
dp.frictionMantle = u.Quantity(0.2)                                           #mantle friction coefficient in Byerlee law (tan(phi))
dp.frictionMantleDepth = dp.frictionMantle*u.pascal/u.meter
dp.diffusionPreExp = 5.34e-10/u.pascal/u.second                   #pre-exp factor for diffusion creep
dp.diffusionEnergy = 3e5*u.joule/(u.mol)
dp.diffusionEnergyDepth = 3e5*u.joule/(u.mol*dp.gasConstant)
dp.diffusionVolume=5e-6*u.meter**3/(u.mol)
dp.diffusionVolumeDepth=5e-6*dp.refDensity.magnitude*dp.refGravity.magnitude*u.joule/(u.mol*dp.gasConstant*u.meter)
dp.viscosityInterface = 5e19*u.pascal   * u.second
dp.adiabaticTempGrad = (dp.refExpansivity*dp.refGravity*dp.potentialTemp)/dp.specificHeat


#####################
#Next, define a standard set of scale factors used to non-dimensionalize the system
#####################


KL = dp.refLength
KT = dp.potentialTemp - dp.surfaceTemp
Kt = KL**2/dp.refDiffusivity
KM = dp.refViscosity * KL * Kt

sub_scaling.scaling["[length]"]      = KL.to_base_units()
sub_scaling.scaling["[temperature]"] = KT.to_base_units()
sub_scaling.scaling["[mass]"]        = KM.to_base_units()
sub_scaling.scaling["[time]"] =        Kt.to_base_units()


#####################
#md is a set of model settings, and parameters that can be used to run basic subduction simulations
#####################

md = edict({})
#Model geometry, and misc Lengths used to control behaviour
md.depth=1000*u.km                                                #Model Depth
md.aspectRatio=5.
md.interfaceViscCutoffDepth = 100*u.km
md.interfaceViscEndWidth = 20*u.km
md.faultThickness = 10.*u.km
md.faultDestroyDepth = 500*u.km                                    #interface material (crust) an top of slabs
md.lowerMantleDepth=660.*u.km
#Slab and plate init. parameters
md.subZoneLoc=-100*u.km                                           #X position of subduction zone...km
md.slabInitMaxDepth=150*u.km
md.radiusOfCurv = 350*u.km                                        #radius of curvature
md.slabAge=70*u.megayears                                      #age of subduction plate at trench
md.opAgeAtTrench=35*u.megayears                                        #age of op
#numerical and computation params
md.res=48
md.ppc=25                                                         #particles per cell
md.elementType="Q1/dQ0"
md.refineMesh = True
md.nltol = 0.01
md.druckerAlpha = 1.
md.penaltyMethod = True
md.buoyancyFac = 1.0
md.faultLocFac = 1. #this is the relative location of the fault in terms of the fault thickess from the top of slab
md.viscosityMin = 1e18* u.pascal * u.second
md.viscosityMax = 1e25* u.pascal * u.second
md.lowerMantleViscFac = 30.0
md.yieldStressMax=200*u.megapascal


#####################
#Now we map dp, md to non-nonDimensionalized dictionaries
#####################

paramDict = edict({})
for key, val in dp.items():
    if val.unitless:
        paramDict[key] = val.magnitude
    else:
        paramDict[key] = sub_scaling.nonDimensionalize(val)


md_ = edict({})
for key, val in md.items():
    md_[key] = sub_scaling.nonDimensionalize(val)

modelDict = md_

#####################
#Finally, define some dimensional numbers and scaling factors
#####################

#Important to remember the to_base_units conversion here
rayleighNumber = ((dp.refExpansivity*dp.refDensity*dp.refGravity*(dp.potentialTemp - dp.surfaceTemp)*dp.refLength**3).to_base_units() \
                  /(dp.refViscosity*dp.refDiffusivity).to_base_units()).magnitude

stressScale = ((dp.refDiffusivity*dp.refViscosity)/dp.refLength**2).magnitude
pressureDepthGrad = ((dp.refDensity*dp.refGravity*dp.refLength**3).to_base_units()/(dp.refViscosity*dp.refDiffusivity).to_base_units()).magnitude
