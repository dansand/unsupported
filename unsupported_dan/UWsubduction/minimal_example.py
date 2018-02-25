from easydict import EasyDict as edict
import math
import numpy as np
import pint
import unsupported_dan.geodynamics.scaling as sub_scaling;

UnitRegistry = sub_scaling.UnitRegistry
u = UnitRegistry

#####################
#A set of physical parameters, that can be used to run basic subduction simulations
#####################

#pd refers to dimensional paramters
pd = edict({})

#Main physical paramters (thermal convection parameters)
pd.refDensity = 3300.* u.kilogram / u.meter**3                    #reference density
pd.refGravity = 9.8* u.metre / u.second**2                        #surface gravity
pd.refDiffusivity = 1e-6 *u.metre**2 / u.second                   #thermal diffusivity
pd.refExpansivity = 3e-5/u.kelvin                                 #surface thermal expansivity
pd.refViscosity = 1e20* u.pascal* u.second
pd.refLength = 2900*u.km
pd.gasConstant = 8.314*u.joule/(u.mol*u.kelvin)                   #gas constant
pd.specificHeat = 1250.4*u.joule/(u.kilogram* u.kelvin)           #Specific heat (Jkg-1K-1)
pd.potentialTemp = 1573.*u.kelvin                                 #mantle potential temp (K)
pd.surfaceTemp = 273.*u.kelvin                                    #surface temp (K)
#these are the shifted temps, which will range from 0 - 1 in the dimensionless system
pd.potentialTemp_ = pd.potentialTemp - pd.surfaceTemp
pd.surfaceTemp_ = pd.surfaceTemp - pd.surfaceTemp
#main rheology parameters
pd.cohesionMantle = 20.*u.megapascal                              #mantle cohesion in Byerlee law
pd.frictionMantle = u.Quantity(0.2)                                           #mantle friction coefficient in Byerlee law (tan(phi))
pd.frictionMantleDepth = pd.frictionMantle*pd.refDensity*pd.refGravity
pd.diffusionPreExp = 5.34e-10/u.pascal/u.second                   #pre-exp factor for diffusion creep
pd.diffusionEnergy = 3e5*u.joule/(u.mol)
pd.diffusionEnergyDepth = 3e5*u.joule/(u.mol*pd.gasConstant)
pd.diffusionVolume=5e-6*u.meter**3/(u.mol)
pd.diffusionVolumeDepth=5e-6*pd.refDensity.magnitude*pd.refGravity.magnitude*u.joule/(u.mol*pd.gasConstant*u.meter)
pd.viscosityInterface = 5e19*u.pascal   * u.second
pd.adiabaticTempGrad = (pd.refExpansivity*pd.refGravity*pd.potentialTemp)/pd.specificHeat
pd.yieldStressMax=200*u.megapascal
pd.lowerMantleViscFac = u.Quantity(30.0)


#####################
#md is a set of model settings, and parameters that can be used to run basic subduction simulations
#####################

md = edict({})
#Model geometry, and misc Lengths used to control behaviour
md.depth=1000*u.km                                                #Model Depth
md.aspectRatio=5.
#lengths, factors relating to subduction interface behaviour
md.interfaceViscDepthTaperStart = 100*u.km
md.interfaceViscDepthTaperWidth = 20*u.km
md.interfaceViscHorizTaperStart = 300*u.k
md.interfaceViscHorizTaperWidth = 300*u.k
md.faultThickness = 10.*u.km
md.faultLocFac = 1.                                                #this is the relative location of the fault in terms of the fault thickess from the top of slab
md.faultDestroyDepth = 500*u.km                                    #interface material (crust) an top of slabs
md.lowerMantleDepth=660.*u.km
md.lowerMantleTransWidth=10.*u.km
#Slab and plate init. parameters
md.subZoneLoc=-100*u.km                                           #X position of subduction zone...km
md.slabInitMaxDepth=150*u.km
md.radiusOfCurv = 350.*u.km                                        #radius of curvature
md.slabAge=70.*u.megayears                                      #age of subduction plate at trench
md.opAgeAtTrench=35.*u.megayears                                        #age of op
#numerical and computation params
md.res=48
md.ppc=25                                                         #particles per cell
md.elementType="Q1/dQ0"
md.refineMesh = True
md.meshRefineFactor = 0.7
md.nltol = 0.01
md.druckerAlpha = 1.
md.penaltyMethod = True
md.buoyancyFac = 1.0
md.viscosityMin = 1e18* u.pascal * u.second
md.viscosityMax = 1e25* u.pascal * u.second

#####################
#Next, define a standard set of scale factors used to non-dimensionalize the system
#####################


KL = pd.refLength
KT = pd.potentialTemp - pd.surfaceTemp
Kt = KL**2/pd.refDiffusivity
KM = pd.refViscosity * KL * Kt

sub_scaling.scaling["[length]"]      = KL.to_base_units()
sub_scaling.scaling["[temperature]"] = KT.to_base_units()
sub_scaling.scaling["[mass]"]        = KM.to_base_units()
sub_scaling.scaling["[time]"] =        Kt.to_base_units()


#####################
#Now we map pd, md to non-nonDimensionalized dictionaries, paramDict, modelDict
#####################

paramDict  = edict({})

for key, val in pd.items():
    if val.unitless:
        paramDict [key] = val.magnitude
    else:
        paramDict [key] = sub_scaling.nonDimensionalize(val)



modelDict= edict({})
for key, val in md.items():
    modelDict[key] = sub_scaling.nonDimensionalize(val)


#####################
#Finally, define some dimensional numbers and scaling factors
#####################

#Important to remember the to_base_units conversion here
rayleighNumber = ((pd.refExpansivity*pd.refDensity*pd.refGravity*(pd.potentialTemp - pd.surfaceTemp)*pd.refLength**3).to_base_units() \
                  /(pd.refViscosity*pd.refDiffusivity).to_base_units()).magnitude

stressScale = ((pd.refDiffusivity*pd.refViscosity)/pd.refLength**2).magnitude
pressureDepthGrad = ((pd.refDensity*pd.refGravity*pd.refLength**3).to_base_units()/(pd.refViscosity*pd.refDiffusivity).to_base_units()).magnitude
