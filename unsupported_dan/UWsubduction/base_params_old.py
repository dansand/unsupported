from easydict import EasyDict as edict
import math
import numpy as np


dp = edict({})
#Main physical paramters
dp.depth=1000e3                         #Depth
dp.refDensity=3300.                        #reference density
dp.refGravity=9.8                          #surface gravity
dp.viscosityScale=1e20                       #reference upper mantle visc.,
dp.refDiffusivity=1e-6                     #thermal diffusivity
dp.refExpansivity=3e-5                     #surface thermal expansivity
dp.gasConstant=8.314                    #gas constant
dp.specificHeat=1250.                   #Specific heat (Jkg-1K-1)
dp.potentialTemp=1573.                  #mantle potential temp (K)
dp.surfaceTemp=273.                     #surface temp (K)
#Rheology - flow law paramters
dp.cohesionMantle=20e6                   #mantle cohesion in Byerlee law
dp.cohesionInterface=1e6                    #crust cohesion in Byerlee law
dp.frictionMantle=0.2                   #mantle friction coefficient in Byerlee law (tan(phi))
dp.frictionInterface=0.02                   #crust friction coefficient
dp.diffusionPreExp=5.34e-10             #1./1.87e9, pre-exp factor for diffusion creep
dp.diffusionEnergy=3e5
dp.diffusionVolume=5e-6
dp.lowerMantlePreExp=4.23e-15           #1./2.36e14
dp.lowerMantleEnergy=2.0e5
dp.lowerMantleVolume=1.5e-6
dp.lowerMantleViscFac = 30.
#Fk approach to interface viscosity
dp.delViscInterface = 1e4
dp.refViscInterface = 1e20
dp.refDepthInterface = 80e3
#power law creep params
dp.powerLawStrain = 1e-15
dp.powerLawExp = 3.5

#Rheology - cutoff values
dp.viscosityMin=1e18
dp.viscosityMax=1e25                #viscosity max in the mantle material
dp.viscosityMinInterface=1e18               #viscosity min in the weak-crust / interface material
dp.viscosityMaxInterface=1e25               #viscosity min in the weak-crust / interface material
dp.ysMaxInterface = 50e6
dp.yieldStressMax=300*1e6              #
dp.interfaceViscCutoffDepth = 100e3
dp.interfaceViscEndWidth = 20e3

#Intrinsic Lengths
dp.mantleCrustDepth=10.*1e3              #Crust depth
dp.faultThickness = 10.*1e3              #interface material (crust) an top of slabs
dp.lowerMantleDepth=660.*1e3
dp.interfaceDestroyDepth=400.*1e3             #Beyond this depth, crust, interface, etc get destroyed
dp.interfaceResetDepth=400.*1e3               #Beyond this depth there is no proximity accumulation()

#Slab and plate init. parameters
dp.subZoneLoc=-100e3                    #X position of subduction zone...km
dp.maxDepth=150e3
dp.radiusOfCurv = 350e3                          #radius of curvature
dp.slabMaxAge=70e6                     #age of subduction plate at trench
dp.opMaxAge=35e6                       #age of op
dp.spAgeGrad = (-1./4)  # (y / cm.) The gradient direction is taken as away from the sz
dp.opAgeGrad = (-1./8)  # (y / cm.) The gradient direction is taken as away from the sz
#Misc
dp.stickyAirDepth=100e3                 #depth of sticky air layer
dp.viscosityStickyAir=1e19              #stick air viscosity, normal
#derived params
dp.deltaTemp = dp.potentialTemp-dp.surfaceTemp
dp.tempGradMantle = (dp.refExpansivity*dp.refGravity*(dp.potentialTemp))/dp.specificHeat
dp.tempGradSlab = (dp.refExpansivity*dp.refGravity*(dp.surfaceTemp + 400.))/dp.specificHeat

#a velocity scale, used in some caes to set an approx viscosity in the interface
dp.subVelocity = 4*(1/100.)*(1./(3600*24*365)) #m/s

#op velocity will be set as a BC, or False
dp.opVelocity = -1.*(1/100.)*(1./(3600*24*365))


#Modelling and Physics switches

md = edict({})
md.refineMeshStatic=True
md.stickyAir=False
md.aspectRatio=5.
md.res=48
md.ppc=35                                 #particles per cell
md.elementType="Q1/dQ0"
#md.elementType="Q2/DPC1"
md.secInvFac=math.sqrt(1.)
md.courantFac=0.5                         #extra limitation on timestepping
md.thermal = True                        #thermal system or compositional
md.swarmInitialFac = 0.6                 #initial swarm layout will be int(md.ppc*md.swarmInitialFac), popControl will densify later
md.compBuoyancy = False
md.nltol = 0.01
md.maxSteps = 20000
md.checkpointEvery = 200
md.swarmUpdate = 10
md.druckerAlpha = 0.
md.druckerAlphaFault = 0.
md.penaltyMethod = True
md.opuniform = False
md.spuniform = False
md.opfixed = False
md.spfixed = False
md.buoyancyFac = 1.0
md.restartParams = True #whether we load params on checkpoint restart. md.restartParams = True #whether we load params on checkpoint restart.
#The following are time-based actions
md.filesMy = 1.0e6 #dimensional time interval to write files
md.diffuseInitial =  0. #5e6 # years to run initial diffusion for. Or set to zero.

#some flags mostlt controlling the interface viscosity.
md.interfaceType = 2 #1 for 'weak' crust, 2 for 'fault ' (isotropic)
md.plasticInterface = 1 #non linear (plastic) of linear (effective viscosity)
md.faultHeal=1
md.faultLocFac = 1. #this is the relative location of the fault in terms of the fault thickess from the top of slab
md.powerLaw = False
md.viscCombine = 1 #zero here will use min(), anythin else toggles harmonic
md.fixedRidges = False
md.phaseBuoyancy = False



sf = edict({})

sf.lengthScale=2900e3
sf.viscosityScale = dp.viscosityScale
sf.stress = (dp.refDiffusivity*sf.viscosityScale)/sf.lengthScale**2
#sf.lithGrad = dp.refDensity*dp.refGravity*(sf.lengthScale)**3/(sf.viscosityScale*dp.refDiffusivity)
sf.lithGrad = (sf.viscosityScale*dp.refDiffusivity) /(dp.refDensity*dp.refGravity*(sf.lengthScale)**3)
sf.velocity = dp.refDiffusivity/sf.lengthScale
sf.strainRate = dp.refDiffusivity/(sf.lengthScale**2)
sf.time = 1./sf.strainRate
sf.actVolume = (dp.gasConstant*dp.deltaTemp)/(dp.refDensity*dp.refGravity*sf.lengthScale)
sf.actEnergy = (dp.gasConstant*dp.deltaTemp)
sf.diffusionPreExp = 1./sf.viscosityScale
sf.deltaTemp  = dp.deltaTemp
sf.pressureDepthGrad = (dp.refDensity*dp.refGravity*sf.lengthScale**3)/(dp.viscosityScale*dp.refDiffusivity)


#dimesionless params
ndp  = edict({})

ndp.rayleigh = md.buoyancyFac*(dp.refExpansivity*dp.refDensity*dp.refGravity*dp.deltaTemp*sf.lengthScale**3)/(dp.viscosityScale*dp.refDiffusivity)

#Take care with these definitions,
ndp.surfaceTemp = dp.surfaceTemp/sf.deltaTemp  #Ts
ndp.potentialTemp = dp.potentialTemp/sf.deltaTemp - ndp.surfaceTemp #Tp' = Tp - TS
ndp.tempGradMantle = dp.tempGradMantle/(sf.deltaTemp/sf.lengthScale)
ndp.tempGradSlab = dp.tempGradSlab/(sf.deltaTemp/sf.lengthScale)

#lengths / distances
ndp.depth = dp.depth/sf.lengthScale
ndp.leftLim = -0.5*ndp.depth*md.aspectRatio
ndp.rightLim = 0.5*ndp.depth*md.aspectRatio
ndp.faultThickness = dp.faultThickness/sf.lengthScale
ndp.mantleCrustDepth =  dp.mantleCrustDepth/sf.lengthScale
ndp.interfaceDestroyDepth = dp.interfaceDestroyDepth/sf.lengthScale
ndp.interfaceResetDepth = dp.interfaceResetDepth/sf.lengthScale
ndp.lowerMantleDepth = dp.lowerMantleDepth/sf.lengthScale


#times - for convenience and sanity the dimensional values are in years, conversion to seconds happens here
ndp.slabMaxAge =  dp.slabMaxAge*(3600*24*365)/sf.time
ndp.opMaxAge = dp.opMaxAge*(3600*24*365)/sf.time

#Rheology - flow law paramters
ndp.cohesionMantle=dp.cohesionMantle/sf.stress                  #mantle cohesion in Byerlee law
ndp.cohesionInterface=dp.cohesionInterface/sf.stress                  #crust cohesion in Byerlee law
ndp.frictionMantle=dp.frictionMantle/sf.lithGrad                  #mantle friction coefficient in Byerlee law (tan(phi))
ndp.frictionInterface=dp.frictionInterface/sf.lithGrad                  #crust friction coefficient
ndp.diffusionPreExp=dp.diffusionPreExp/sf.diffusionPreExp                #pre-exp factor for diffusion creep
ndp.diffusionEnergy=dp.diffusionEnergy/sf.actEnergy
ndp.diffusionVolume=dp.diffusionVolume/sf.actVolume
ndp.lowerMantlePreExp=dp.lowerMantlePreExp/sf.diffusionPreExp
ndp.lowerMantleEnergy=dp.lowerMantleEnergy/sf.actEnergy
ndp.lowerMantleVolume=dp.lowerMantleVolume/sf.actVolume
ndp.yieldStressMax=dp.yieldStressMax/sf.stress
#Fk approach to interface viscosity
ndp.logDelVisc = np.log(dp.delViscInterface)
ndp.refViscInterface = dp.refViscInterface/sf.viscosityScale
ndp.refDepthInterface = dp.refDepthInterface/sf.lengthScale
#
ndp.powerLawStrain = dp.powerLawStrain/sf.strainRate
ndp.powerLawExp = dp.powerLawExp

#Rheology - cutoff values
ndp.viscosityMin= dp.viscosityMin /sf.viscosityScale
ndp.viscosityMax=dp.viscosityMax/sf.viscosityScale
ndp.viscosityMinInterface= dp.viscosityMinInterface /sf.viscosityScale
ndp.viscosityMaxInterface= dp.viscosityMaxInterface /sf.viscosityScale
ndp.lowerMantleViscFac = dp.lowerMantleViscFac
ndp.interfaceViscCutoffDepth = dp.interfaceViscCutoffDepth/sf.lengthScale
ndp.interfaceViscEndWidth = dp.interfaceViscEndWidth/sf.lengthScale
ndp.ysMaxInterface  = dp.ysMaxInterface/sf.stress


#Slab and plate init. parameters
ndp.subZoneLoc = dp.subZoneLoc/sf.lengthScale
ndp.maxDepth = dp.maxDepth/sf.lengthScale
ndp.radiusOfCurv = dp.radiusOfCurv/sf.lengthScale
ndp.subVelocity = dp.subVelocity/sf.velocity
ndp.opVelocity = dp.opVelocity/sf.velocity
ndp.spAgeGrad = dp.spAgeGrad* (3600*24*365)*(100.)*(sf.lengthScale/sf.time)
ndp.opAgeGrad = dp.opAgeGrad* (3600*24*365)*(100.)*(sf.lengthScale/sf.time)
