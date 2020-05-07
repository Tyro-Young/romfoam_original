#!/usr/bin/env python
"""
ROMFoam run script for the Ahmed body case
"""

# =================================================================================================
# Imports
# =================================================================================================
import os,time
import argparse
import sys
import numpy as np
from mpi4py import MPI
from baseclasses import *
from romfoam import *
from pygeo import *
from pyspline import *
from idwarp import *


# =============================================================================
# Input Parameters
# =============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("--output", help='Output directory', type=str,default='../output')
parser.add_argument("--task", help="type of run to do", type=str, default='run')
parser.add_argument('--optVars',type=str,help='Vars for the optimizer',default="['rampAngle']")
parser.add_argument('--sample',type=int,help='which sample DV to run',default=1)
parser.add_argument('--mode',type=str,help='can be either train or predict',default='train')
parser.add_argument('--runEndTime',type=int,help='number of time steps for flow',default=500)
parser.add_argument('--nSamples',type=int,help='number of samples',default=1)
args = parser.parse_args()
exec('optVars=%s'%args.optVars)
task = args.task
outputDirectory = args.output
gcomm = MPI.COMM_WORLD

sample=args.sample-1

if gcomm.rank==0:
    print((sample,optVars))

# NOTE: put the ref angle at the end of this list because we will use startFrom latestTime in system/controlDict 
if optVars[0]=='rampAngle':
    # 5 samples, 1 DV
    DVs_Train=   [[5.0],
                  [10.0],
                  [15.0],
                  [25.0],
                  [20.0]]
    DVs_Predict= [[21.0]]
elif optVars[0]=='rideHeight':
    # 5 samples, 1 DV
    DVs_Train=   [[0.03],
                  [0.04],
                  [0.06],
                  [0.07],
                  [0.05]]
    DVs_Predict= [[0.07],
                  [0.03]]
elif optVars[0]=='rampAngleAndRideHeight':
    # 10 samples, 2 DVs
    DVs_Train=   [[24.1,0.045],
                  [25.3,0.049],
                  [25.7,0.043],
                  [25.1,0.059],
                  [25.5,0.055],
                  [24.5,0.041],
                  [24.3,0.053],
                  [24.7,0.057],
                  [25.9,0.047],
                  [24.9,0.051]]
    DVs_Predict= [[25.7,0.043],
                  [24.1,0.045]]
elif optVars[0]=='shape':
    # 20 samples, 4 DVs
    DVs_Train=   [[0.00625,0.04625,0.04875,0.02875],
                  [0.01625,0.02625,0.04375,0.03375],
                  [0.03125,0.02125,0.02375,0.01875],
                  [0.03375,0.01375,0.03625,0.00125],
                  [0.00125,0.00875,0.02125,0.00375],
                  [0.04625,0.00625,0.00875,0.01125],
                  [0.01125,0.03625,0.01125,0.03625],
                  [0.03875,0.00375,0.04625,0.04375],
                  [0.02625,0.02375,0.03375,0.04875],
                  [0.04375,0.03375,0.02625,0.01625],
                  [0.02125,0.03875,0.01625,0.02375],
                  [0.04125,0.04375,0.03875,0.03875],
                  [0.00375,0.04125,0.01875,0.00875],
                  [0.04875,0.01875,0.03125,0.04625],
                  [0.02875,0.00125,0.04125,0.04125],
                  [0.02375,0.04875,0.01375,0.00625],
                  [0.03625,0.01625,0.00375,0.01375],
                  [0.01875,0.03125,0.00625,0.03125],
                  [0.00875,0.01125,0.00125,0.02625],
                  [0.01375,0.02875,0.02875,0.02125],
                  [0.00125,0.00875,0.02125,0.00375]]
    DVs_Predict= [[0.00875,0.01125,0.00125,0.02625],
                  [0.02125,0.03875,0.01625,0.02375]]
else:
    print("optVars not valid")
    exit(1)

DVs_Train=np.asarray(DVs_Train)
DVs_Predict=np.asarray(DVs_Predict)

if args.mode=='train':
    DVs0=DVs_Train[sample]
elif args.mode=='predict':
    DVs0=DVs_Predict[sample]
else:
    print("mode should be either train or predict")
    exit(1)

if args.mode=='predict':
    deltaDVs=DVs_Predict[sample]-DVs_Train[-1]
    deltaDVs=deltaDVs.tolist()
elif args.mode=='train': 
    deltaDVs=[0.0]


# Set the parameters for optimization
aeroOptions = {
    # output options
    'casename':                 'AhmedBody_'+task+'_'+optVars[0],
    'outputdirectory':          outputDirectory,
    'writesolution':            True,
    'usecoloring':              False,
    'printalloptions':          False,

    # design surfaces and cost functions 
    'designsurfacefamily':     'designSurfaces', 
    'designsurfaces':          ['body'], 
    'objfuncs':                ['CD'],
    'objfuncgeoinfo':          [['body']],
    'referencevalues':         {'magURef':20.0,'ARef':0.056016,'rhoRef':1.0,'pRef':0.0,'LRef':1.0},

    # flow setup
    'adjointsolver':           'simpleROMFoam',
    'rasmodel':                'SpalartAllmarasFv3',
    'flowcondition':           'Incompressible',
    'maxflowiters':            args.runEndTime, 
    'writeinterval':           args.runEndTime,
    'avgobjfuncs':             False,
    'avgobjfuncsstart':        2000,
    'setflowbcs':              False, 
    'inletpatches':            ['inlet'],
    'outletpatches':           ['outlet'],
    'flowbcs':                 {'bc0':{'patch':'inlet','variable':'U','value':[20.0,0.0,0.0]},
                                'useWallFunction':'true'},                
    # romDict
    'nsamples':                args.nSamples,
    'deltaffd':                deltaDVs,
    'svdtype':                 'cross',
    'svdtol':                  1e-8,
    'svdmaxits':               100,
    'svdrequestedn':           args.nSamples,
    'usemf':                   1,
    'mfstep':                  1e-6,
    'debugmode':               0,
    'uselspg':                 0,
    'romnkabstol':             1e-8,

    # adjoint setup
    'adjdvtypes':              ['FFD'], 
    'epsderivffd':             1.0e-3,
    'adjjacmatordering':       'state',
    'adjjacmatreordering':     'rcm',
    'normalizestates':         [],
    'normalizeresiduals':      ['URes','pRes','nuTildaRes'],    
    
    ########## misc setup ##########
    'mpispawnrun':             False,

}

# mesh warping parameters, users need to manually specify the symmetry plane
meshOptions = {
    'gridFile':                os.getcwd(),
    'fileType':                'openfoam',
    # point and normal for the symmetry plane
    'symmetryPlanes':          [[[0.,0., 0.],[0., 0., 0.]]], 
}

# =================================================================================================
# DVGeo
# =================================================================================================
FFDFile = './FFD/bodyFittedFFD.xyz'
DVGeo = DVGeometry(FFDFile)
x = [0.000,0.100,0.843,1.044]
y = [0.100,0.100,0.100,0.100]
z = [0.194,0.194,0.194,0.147]
c1 = pySpline.Curve(x=x, y=y, z=z, k=2)
DVGeo.addRefAxis('bodyAxis', curve = c1,axis='z')


def rampAngle(val,geo):

    C = geo.extractCoef('bodyAxis')

    # the value will be ramp angle in degree.
    # start with a conversion to rads
    angle = (val[0])*np.pi/180.0

    # Overall length needs to stay a 1.044, so use that as a ref for
    # the final mesh point

    # set the target length
    lTarget = 0.222
    hInit = 0.246 - 0.05

    # compute the coefficient deltas
    dx = lTarget*np.cos(angle)
    dz = lTarget*np.sin(angle)

    topEdge = 0.338-dz
    rearHeight = topEdge-0.05
    coefPoint = rearHeight/2.0 +0.05
    scalez = rearHeight/hInit

    # Set the coefficients
    C[3,0] = 1.044
    C[2,0] = C[3,0]-dx
    C[2,2] = 0.194
    C[3,2] = coefPoint

    geo.restoreCoef(C, 'bodyAxis')

    geo.scale_z['bodyAxis'].coef[3] = scalez

    return

def rampAngleAndRideHeight(val,geo):
    
    C = geo.extractCoef('bodyAxis')

    # the value will be ramp angle in degree.
    # start with a conversion to rads
    angle = (val[0])*np.pi/180.0

    # Overall length needs to stay a 1.044, so use that as a ref for
    # the final mesh point

    # set the target length
    lTarget = 0.222
    hInit = 0.246 - 0.05

    # compute the coefficient deltas
    dx = lTarget*np.cos(angle)
    dz = lTarget*np.sin(angle)

    topEdge = 0.338-dz
    rearHeight = topEdge-0.05
    coefPoint = rearHeight/2.0 +0.05
    scalez = rearHeight/hInit

    # Set the coefficients
    C[3,0] = 1.044
    C[2,0] = C[3,0]-dx
    C[2,2] = 0.194
    C[3,2] = coefPoint

    geo.restoreCoef(C, 'bodyAxis')

    geo.scale_z['bodyAxis'].coef[3] = scalez

    rideHeight=val[1]
    for i in range(len(C)):
        if i>=2:
            C[i,2]=C[i,2]+(rideHeight-0.05)*scalez
        else:
            C[i,2]=C[i,2]+rideHeight-0.05
    geo.restoreCoef(C, 'bodyAxis')

    return

def rideHeight(val,geo):

    C = geo.extractCoef('bodyAxis')
    rideHeight=val[0]
    for i in range(len(C)):
        C[i,2]=C[i,2]+rideHeight-0.05
    geo.restoreCoef(C, 'bodyAxis')

    return


if optVars[0]=='rampAngle':
    DVGeo.addGeoDVGlobal(optVars[0], DVs0, rampAngle,lower=1.0, upper=50.0, scale=1.0)
elif optVars[0]=='rideHeight':
    DVGeo.addGeoDVGlobal(optVars[0], DVs0, rideHeight,lower=0.0, upper=1.0, scale=1.0)
elif optVars[0]=='rampAngleAndRideHeight':
    DVGeo.addGeoDVGlobal(optVars[0], DVs0, rampAngleAndRideHeight,lower=0.0, upper=50.0, scale=1.0)
elif optVars[0]=='shape':
    # Select points
    iVol=2 # iVol=2; ramp of the Ahmed body
    pts=DVGeo.getLocalIndex(iVol)
    indexList=pts[1:,1:,-1].flatten()  # select the top layer FFD starts with i=1
    PS=geo_utils.PointSelect('list',indexList)
    # setup local design variables, lower and upper are the bounds for the FFD points
    DVGeo.addGeoDVLocal(optVars[0], lower=0.0, upper=0.05, axis='z', scale=1.0, pointSelect=PS)


# =================================================================================================
# DAFoam
# =================================================================================================
CFDSolver = PYDAFOAM(options=aeroOptions, comm=gcomm)
CFDSolver.setDVGeo(DVGeo)
mesh = USMesh(options=meshOptions,comm=gcomm)
CFDSolver.addFamilyGroup(CFDSolver.getOption('designsurfacefamily'),CFDSolver.getOption('designsurfaces'))
if MPI.COMM_WORLD.rank == 0:
    CFDSolver.printFamilyList()
CFDSolver.setMesh(mesh)
CFDSolver.computeAdjointColoring()
evalFuncs = CFDSolver.getOption('objfuncs')
xDVs = DVGeo.getValues()
nDVs = 0
for key in xDVs.keys():
    nDVs += len( xDVs[key] )
CFDSolver.setOption('nffdpoints',nDVs)


# =================================================================================================
# DVCon
# =================================================================================================
DVCon = DVConstraints()
DVCon.setDVGeo(DVGeo)
[p0, v1, v2] = CFDSolver.getTriangulatedMeshSurface(groupName=CFDSolver.getOption('designsurfacefamily'))
surf = [p0, v1, v2]
DVCon.setSurface(surf)
#DVCon.writeSurfaceTecplot('trisurface.dat') 


# =================================================================================================
# optFuncs
# =================================================================================================
optFuncs.CFDSolver = CFDSolver
optFuncs.DVGeo = DVGeo
optFuncs.DVCon = DVCon
optFuncs.evalFuncs = evalFuncs
optFuncs.gcomm = gcomm

# =================================================================================================
# Task
# =================================================================================================
if task.lower() == 'run':

    xDV = DVGeo.getValues()
    try:
        for idxI in range(len(xDV[optVars[0]])):
            xDV[optVars[0]][idxI] = float(DVs0[idxI])
    except:
        xDV[optVars[0]] = float(DVs0[0])
    # Evaluate the functions
    funcs = {}
    funcs,fail = optFuncs.aeroFuncs(xDV)

    if gcomm.rank == 0:
        print(funcs)

elif task.lower() == 'writedelmat':

    xDV = DVGeo.getValues()

    try:
        for idxI in range(len(xDV[optVars[0]])):
            xDV[optVars[0]][idxI] = float(DVs0[idxI])
    except:
        xDV[optVars[0]] = float(DVs0[0])


    if gcomm.rank == 0:
        print(("write deltaVolPointsMat at sample=%d"%sample))
        print(xDV)

    DVGeo.setDesignVars(xDV)
    CFDSolver.updateVolumePoints()
    CFDSolver.writeUpdatedVolumePoints()
    CFDSolver._writeDeltaVolPointMat()

elif task.lower() == 'deform':

    xDV = DVGeo.getValues()

    try:
        for idxI in range(len(xDV[optVars[0]])):
            xDV[optVars[0]][idxI] = float(DVs0[idxI])
    except:
        xDV[optVars[0]] = float(DVs0[0])


    if gcomm.rank == 0:
        print(("Deforming at sample=%d"%sample))
        print(xDV)

    DVGeo.setDesignVars(xDV)
    CFDSolver.updateVolumePoints()
    CFDSolver.writeUpdatedVolumePoints()

elif task.lower() == 'xdv2xv':

    optFuncs.xDV2xV()

else:
    print("task arg not found!")
    exit(0)


