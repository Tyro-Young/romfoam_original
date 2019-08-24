#!/usr/bin/env python
"""
DAFoam run script for the Ahmed body case
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
from dafoam import *
from pygeo import *
from pyspline import *
from idwarp import *


# =============================================================================
# Input Parameters
# =============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("--output", help='Output directory', type=str,default='../output')
parser.add_argument("--opt", help="optimizer to use", type=str, default='slsqp')
parser.add_argument("--task", help="type of run to do", type=str, default='run')
parser.add_argument('--optVars',type=str,help='Vars for the optimizer',default="['rampAngle']")
parser.add_argument('--sample',type=int,help='which sample DV to run',default=1)
parser.add_argument('--mode',type=str,help='can be either train or predict',default='train')
parser.add_argument('--nSamples',type=int,help='number of samples',default=1)
args = parser.parse_args()
exec('optVars=%s'%args.optVars)
task = args.task
outputDirectory = args.output
gcomm = MPI.COMM_WORLD

sample=args.sample-1

if gcomm.rank==0:
    print(sample,optVars)

DVs_Train=[[5.0],[25.0],[15.0]]
DVs_Predict=[[5.0],[25.0],[15.0]]
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
    if gcomm.rank==0:
        f=open('system/romDict','w')
        f.write('FoamFile{version 2.0;format ascii;class dictionary;location %s;object adjointDict;}\n'%"system")
        f.write('nSamples %d; deltaFFD (%s); svdType cross; svdTol 1e-8; svdMaxIts 100; svdRequestedN %d; useMF 1; mfStep 1e-6;\n'%(args.nSamples,' '.join(map(str,deltaDVs)),args.nSamples))
elif args.mode=='train': 
    if gcomm.rank==0:
        f=open('system/romDict','w')
        f.write('FoamFile{version 2.0;format ascii;class dictionary;location %s;object adjointDict;}\n'%"system")
        f.write('nSamples %d; deltaFFD (1.0); svdType cross; svdTol 1e-8; svdMaxIts 100; svdRequestedN %d; useMF 1; mfStep 1e-6;\n'%(args.nSamples,args.nSamples))


# Set the parameters for optimization
aeroOptions = {
    # output options
    'casename':                 'AhmedBody_'+task+'_'+optVars[0],
    'outputdirectory':          outputDirectory,
    'writesolution':            True,
    'usecoloring':              False,


    # design surfaces and cost functions 
    'designsurfacefamily':     'designSurfaces', 
    'designsurfaces':          ['body'], 
    'objfuncs':                ['CD'],
    'objfuncgeoinfo':          [['body']],
    'referencevalues':         {'magURef':20.0,'ARef':0.056016,'rhoRef':1.0,'pRef':0.0,'LRef':1.0},

    # flow setup
    'adjointsolver':           'simpleDAFoam',
    'rasmodel':                'SpalartAllmarasFv3',
    'flowcondition':           'Incompressible',
    'maxflowiters':            500, 
    'writeinterval':           500,
    'avgobjfuncs':             False,
    'avgobjfuncsstart':        2000,
    'setflowbcs':              False, 
    'inletpatches':            ['inlet'],
    'outletpatches':           ['outlet'],
    'flowbcs':                 {'bc0':{'patch':'inlet','variable':'U','value':[20.0,0.0,0.0]},
                                'bc1':{'patch':'outlet','variable':'p','value':[0.0]},
                                'bc2':{'patch':'inlet','variable':'k','value':[0.06]},
                                'bc3':{'patch':'inlet','variable':'omega','value':[400.0]},
                                'bc4':{'patch':'inlet','variable':'epsilon','value':[2.16]},
                                'bc5':{'patch':'inlet','variable':'nuTilda','value':[1.5e-4]},
                                'useWallFunction':'true'},                
    'transproperties':         {'nu':1.5E-5,
                                'TRef':300.0,
                                'beta':3e-3,
                                'Pr':0.7,
                                'Prt':0.85}, 

    # adjoint setup
    'adjdvtypes':              ['FFD'], 
    'epsderivffd':             1.0e-3,
    'adjjacmatordering':       'state',
    'adjjacmatreordering':     'rcm',
    'normalizestates':         [],
    'normalizeresiduals':      [],    
    
    ########## misc setup ##########
    'mpispawnrun':             False,
    'restartopt':              False,

}

# mesh warping parameters, users need to manually specify the symmetry plane
meshOptions = {
    'gridFile':                os.getcwd(),
    'fileType':                'openfoam',
    # point and normal for the symmetry plane
    'symmetryPlanes':          [[[0.,0., 0.],[0., 1., 0.]]], 
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

if 'rampAngle' in optVars:
    DVGeo.addGeoDVGlobal(optVars[0], 25.0, rampAngle,lower=1.0, upper=50.0, scale=1.0)


if 'shape' in optVars:
    # Select points
    iVol=2 # iVol=2; ramp of the Ahmed body
    pts=DVGeo.getLocalIndex(iVol)
    indexList=pts[1:,:,-1].flatten()  # select the top layer FFD starts with i=1
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
   
    if isinstance(xDV[optVars[0]],list):
        for idxI in range(len(xDV[optVars[0]])):
            xDV[optVars[0]][idxI] = float(DVs0[idxI])
    else:
        xDV[optVars[0]] = float(DVs0[0])

    # Evaluate the functions
    funcs = {}
    funcs,fail = optFuncs.aeroFuncs(xDV)

    if gcomm.rank == 0:
        print funcs
    
    # Evaluate the sensitivities
    #funcsSens = {}
    #funcsSens,fail = optFuncs.aeroFuncsSens(xDV,funcs)
    
    #if gcomm.rank == 0:
    #    print funcsSens

elif task.lower() == 'writedelmat':

    xDV = DVGeo.getValues()

    if isinstance(xDV[optVars[0]],list):
        for idxI in range(len(xDV[optVars[0]])):
            xDV[optVars[0]][idxI] = float(DVs0[idxI])
    else:
        xDV[optVars[0]] = float(DVs0[0])


    if gcomm.rank == 0:
        print ("write deltaVolPointsMat at sample=%d"%sample)
        print(xDV)

    DVGeo.setDesignVars(xDV)
    CFDSolver.updateVolumePoints()
    CFDSolver.writeUpdatedVolumePoints()
    CFDSolver._writeDeltaVolPointMat()

elif task.lower() == 'deform':

    xDV = DVGeo.getValues()

    if isinstance(xDV[optVars[0]],list):
        for idxI in range(len(xDV[optVars[0]])):
            xDV[optVars[0]][idxI] = float(DVs0[idxI])
    else:
        xDV[optVars[0]] = float(DVs0[0])


    if gcomm.rank == 0:
        print ("Deforming at sample=%d"%sample)
        print(xDV)

    DVGeo.setDesignVars(xDV)
    CFDSolver.updateVolumePoints()
    CFDSolver.writeUpdatedVolumePoints()

elif task.lower() == 'xdv2xv':

    optFuncs.xDV2xV()

else:
    print("task arg not found!")
    exit(0)


