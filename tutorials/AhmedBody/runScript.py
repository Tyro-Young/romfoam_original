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
from pyoptsparse import Optimization, OPT


# =============================================================================
# Input Parameters
# =============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("--output", help='Output directory', type=str,default='../optOutput/')
parser.add_argument("--opt", help="optimizer to use", type=str, default='slsqp')
parser.add_argument("--task", help="type of run to do", type=str, default='run')
parser.add_argument('--optVars',type=str,help='Vars for the optimizer',default="['rampAngle']")
args = parser.parse_args()
exec('optVars=%s'%args.optVars)
task = args.task
outputDirectory = args.output
gcomm = MPI.COMM_WORLD

rampAngle0=17.0

# Set the parameters for optimization
aeroOptions = {
    # output options
    'casename':                 'AhmedBody_'+task+'_'+optVars[0],
    'outputdirectory':          outputDirectory,
    'writesolution':            True,
    'usecoloring':              True,


    # design surfaces and cost functions 
    'designsurfacefamily':     'designSurfaces', 
    'designsurfaces':          ['body'], 
    'objfuncs':                ['CD'],
    'objfuncgeoinfo':          [['body']],
    'referencevalues':         {'magURef':20.0,'ARef':0.056016,'LRef':1.0,'pRef':0.0,'rhoRef':1.0},
    'liftdir':                 [0.0,0.0,1.0],
    'dragdir':                 [1.0,0.0,0.0],


    # flow setup
    'adjointsolver':           'simpleDAFoam',
    'rasmodel':                'SpalartAllmarasFv3',
    'flowcondition':           'Incompressible',
    'maxflowiters':            500, 
    'writeinterval':           500,
    'setflowbcs':              True,  
    'inletpatches':            ['inlet'],
    'outletpatches':           ['outlet'],
    'flowbcs':                 {'bc0':{'patch':'inlet','variable':'U','value':[20.0,0.0,0.0]},
                                'bc1':{'patch':'outlet','variable':'p','value':[0.0]},
                                'bc2':{'patch':'inlet','variable':'k','value':[0.06]},
                                'bc3':{'patch':'inlet','variable':'omega','value':[400.0]},
                                'bc4':{'patch':'inlet','variable':'epsilon','value':[2.16]},
                                'bc5':{'patch':'inlet','variable':'nuTilda','value':[1.5e-4]},
                                'bc6':{'patch':'inlet','variable':'T','value':[300.0]},
                                'useWallFunction':'true'},                
    'transproperties':         {'nu':1.5E-5,
                                'TRef':300.0,
                                'beta':3e-3,
                                'Pr':0.7,
                                'Prt':0.85}, 


    # adjoint setup
    'adjgmresmaxiters':        500,
    'adjgmresrestart':         500,
    'adjgmresreltol':          1e-6,
    'stateresettol':           1e-3,
    'adjdvtypes':              ['FFD'], 
    'epsderiv':                1.0e-6, 
    'epsderivffd':             1.0e-3,
    'adjpcfilllevel':          0, 
    'adjjacmatordering':       'state',
    'adjjacmatreordering':     'rcm',
    'normalizestates':         [],
    'normalizeresiduals':      [],
    'maxresconlv4jacpcmat':    {'URes':2,'pRes':2,'phiRes':1,'nuTildaRes':2,'kRes':2,'omegaRes':2,'epsilonRes':2},
    'statescaling':            {'UScaling':20.0,
                                'pScaling':200.0,
                                'phiScaling':1.0,
                                'nuTildaScaling':1.5e-4,
                                'kScaling':0.06,
                                'epsilonScaling':2.16,
                                'omegaScaling':400.0},
    
    
    ########## misc setup ##########
    'mpispawnrun':             True,
    'restartopt':              False,

}

# mesh warping parameters, users need to manually specify the symmetry plane
meshOptions = {
    'gridFile':                os.getcwd(),
    'fileType':                'openfoam',
    # point and normal for the symmetry plane
    'symmetryPlanes':          [[[0.,0., 0.],[0., 1., 0.]]], 
}

# options for optimizers
outPrefix = outputDirectory+task+optVars[0]
if args.opt == 'snopt':
    optOptions = {
        'Major feasibility tolerance':  1.0e-7,   # tolerance for constraint
        'Major optimality tolerance':   1.0e-7,   # tolerance for gradient 
        'Minor feasibility tolerance':  1.0e-7,   # tolerance for constraint
        'Verify level':                 -1,
        'Function precision':           1.0e-7,
        'Nonderivative linesearch':     None, 
        'Print file':                   os.path.join(outPrefix+'_SNOPT_print.out'),
        'Summary file':                 os.path.join(outPrefix+'_SNOPT_summary.out')
    }
elif args.opt == 'psqp':
    optOptions = {
        'TOLG':                         1.0e-7,   # tolerance for gradient 
        'TOLC':                         1.0e-7,   # tolerance for constraint
        'MIT':                          25,       # max optimization iterations
        'IFILE':                        os.path.join(outPrefix+'_PSQP.out')
    }
elif args.opt == 'slsqp':
    optOptions = {
        'ACC':                          1.0e-7,   # convergence accuracy
        'MAXIT':                        25,       # max optimization iterations
        'IFILE':                        os.path.join(outPrefix+'_SLSQP.out')
    }
else:
    print("opt arg not valid!")
    exit(0)


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

DVGeo.addGeoDVGlobal('rampAngle', rampAngle0, rampAngle,lower=5.0, upper=50.0, scale=1.0)

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
if task.lower()=='opt':
    optProb = Optimization('opt', optFuncs.aeroFuncs, comm=gcomm)
    DVGeo.addVariablesPyOpt(optProb)
    DVCon.addConstraintsPyOpt(optProb)

    # Add objective
    optProb.addObj('CD', scale=1)
    # Add physical constraints
    #optProb.addCon('CL',lower=0.5,upper=0.5,scale=1)

    if gcomm.rank == 0:
        print optProb

    opt = OPT(args.opt, options=optOptions)
    histFile = os.path.join(outputDirectory, '%s_hist.hst'%args.opt)
    sol = opt(optProb, sens=optFuncs.aeroFuncsSens, storeHistory=histFile)
    if gcomm.rank == 0:
        print sol

elif task.lower() == 'run':

    xDV = DVGeo.getValues()

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

elif task.lower() == 'writedeltavolmat':

    CFDSolver._writeDeltaVolPointMat()

elif task.lower() == 'testsensuin':

    optFuncs.testSensUIn(normStatesList=[True],deltaUList=[1e-8])
        
elif task.lower() == 'testsensshape':

    optFuncs.testSensShape(normStatesList=[True],deltaUList=[1e-7],deltaXList=[1e-4])

elif task.lower() == 'xdv2xv':

    optFuncs.xDV2xV()

else:
    print("task arg not found!")
    exit(0)


