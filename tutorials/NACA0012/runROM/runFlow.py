#!/usr/bin/env python
"""
ROMFoam run script for the NACA0012 airfoil at low-speed
"""

# =================================================================================================
# Imports
# =================================================================================================
import os
import argparse
import numpy as np
from mpi4py import MPI
from romfoam import *
from pygeo import *
from pyspline import *
from idwarp import *


# =============================================================================
# Input Parameters
# =============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("--output", help="Output directory", type=str, default="../output")
parser.add_argument("--task", help="type of run to do", type=str, default="run")
parser.add_argument("--optVars", type=str, help="Vars for the optimizer", default="['twist']")
parser.add_argument("--sample", type=int, help="which sample DV to run", default=1)
parser.add_argument("--mode", type=str, help="can be either train or predict", default="train")
parser.add_argument("--runEndTime", type=int, help="number of time steps for flow", default=500)
parser.add_argument("--avgFieldEvery", type=int, help="average obj and field every ? steps", default=-1)
parser.add_argument("--nSamples", type=int, help="number of samples", default=1)
args = parser.parse_args()
exec("optVars=%s" % args.optVars)
task = args.task
outputDirectory = args.output
gcomm = MPI.COMM_WORLD

sample = args.sample - 1

if gcomm.rank == 0:
    print((sample, optVars))

# NOTE: put the ref angle at the end of this list because we will use startFrom latestTime in system/controlDict
if optVars[0] == "twist":
    # 5 samples, 1 DV
    DVs_Train = [[1.0], [3.0], [7.0], [9.0], [5.0]]
    DVs_Predict = [[4.0], [6.0]]
else:
    print("optVars not valid")
    exit(1)

DVs_Train = np.asarray(DVs_Train)
DVs_Predict = np.asarray(DVs_Predict)

if args.mode == "train":
    DVs0 = DVs_Train[sample]
elif args.mode == "predict":
    DVs0 = DVs_Predict[sample]
else:
    print("mode should be either train or predict")
    exit(1)

if args.mode == "predict":
    deltaDVs = DVs_Predict[sample] - DVs_Train[-1]
    deltaDVs = deltaDVs.tolist()
elif args.mode == "train":
    deltaDVs = [0.0]

UmagIn = 10.0
CL_star = 0.5
twist0 = 4.193339

evalFuncs = ["CD"]

# Set the parameters
aeroOptions = {
    # output options
    "casename": "NACA0012_" + task + "_" + optVars[0],
    "outputdirectory": outputDirectory,
    "writesolution": True,
    "usecoloring": False,
    "printalloptions": False,
    # design surfaces and cost functions
    "designsurfacefamily": "designSurfaces",
    "designsurfaces": ["wing"],
    "objfuncs": ["CD", "CL"],
    "objfuncgeoinfo": [["wing"], ["wing"]],
    "referencevalues": {"magURef": UmagIn, "ARef": 0.1, "LRef": 1.0, "pRef": 0.0, "rhoRef": 1.0},
    "liftdir": [0.0, 1.0, 0.0],
    "dragdir": [1.0, 0.0, 0.0],
    # flow setup
    "adjointsolver": "simpleROMFoam",
    "rasmodel": "SpalartAllmarasFv3",
    "flowcondition": "Incompressible",
    "maxflowiters": args.runEndTime,
    "writeinterval": args.runEndTime,
    "setflowbcs": True,
    "inletpatches": ["inout"],
    "outletpatches": ["inout"],
    "flowbcs": {"bc0": {"patch": "inout", "variable": "U", "value": [UmagIn, 0.0, 0.0]}, "useWallFunction": "true"},
    # romDict
    "nsamples": args.nSamples,
    "deltaffd": deltaDVs,
    "svdtype": "cross",
    "svdtol": 1e-8,
    "svdmaxits": 100,
    "svdrequestedn": args.nSamples,
    "usemf": 1,
    "mfstep": 1e-6,
    "debugmode": 0,
    "uselspg": 1,
    "romnkgmresmf": 0,
    "romnkabstol": 1e-1,
    "romnkmaxits": 20,
    "romnkgmresmaxls": 5,
    # adjoint setup
    "adjdvtypes": ["FFD"],
    "epsderivffd": 1.0e-3,
    "adjjacmatordering": "state",
    "adjjacmatreordering": "rcm",
    "normalizestates": [],
    "normalizeresiduals": ["URes", "pRes", "nuTildaRes"],
    ########## misc setup ##########
    "mpispawnrun": False,
}

# mesh warping parameters, users need to manually specify the symmetry plane
meshOptions = {
    "gridFile": os.getcwd(),
    "fileType": "openfoam",
    'useRotations': False,
    # point and normal for the symmetry plane
    "symmetryPlanes": [[[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]], [[0.0, 0.0, 0.1], [0.0, 0.0, 1.0]]],
}

# =================================================================================================
# DVGeo
# =================================================================================================
FFDFile = "./FFD/wingFFD.xyz"
DVGeo = DVGeometry(FFDFile)

# ref axis
x = [0.25, 0.25]
y = [0.00, 0.00]
z = [0.00, 0.10]
c1 = pySpline.Curve(x=x, y=y, z=z, k=2)
DVGeo.addRefAxis("bodyAxis", curve=c1, axis="z")


def twist(val, geo):
    # Set all the twist values
    for i in range(2):
        geo.rot_z["bodyAxis"].coef[i] = -val[0]


# select points
pts = DVGeo.getLocalIndex(0)
indexList = pts[:, :, :].flatten()
PS = geo_utils.PointSelect("list", indexList)
# DVGeo.addGeoDVLocal("shapey", lower=-1.0, upper=1.0, axis="y", scale=1.0, pointSelect=PS)
DVGeo.addGeoDVGlobal("twist", twist0, twist, lower=-10.0, upper=10.0, scale=1.0)

# =================================================================================================
# DAFoam
# =================================================================================================
CFDSolver = PYDAFOAM(options=aeroOptions, comm=gcomm)
CFDSolver.setDVGeo(DVGeo)
mesh = USMesh(options=meshOptions, comm=gcomm)
CFDSolver.addFamilyGroup(CFDSolver.getOption("designsurfacefamily"), CFDSolver.getOption("designsurfaces"))
if MPI.COMM_WORLD.rank == 0:
    CFDSolver.printFamilyList()
CFDSolver.setMesh(mesh)
CFDSolver.computeAdjointColoring()
evalFuncs = CFDSolver.getOption("objfuncs")
xDVs = DVGeo.getValues()
nDVs = 0
for key in xDVs.keys():
    nDVs += len(xDVs[key])
CFDSolver.setOption("nffdpoints", nDVs)

# =================================================================================================
# DVCon
# =================================================================================================
DVCon = DVConstraints()
DVCon.setDVGeo(DVGeo)
[p0, v1, v2] = CFDSolver.getTriangulatedMeshSurface(groupName=CFDSolver.getOption("designsurfacefamily"))
surf = [p0, v1, v2]
DVCon.setSurface(surf)
# DVCon.writeSurfaceTecplot('trisurface.dat')

# =================================================================================================
# Task
# =================================================================================================
if task.lower() == "run":

    xDV = DVGeo.getValues()
    try:
        for idxI in range(len(xDV[optVars[0]])):
            xDV[optVars[0]][idxI] = float(DVs0[idxI])
    except:
        xDV[optVars[0]] = float(DVs0[0])
    # Evaluate the functions
    funcs = {}
    funcs, fail = optFuncs.aeroFuncs(xDV)

    if gcomm.rank == 0:
        print(funcs)

elif task.lower() == "writedelmat":

    xDV = DVGeo.getValues()

    try:
        for idxI in range(len(xDV[optVars[0]])):
            xDV[optVars[0]][idxI] = float(DVs0[idxI])
    except:
        xDV[optVars[0]] = float(DVs0[0])

    if gcomm.rank == 0:
        print(("write deltaVolPointsMat at sample=%d" % sample))
        print(xDV)

    DVGeo.setDesignVars(xDV)
    CFDSolver.updateVolumePoints()
    CFDSolver.writeUpdatedVolumePoints()
    CFDSolver._writeDeltaVolPointMat()

elif task.lower() == "deform":

    xDV = DVGeo.getValues()

    try:
        for idxI in range(len(xDV[optVars[0]])):
            xDV[optVars[0]][idxI] = float(DVs0[idxI])
    except:
        xDV[optVars[0]] = float(DVs0[0])

    if gcomm.rank == 0:
        print(("Deforming at sample=%d" % sample))
        print(xDV)

    DVGeo.setDesignVars(xDV)
    CFDSolver.updateVolumePoints()
    CFDSolver.writeUpdatedVolumePoints()

elif task.lower() == "xdv2xv":

    optFuncs.xDV2xV()

else:
    print("task arg not found!")
    exit(0)
