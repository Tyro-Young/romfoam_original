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
import csv

# =============================================================================
# Input Parameters
# =============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("--output", help="Output directory", type=str, default="../output")
parser.add_argument("--task", help="type of run to do", type=str, default="run")
parser.add_argument("--optVars", type=str, help="Vars for the optimizer", default="['shape']")
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

# NOTE: put the ref angle at the end of this list because we will use startFrom latestTime in system/controlDict
if optVars[0] == "twist":
    # 20 samples, 3 DVs twist range [-1:1] degrees on top of the 4 degree ref
    DVs_Train = []
    with open("twistSamples.csv") as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=",")
        for row in csv_reader:
            DVs_Train.append(list(map(float, row)))

    DVs_Predict = []
    with open("twistPredictions.csv") as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=",")
        for row in csv_reader:
            DVs_Predict.append(list(map(float, row)))

elif optVars[0] == "shape":
    DVs_Train = []
    with open("shapeSamples.csv") as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=",")
        for row in csv_reader:
            DVs_Train.append(list(map(float, row)))
    DVs_Predict = []
    with open("shapePredictions.csv") as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=",")
        for row in csv_reader:
            DVs_Predict.append(list(map(float, row)))
else:
    print("optVars not valid")
    exit(1)

DVs_Train = np.asarray(DVs_Train)
DVs_Predict = np.asarray(DVs_Predict)

if gcomm.rank == 0:
    print((sample, optVars))

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
alpha0 = 5.0

evalFuncs = ["CD", "CL"]


def calcUAndDir(UIn, alpha1):
    dragDir = [np.cos(alpha1 * np.pi / 180), np.sin(alpha1 * np.pi / 180), 0]
    liftDir = [-np.sin(alpha1 * np.pi / 180), np.cos(alpha1 * np.pi / 180), 0]
    inletU = [float(UIn * np.cos(alpha1 * np.pi / 180)), float(UIn * np.sin(alpha1 * np.pi / 180)), 0]
    return inletU, dragDir, liftDir


inletu0, dragdir0, liftdir0 = calcUAndDir(UmagIn, alpha0)

# Set the parameters
aeroOptions = {
    # output options
    "casename": "ADODG3_" + task + "_" + optVars[0],
    "outputdirectory": outputDirectory,
    "writesolution": True,
    "usecoloring": False,
    "printalloptions": False,
    # design surfaces and cost functions
    "designsurfacefamily": "designSurfaces",
    "designsurfaces": ["wing"],
    "objfuncs": ["CD", "CL"],
    "objfuncgeoinfo": [["wing"], ["wing"]],
    "referencevalues": {"magURef": UmagIn, "ARef": 0.2754, "LRef": 0.3, "pRef": 0.0, "rhoRef": 1.0},
    "liftdir": liftdir0,
    "dragdir": dragdir0,
    # flow setup
    "adjointsolver": "simpleROMFoam",
    "rasmodel": "SpalartAllmarasFv3",
    "flowcondition": "Incompressible",
    "maxflowiters": args.runEndTime,
    "writeinterval": args.runEndTime,
    "setflowbcs": True,
    "inletpatches": ["inout"],
    "outletpatches": ["inout"],
    "flowbcs": {"bc0": {"patch": "inout", "variable": "U", "value": inletu0}, "useWallFunction": "true"},
    # romDict
    "nsamples": args.nSamples,
    "deltaffd": deltaDVs,
    "svdtype": "cross",
    "svdtol": 1e-8,
    "svdmaxits": 100,
    "svdrequestedn": 20,
    "usemf": 1,
    "mfstep": 1e-6,
    "debugmode": 0,
    "uselspg": 1,
    "romnkgmresmf": 0,
    "romnkabstol": 1e2,
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
    "useRotations": False,
    # point and normal for the symmetry plane
    "symmetryPlanes": [[[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]]],
}

# =================================================================================================
# DVGeo
# =================================================================================================
FFDFile = "./FFD/wingFFD.xyz"
DVGeo = DVGeometry(FFDFile)

# nTwists is the number of FFD points in the spanwise direction
nTwists = DVGeo.addRefAxis("bodyAxis", xFraction=0.25, alignIndex="k")


def twist(val, geo):
    # Set all the twist values
    for i in range(nTwists):
        geo.rot_z["bodyAxis"].coef[i] = -val[i]


# select points
pts = DVGeo.getLocalIndex(0)
indexList = pts[:, :, :].flatten()
PS = geo_utils.PointSelect("list", indexList)
if "shape" in optVars[0]:
    DVGeo.addGeoDVLocal("shape", lower=-1.0, upper=1.0, axis="y", scale=1.0, pointSelect=PS)
if "twist" in optVars[0]:
    twist0 = np.zeros(nTwists)
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
    except Exception:
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
    except Exception:
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
    except Exception:
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
