#!/usr/bin/env python
"""
ROMFoam run script for the Ahmed body case
"""

# =================================================================================================
# Imports
# =================================================================================================
import os, time
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
parser.add_argument("--output", help="Output directory", type=str, default="../output")
parser.add_argument("--task", help="type of run to do", type=str, default="run")
parser.add_argument("--optVars", type=str, help="Vars for the optimizer", default="['shape']")
parser.add_argument("--sample", type=int, help="which sample DV to run", default=1)
parser.add_argument("--mode", type=str, help="can be either train or predict", default="train")
parser.add_argument("--runEndTime", type=int, help="number of time steps for flow", default=500)
parser.add_argument("--avgFieldEvery", type=int, help="average obj and field every ? steps", default=-1)
parser.add_argument("--nSamples", type=int, help="number of samples", default=1)
parser.add_argument("--resSamples", type=int, help="number of res samples", default=23)
parser.add_argument("--nSampleSV", type=int, help="number of sample SVs", default=5)
args = parser.parse_args()
exec("optVars=%s" % args.optVars)
task = args.task
outputDirectory = args.output
gcomm = MPI.COMM_WORLD

sample = args.sample - 1

if args.avgFieldEvery > 0:
    useAvg = True
    avgFrom = args.runEndTime - args.avgFieldEvery
    avgFields = ["U", "p", "phi"]
else:
    useAvg = False
    avgFrom = 0
    avgFields = []

if gcomm.rank == 0:
    print((sample, optVars))

# NOTE: put the ref angle at the end of this list because we will use startFrom latestTime in system/controlDict
if optVars[0] == "rampAngle":
    # 5 samples, 1 DV
    DVs_Train = [[5.0], [10.0], [15.0], [20.0], [25.0]]
    DVs_Predict = [[9.0],[13.0],[17.0],[21.0]]
elif optVars[0] == "rideHeight":
    # 5 samples, 1 DV
    DVs_Train = [[0.03], [0.04], [0.06], [0.07], [0.05]]
    DVs_Predict = [[0.07], [0.03]]
elif optVars[0] == "rampAngleAndRideHeight":
    # 10 samples, 2 DVs
      DVs_Train= [[20.0,0.045],
                  [25.0,0.045],
                  [30.0,0.045],
                  [35.0,0.045],
                  [40.0,0.045],
                  [20.0,0.050],
                  [25.0,0.050],
                  [30.0,0.050],
                  [35.0,0.050],
                  [40.0,0.050],
                  [20.0,0.055],
                  [25.0,0.055],
                  [30.0,0.055],
                  [35.0,0.055],
                  [40.0,0.055],
                  [23.5,0.052],
                  [28.5,0.048],
                  [33.5,0.054],
                  [38.5,0.046]]
      DVs_Predict = [[25.7, 0.043], [24.1, 0.045]]

elif optVars[0] == "shape":
    DVs_Train = [[0.019546,0.025956,0.018573],
    [0.036348,0.019889,0.011424],
    [0.034073,3.8459e-05,0.032598],
    [0.012593,0.029996,0.030649],
    [0.00076034,0.038672,0.026899],
    [0.027253,0.035304,0.015976],
    [0.031011,0.0054975,0.00079235],
    [0.0055642,0.0089961,0.037849],
    [0.020891,0.023386,0.0045163],
    [-0.0027643,-0.050315,-0.068595],
    [-0.029095,-0.031552,-0.060527],
    [-0.050164,-0.064028,-0.053321],
    [-0.069289,-0.069093,-0.074236],
    [-0.06595,-0.029963,-2.4963e-05],
    [-0.059183,-0.00071484,-0.029649],
    [-0.03428,-0.046467,-0.03506],
    [-0.053021,-0.05445,-0.044193],
    [-0.073935,-0.025243,-0.010335],
    [-0.038942,-0.040282,-0.058726],
    [-0.023327,-0.020696,-0.006758],
    [-0.060169,-0.074147,-0.039349],
    [-0.015301,-0.037318,-0.063935],
    [-0.013557,-0.044365,-0.03227],
    [-0.03314,-0.017781,-0.046558],
    [-0.045762,-0.010508,-0.016475],
    [-0.044573,-0.0052287,-0.04982],
    [-0.0069044,-0.060661,-0.024744],
    [-0.0099156,-0.05696,-0.01178],
    [-0.019321,-0.011992,-0.019048],
    [-0.075,-0.075,-0.075],
    [-0.07,-0.07,-0.07],
    [-0.065,-0.065,-0.065],
    [0.03,0.03,0.03],
    [0.035,0.035,0.035],
    [0.04,0.04,0.04],
    [-0.06,-0.06,-0.06],
    [-0.05,-0.05,-0.05],
    [-0.04,-0.04,-0.04],
    [-0.03,-0.03,-0.03],
    [-0.02,-0.02,-0.02],
    [-0.01,-0.01,-0.01],
    [0,0,0],
    [0.01,0.01,0.01],
    [0.02,0.02,0.02],
    [0.037539,0.02741,0.030662],
    [0.014791,0.020132,0.037229],
    [0.0027676,0.015129,0.014032],
    [0.024251,0.035776,0.0034007],
    [-0.065693,-0.047967,-0.053659],
    [-0.025885,-0.03523,-0.065152],
    [-0.034088,-0.0038042,-0.034091],
    [-0.0048433,-0.026475,-0.024555],
    [-0.042438,-0.062607,-0.0059511]]

    DVs_Predict = [[0.04, 0.04, 0.04]]

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


# Set the parameters for optimization
aeroOptions = {
    # output options
    "casename": "AhmedBody_" + task + "_" + optVars[0],
    "outputdirectory": outputDirectory,
    "writesolution": True,
    "usecoloring": False,
    "printalloptions": False,
    # design surfaces and cost functions
    "designsurfacefamily": "designSurfaces",
    "designsurfaces": ["body"],
    "objfuncs": ["CD"],
    "objfuncgeoinfo": [["body"]],
    "referencevalues": {"magURef": 20.0, "ARef": 0.056016, "rhoRef": 1.0, "pRef": 0.0, "LRef": 1.0},
    # flow setup
    "adjointsolver": "simpleROMFoam",
    "rasmodel": "SpalartAllmarasFv3",
    "flowcondition": "Incompressible",
    "maxflowiters": args.runEndTime,
    "writeinterval": 10,
    "avgobjfuncs": useAvg,
    "avgobjfuncsstart": avgFrom,
    "avgfields": avgFields,
    "avgfieldsevery": args.avgFieldEvery,
    "setflowbcs": False,
    "inletpatches": ["inlet"],
    "outletpatches": ["outlet"],
    "flowbcs": {"bc0": {"patch": "inlet", "variable": "U", "value": [20.0, 0.0, 0.0]}, "useWallFunction": "true"},
    # romDict
    "nsamples": args.nSamples,
    "ressamples": args.resSamples,
    "nsamplesv": args.nSampleSV,
    "deltaffd": deltaDVs,
    "svdtype": "cross",
    "svdtol": 1e-8,
    "svdmaxits": 1000,
    "svdrequestedn": 50,
    "usemf": 1,
    "romusesvdres": 1,
    "mfstep": 1e-6,
    "debugmode": 0,
    "uselspg": 0,
    "romnkabstol": 1e-8,
    "romnkmaxits": 200,
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
    # point and normal for the symmetry plane
    "symmetryPlanes": [[[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]],
}

# =================================================================================================
# DVGeo
# =================================================================================================
FFDFile = "./FFD/bodyFittedFFD.xyz"
DVGeo = DVGeometry(FFDFile)
x = [0.000, 0.100, 0.843, 1.044]
y = [0.100, 0.100, 0.100, 0.100]
z = [0.194, 0.194, 0.194, 0.147]
c1 = pySpline.Curve(x=x, y=y, z=z, k=2)
DVGeo.addRefAxis("bodyAxis", curve=c1, axis="z")


def rampAngle(val, geo):

    C = geo.extractCoef("bodyAxis")

    # the value will be ramp angle in degree.
    # start with a conversion to rads
    angle = (val[0]) * np.pi / 180.0

    # Overall length needs to stay a 1.044, so use that as a ref for
    # the final mesh point

    # set the target length
    lTarget = 0.222
    hInit = 0.246 - 0.05

    # compute the coefficient deltas
    dx = lTarget * np.cos(angle)
    dz = lTarget * np.sin(angle)

    topEdge = 0.338 - dz
    rearHeight = topEdge - 0.05
    coefPoint = rearHeight / 2.0 + 0.05
    scalez = rearHeight / hInit

    # Set the coefficients
    C[3, 0] = 1.044
    C[2, 0] = C[3, 0] - dx
    C[2, 2] = 0.194
    C[3, 2] = coefPoint

    geo.restoreCoef(C, "bodyAxis")

    geo.scale_z["bodyAxis"].coef[3] = scalez

    return


def rampAngleAndRideHeight(val, geo):

    C = geo.extractCoef("bodyAxis")

    # the value will be ramp angle in degree.
    # start with a conversion to rads
    angle = (val[0]) * np.pi / 180.0

    # Overall length needs to stay a 1.044, so use that as a ref for
    # the final mesh point

    # set the target length
    lTarget = 0.222
    hInit = 0.246 - 0.05

    # compute the coefficient deltas
    dx = lTarget * np.cos(angle)
    dz = lTarget * np.sin(angle)

    topEdge = 0.338 - dz
    rearHeight = topEdge - 0.05
    coefPoint = rearHeight / 2.0 + 0.05
    scalez = rearHeight / hInit

    # Set the coefficients
    C[3, 0] = 1.044
    C[2, 0] = C[3, 0] - dx
    C[2, 2] = 0.194
    C[3, 2] = coefPoint

    geo.restoreCoef(C, "bodyAxis")

    geo.scale_z["bodyAxis"].coef[3] = scalez

    rideHeight = val[1]
    for i in range(len(C)):
        if i >= 2:
            C[i, 2] = C[i, 2] + (rideHeight - 0.05) * scalez
        else:
            C[i, 2] = C[i, 2] + rideHeight - 0.05
    geo.restoreCoef(C, "bodyAxis")

    return


def rideHeight(val, geo):

    C = geo.extractCoef("bodyAxis")
    rideHeight = val[0]
    for i in range(len(C)):
        C[i, 2] = C[i, 2] + rideHeight - 0.05
    geo.restoreCoef(C, "bodyAxis")

    return


if optVars[0] == "rampAngle":
    DVGeo.addGeoDVGlobal(optVars[0], DVs0, rampAngle, lower=1.0, upper=50.0, scale=1.0)
elif optVars[0] == "rideHeight":
    DVGeo.addGeoDVGlobal(optVars[0], DVs0, rideHeight, lower=0.0, upper=1.0, scale=1.0)
elif optVars[0] == "rampAngleAndRideHeight":
    DVGeo.addGeoDVGlobal(optVars[0], DVs0, rampAngleAndRideHeight, lower=0.0, upper=50.0, scale=1.0)
elif optVars[0] == "shape":
    # Select points
    iVol = 2  # iVol=2; ramp of the Ahmed body
    pts = DVGeo.getLocalIndex(iVol)
    print("printing points!")
    print(pts)
    indexList = pts[1, 0:, -1].flatten()  # select the top layer FFD starts with i=1
    print(indexList)
    PS = geo_utils.PointSelect("list", indexList)
    print(PS)
    # setup local design variables, lower and upper are the bounds for the FFD points
    DVGeo.addGeoDVLocal(optVars[0], lower=-0.1, upper=0.05, axis="z", scale=1.0, pointSelect=PS)


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
