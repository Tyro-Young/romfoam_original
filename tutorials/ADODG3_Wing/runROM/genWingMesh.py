from pyhyp import pyHyp

fileName = "surfaceMesh.cgns"

options = {
    # ---------------------------
    #        Input Parameters
    # ---------------------------
    "inputFile": fileName,
    "fileType": "cgns",
    "unattachedEdgesAreSymmetry": True,
    "outerFaceBC": "farField",
    "autoConnect": True,
    "BC": {},
    "families": "wall",
    # ---------------------------
    #        Grid Parameters
    # ---------------------------
    "N": 49,
    "s0": 1.0e-3,
    "marchDist": 6,
    #'nConstantStart':1,
    # ---------------------------
    #   Pseudo Grid Parameters
    # ---------------------------
    "ps0": -1,
    "pGridRatio": -1,
    "cMax": 0.3,
    # ---------------------------
    #   Smoothing parameters
    # ---------------------------
    "epsE": 2.0,
    "epsI": 4.0,
    "theta": 2.0,
    "volCoef": 0.20,
    "volBlend": 0.0005,
    "volSmoothIter": 20,
    #'kspreltol':1e-4,
}

hyp = pyHyp(options=options)
hyp.run()
hyp.writePlot3D("volumeMesh.xyz")

