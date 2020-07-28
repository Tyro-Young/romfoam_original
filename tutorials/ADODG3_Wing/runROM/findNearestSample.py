import sys
import csv
import numpy as np

def findNearestSample(sampleFileName, predictFileName, predictIndex):

    DVs_Train = []
    with open(sampleFileName) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=",")
        for row in csv_reader:
            DVs_Train.append(list(map(float, row)))

    DVs_Predict = []
    with open(predictFileName) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=",")
        for row in csv_reader:
            DVs_Predict.append(list(map(float, row)))

    predict = DVs_Predict[predictIndex]

    DVs_Train = np.asarray(DVs_Train)
    predict = np.asarray(predict)

    minNorm = 1e10
    minNormIdx = -999
    for idxI, train in enumerate(DVs_Train):
        diffNorm = np.linalg.norm(train - predict)
        if diffNorm < minNorm:
            minNorm = diffNorm
            minNormIdx = idxI

    #print("Finding the nearest sample for initial conditions...\n")
    #print("The prediction point is: ", predict, "\n")
    # sample starts with 1 while the python list starts with 0
    #print("The nearest sample is sample%d \n" % (minNormIdx + 1))
    #print("Its design variables are: ", DVs_Train[minNormIdx], "\n")

    return str(minNormIdx + 1)


if __name__ == "__main__":

    sampleFileName = sys.argv[1]
    predictFileName = sys.argv[2]
    predictIndex = int(sys.argv[3])

    returnValue = findNearestSample(sampleFileName, predictFileName, predictIndex)

    sys.exit(returnValue)
