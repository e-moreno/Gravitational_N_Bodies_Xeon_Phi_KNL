#!/usr/bin/python3

def myModel(N, iterations, time):
    # 16N² + 18N
    return (((16 * N * N) + 18 * N)*iterations*1E-9)/time

def standardModel(N, iterations, time):
    # 20N²
    return ((20*N*N)*iterations*1E-9)/time

filename = input()
dataFile = open(filename, "r")

simulationParams = dataFile.readline().split(",")
N = int(simulationParams[0])
iterations = int(simulationParams[1])

for line in dataFile:
    time = float(line.strip())
    print(str(myModel(N, iterations, time)) + ", " + str(standardModel(N,
                                                         iterations, time)))
