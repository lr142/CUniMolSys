import os
import numpy as np

def ListPath():
    li = os.listdir()
    seeds = []
    for s in li:
        if s.lower().find("seed") != -1:
            seeds.append(s)
    seeds.sort()
    li = os.listdir(seeds[0])
    cases = []
    for s in li:
        if s.lower().find("poly") != -1:
            cases.append(s)
    cases.sort()
    combinations = []
    for c in cases:
        combinations.append([])
        for s in seeds:
            combinations[-1].append(os.path.join('.',s,c))
    return combinations

def Density():
    lists = ListPath()
    outfile = open("analyze/RHO.csv","w")
    outfile.write("# Item, Rho, Std\n")
    for l in lists:
        values = np.zeros(len(l))
        for i,c in enumerate(l):
            filename = os.path.join(c,"analyze","RHO.csv")
            with open(filename,"r") as file:
                values[i] = float(file.readlines()[-1])
        mean = np.average(values)
        std = np.std(values)
        outfile.write("{},{},{}\n".format(os.path.basename(c),mean,std))

    outfile.close()

def DiffusionCoeff():
    lists = ListPath()
    filename = os.path.join(lists[0][0],"analyze","MSD.csv")
    with open(filename,"r") as file:
        items_raw = file.readlines()[1].strip().rstrip().split(",")
        items = [i for i in items_raw if i != "#" and i != ""]
    nItems = len(items)
    nSeeds = len(lists[0])
    data = np.zeros((nItems,nSeeds))

    outfile = open("analyze/DIFCOEF.csv", "w")
    outfile.write("# Casename, Items (given in 2nd line), Mean, Std(in pairs)\n")
    outfile.write("#,")
    for item in items:
        outfile.write("{} avg,{} std,".format(item,item))
    outfile.write("\n")
    for l in lists:
        values = np.zeros(len(l))
        for iSeed,c in enumerate(l):
            filename = os.path.join(c,"analyze","MSD.csv")
            # The last 3 lines are slope, intercept, and Diffusion Coefficients. We only need the last line
            lastLine = None
            with open(filename,"r") as file:
                lastLine = file.readlines()[-1].strip()
                if len(lastLine) == 0:
                    lastLine = file.readlines()[-2].strip()
            parts_raw = lastLine.split(",")[1:]
            parts = [float(i) for i in parts_raw if i != ""]
            for index in range(nItems):
                data[index][iSeed] = parts[index]

        outfile.write("{},".format(os.path.basename(c)))
        for iItem in range(nItems):
            mean = np.average(data[iItem])
            std = np.std(data[iItem])
            outfile.write("{:.4f},{:.6f},".format(mean,std))
        outfile.write("\n")
    outfile.close()

if __name__ == '__main__':
    Density()
    DiffusionCoeff()