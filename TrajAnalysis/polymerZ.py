import sys
import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

def Smooth(array, strength=0.9):
    result = np.zeros_like(array)
    result[0] = array[0]
    for iVal in range(1,len(array)):
        result[iVal] = result[iVal-1]*strength + array[iVal] * (1-strength)
    return result


# This class created an index of a whole folder, can return the correct path given a series of keywords
# note that all keyword-based searching are case-insensitive
class PathFinder:
    # pwd is the working dir (default ./), if key is present, only folders contains this key will be indexed
    # currently we use the keyword 'analyze' (subfolder under each folder) to store all analysis data
    def __init__(self,pwd="./",key='analyze'):
        self.pwd = pwd;
        self.files = []
        for dir in os.walk(pwd,followlinks=True):
            foldername = dir[0]
            if key == "" or foldername.upper().find(key.upper())!= -1:
                for file in dir[-1]:
                    # index all files (not subfolders) under this folder
                    file = os.path.join(dir[0],file)
                    self.files.append(file)

    def Show(self):
        for d in self.files:
            print(d)
    def Search(self,keywords:[str]):
        for file in self.files:
            match = True
            # to find a match, all keywords must be present
            for key in keywords:
                if file.upper().find(key.upper()) == -1:
                    match = False
                    break
            # if there are multiple matches, only the first match will be returned!
            if match:
                return file
        # if reaches here, there is no match
        print("Not found : {}".format(keywords))
        return ""
def DataFileReader(filename:str,transpose:bool=True):
    if(filename==""):
        print("filename is empty, re-check!")
        return
    else:
        print("Reading [{}]...".format(filename))

    file = None
    try:
        file = open(filename,"r",encoding='UTF-8')
    except:
        print("Can't open {}".format(file))
    # first few lines with # are comments
    line = ""
    while(True):
        line = file.readline()
        if len(line.strip()) == 0 or line.startswith("#"):
            continue
        else:
            break
    parts = line.strip().split(',')
    nRow = int(parts[0])
    nCol = int(parts[1])
    data = np.zeros((nRow,nCol))
    for i in range(nRow):
        parts = file.readline().strip().split(',')
        for j in range(nCol):
            data[i,j] = float(parts[j])
    if(transpose):
        data = data.transpose()
    return data

class PlotSeries:
    def __init__(self,x,y,label,c='k',ls='-',lw=1,m=None,ms=0):
        self.x = x
        self.y = y
        self.color = c
        self.linewidth = lw
        self.linestyle = ls
        self.markersize = ms
        self.marker = m
        self.label = label

class Plotter:
    f: plt.Figure
    ax: [plt.Axes] # ax[0] is a bar chart, ylabel on the left. ax[1] is a line plot, ylabel on the right
    outputfile: str
    dpi:int
    series:[[PlotSeries],[PlotSeries]]

    def __init__(self,outputfile,figsize=(4,3),dpi=600,bar_alpha=0.6):
        mpl.rcParams["font.sans-serif"] = ["DejaVu Sans"]
        mpl.rcParams["figure.max_open_warning"] = False
        self.f = plt.figure(figsize=figsize)
        self.dpi = dpi
        self.ax = [None,None]
        self.ax[0] = self.f.add_subplot(111)
        self.ax[1] = self.ax[0].twinx()
        self.outputfile = outputfile
        self.series = [[],[]]
        self.bar_alpha = bar_alpha
    def AddSeries(self,series:PlotSeries,axId):
        self.series[axId].append(series)
    def Plot(self):
        for i in range(1):
            for j,s in enumerate(self.series[i]):
                NSeries = len(self.series[i])  # number of series, also number of bars cluttered together, eg. 3
                NSeriesMid = (NSeries-1)/2.0  # the mid number of series, eg. if NSeries = 3, this number is 1. This # is used to offset series to the center
                width = 0.8*s.x[0]*2/NSeries
                self.ax[i].bar(s.x+width*(j-NSeriesMid),s.y,label=s.label,
                               width=width,alpha=self.bar_alpha,edgecolor=None,facecolor=s.color,fill=True)
        for i in range(1,2):
            for s in self.series[i]:
                self.ax[i].plot(s.x,s.y,color=s.color,
                             linestyle=s.linestyle,linewidth=s.linewidth,
                             markersize=s.markersize,marker=s.marker,
                             label=s.label)

    def SetAxes(self,xmin,xmax,ymin1,ymax1,xlabel,ylabel1,ymin2=None,ymax2=None,ylabel2=None):
        self.ax[0].set_xlim(xmin, xmax)
        self.ax[0].set_ylim(ymin1, ymax1)
        self.ax[0].set_xlabel(xlabel)
        self.ax[0].set_ylabel(ylabel1)
        if ymin2 != None:
            self.ax[1].tick_params('y')
            self.ax[1].set_ylim(ymin2,ymax2)
        if ylabel2 != None:
            self.ax[1].set_ylabel(ylabel2)
    def SetLegend(self,whichaxis=0,loc='lower right'):
        self.ax[whichaxis].legend(loc=loc)
    def SetTitle(self,title:str):
        self.ax[0].set_title(title)

    def __del__(self):
        plt.savefig(self.outputfile,dpi=self.dpi)

# title is the string shown on the figure
# prefix is prepend to the filename, such as a numbering
# seeds means which seed are considered
# fixed are keywords same for all series
# vars are keywords that are different for each series
# plotSchemes are the colors, linestyles, etc, for each series
# if vars is empty [], it means plot each seed separately.
def PolymerZSingleCase(prefix:str, title:str, pathFinder:PathFinder, seeds:[str], fixed:[str], vars:[str], plotSchemes:[PlotSeries]):
    # the name of the output jpg file. Split and Rejoin to replace ' ' with '_'
    filename = os.path.join(pathFinder.pwd,"analyze","{}{}.jpg".format(prefix,"_".join(title.split())))
    plot = Plotter(filename,figsize=(5,5),dpi=1200,bar_alpha=0.5)
    dists = {}
    if vars == None or len(vars) == 0:
        # This means only the caller wants each seeds to be plotted separated.
        vars = [""]
    for v in vars:
        for s in seeds:
            for type in ["avg","min"]:
                keys = fixed.copy()
                keys.extend([v,s,type])
                filename = pathFinder.Search(keys)
                #print("Reading [{}]...".format(filename))
                dists[v+s+type] = DataFileReader(filename)

    def __calculate_average_dist_of_chains_(raw_data,smooth=True):
        new_data = np.zeros((2, raw_data.shape[1]))
        new_data[0, :] = (raw_data[0, :]) / 100.0;  # frame index to ns
        nPolymers = raw_data.shape[0] - 1
        # calculate the average of all polymer chains, also convert angstrom to nm
        chunk = raw_data[1:nPolymers + 1, :]
        new_data[1, :] = np.average(chunk, axis=0) /  10.0
        if smooth:
            new_data[1, :] = Smooth(new_data[1, :])
        # dists[variable + seed] is now the average distance from the surface all all polymer chains, subject to smoothing
        return new_data
    def __calculate_percent_of_adsorbed_chains_(raw_data,threshhold_nm=0.4,binsize_ns=0.25,smooth=True):
        nPolymers, nFrames = raw_data.shape
        nPolymers -= 1
        # each timestep is 0.01 nm
        ts = 0.01
        import math
        binCount = int(math.floor(nFrames * ts / binsize_ns))
        result_data = np.zeros((2,binCount))
        for iBin in range(binCount):
            binStart = int(iBin * binsize_ns / ts)
            binEnd = min(  int( (iBin+1) *binsize_ns / ts), nFrames )
            #print(binStart,binEnd)
            totalCount = 0
            adsorbedCount = 0
            for iFrame in range(binStart,binEnd):
                for iPoly in range(1, nPolymers + 1):
                    totalCount += 1
                    if raw_data[iPoly,iFrame] < threshhold_nm*10.0:
                        adsorbedCount += 1
            result_data[0,iBin] = (binEnd + binStart) * 0.5 * ts
            result_data[1,iBin] = 100.0 * adsorbedCount / totalCount
        if smooth:
            result_data[1,:] = Smooth(result_data[1,:],0.75)
        return result_data


    # Data processing
    adsorption_dists = {}
    adsorption_rates = {}
    for v in vars:
        first_flag = True
        for s in seeds:
            adsorption_dists[v+s] = __calculate_average_dist_of_chains_(dists[v+s+"avg"])
            adsorption_rates[v+s] = __calculate_percent_of_adsorbed_chains_(dists[v+s+"min"])
            if first_flag:
                first_flag = False
                adsorption_dists[v + "seeds_avg"] = adsorption_dists[v+s].copy()
                adsorption_rates[v + "seeds_avg"] = adsorption_rates[v+s].copy()
            else:
                adsorption_dists[v + "seeds_avg"] += adsorption_dists[v+s].copy()
                adsorption_rates[v + "seeds_avg"] += adsorption_rates[v+s].copy()
        adsorption_dists[v + "seeds_avg"] /= len(seeds)
        adsorption_rates[v + "seeds_avg"] /= len(seeds)

    max_z = - 9999 # this is to determine the highest point the curves can reach. It affects where we place the legend
    # there are two cases,
    if vars==[""]:
        # case 1, compare results from different seeds:
        indexes = ["seeds_avg"]
        indexes.extend(seeds)
    else:
        # case 2, We plot only the seeds_avg data for each variable
        indexes = [v+"seeds_avg" for v in vars]
    # Add data
    for i,index in enumerate(indexes):
        dist_data = adsorption_dists[index]
        ratio_data = adsorption_rates[index]
        s = plotSchemes[i]
        if index!="seeds_avg":
            # special rule, do not show average adsorption percentage in the bar chart to avoid clutter
            plot.AddSeries(PlotSeries(ratio_data[0],ratio_data[1],s.label,c=s.color,ls=s.linestyle,lw=s.linewidth),0)
        plot.AddSeries(PlotSeries(dist_data[0],dist_data[1],s.label,c=s.color,ls=s.linestyle,lw=s.linewidth),1)
        max_z = max(np.max(dist_data[1]), max_z)

    plot.SetAxes(xmin=0,xmax=5.0,ymin1=0,ymax1=105,
                 xlabel="Simulation Time (ns)",
                 ylabel1="Adsorption Percentage (%)",
                 ymin2=0,ymax2=4.0,
                 ylabel2="Polymer/MMT Distance (nm)")
    plot.Plot()

    if max_z > 3.3:
        loc = 'center right'
    else:
        loc = 'upper right'
    plot.SetLegend(whichaxis=1,loc=loc)
    plot.SetTitle(title)

def ParticleZDistributionSingleCase(prefix:str, title:str, pathFinder:PathFinder, seeds:[str], fixed:[str], vars:[str], plotSchemes:[PlotSeries]):
    # the name of the output jpg file. Split and Rejoin to replace ' ' with '_'
    filename = os.path.join(pathFinder.pwd,"analyze","{}{}.jpg".format(prefix,"_".join(title.split())))
    plot = Plotter(filename,figsize=(5,5),dpi=1200,bar_alpha=0.5)
    dists = {}

    # Read data
    for s in seeds:
        keys = fixed.copy()
        keys.append(s)
        filename = pathFinder.Search(keys)
        dists[s] = DataFileReader(filename)
    # Calculate Seeds Average
    first_flag = True
    for s in seeds:
        if first_flag:
            first_flag = False
            dists["seeds_avg"] = dists[s].copy()
        else:
            dists["seeds_avg"] += dists[s].copy()
    dists["seeds_avg"] /= len(seeds)

    # Add data
    ColumnMeanings = {'CATIONS':1, 'ANIONS':2, 'COO':3, 'AMMONIUM':4, 'CONH2':5}
    for i,var in enumerate(vars):
        if var.upper() == "PLACEHOLDER":
            continue
        if var not in ColumnMeanings:
            print("{} not in the reported columns".format(var))
            exit()
        data = dists["seeds_avg"]
        iCol = ColumnMeanings[var]
        s = plotSchemes[i]
        # Smoothing here
        smoothed = Smooth(data[iCol],0.8)
        plot.AddSeries(PlotSeries(data[0],smoothed,s.label,c=s.color,ls=s.linestyle,lw=s.linewidth),1)

    ymax = 4.0
    plot.SetAxes(xmin=0,xmax=3.0,ymin1=0,ymax1=ymax,
                 xlabel="Distance from MMT Surface (nm)",
                 ylabel1="Concentration (mol/L)",
                 ymin2=0,ymax2=ymax,
                 ylabel2=None)
    plot.Plot()

    plot.SetLegend(whichaxis=1,loc="upper right")
    plot.SetTitle(title)

def RDFSingleCase(prefix:str, title:str, pathFinder:PathFinder, seeds:[str], fixed:[str], vars:[str], plotSchemes:[PlotSeries]):
    # the name of the output jpg file. Split and Rejoin to replace ' ' with '_'
    filename = os.path.join(pathFinder.pwd,"analyze","{}{}.jpg".format(prefix,"_".join(title.split())))
    plot = Plotter(filename,figsize=(5,5),dpi=1200,bar_alpha=0.5)
    rdfs = {}

    # Read data
    column_meanings = {}
    for s in seeds:
        keys = fixed.copy()
        keys.append(s)
        filename = pathFinder.Search(keys)
        rdfs[s] = DataFileReader(filename)
        # the 2nd line of file contains the meaning of each column
        if len(column_meanings) != 0:
            continue
        with open(filename,"r",encoding='UTF-8') as datafile:
            datafile.readline()
            line = datafile.readline().strip()
            items = line.split(",")[1:]
            colNo = 1
            for item in items:
                item = item.strip().rstrip()
                if len(item) > 0:
                    column_meanings[item] = colNo
                    colNo += 1

    # Calculate Seeds Average
    first_flag = True
    for s in seeds:
        if first_flag:
            first_flag = False
            rdfs["seeds_avg"] = rdfs[s].copy()
        else:
            rdfs["seeds_avg"] += rdfs[s].copy()
    rdfs["seeds_avg"] /= len(seeds)

    # Add data
    for i,var in enumerate(vars):
        if var.upper() == "PLACEHOLDER":
            continue
        if var not in column_meanings:
            print("{} not in the reported columns".format(var))
            exit()
        data = rdfs["seeds_avg"]
        iCol = column_meanings[var]
        s = plotSchemes[i]
        # Smoothing here
        smoothed = Smooth(data[iCol],0.8)
        plot.AddSeries(PlotSeries(data[0],smoothed,s.label,c=s.color,ls=s.linestyle,lw=s.linewidth),1)

    ymax = 3.0
    plot.SetAxes(xmin=0,xmax=3.0,ymin1=0,ymax1=ymax,
                 xlabel="Distance from MMT Surface (nm)",
                 ylabel1="Radial Distribution Function",
                 ymin2=0,ymax2=ymax,
                 ylabel2=None)
    plot.Plot()

    plot.SetLegend(whichaxis=1,loc="upper right")
    plot.SetTitle(title)

def PlotAllPolymerZCompareSeeds(prefix:str,pf:PathFinder,line):
    schemes = []
    lw = 1.4
    m = None
    ms = 0
    schemes.append(PlotSeries(None, None, None, c='k',  ls='-',lw=lw, m=m, ms=ms))
    schemes.append(PlotSeries(None, None, None, c='r',  ls='-', lw=lw, m=m, ms=ms))
    schemes.append(PlotSeries(None, None, None, c='b',  ls='-', lw=lw, m=m, ms=ms))
    schemes.append(PlotSeries(None, None, None, c='orange',  ls='-', lw=lw, m=m, ms=ms))

    parts = line.split(",")
    title = parts[1].strip().rstrip()
    fixes = parts[2].split()
    vars = parts[3].split()
    seeds = ["seed1","seed2","seed3"]
    if len(vars) == 1 and vars[0] == '*':
        # compare different seeds
        vars = None

    # Now determine the labels, written to schemes.
    # If users have provided the labels, that will be preferred
    if(len(parts)>4):
        labels = parts[4].split()
        for i in range(len(labels)):
            label = labels[i].strip().rstrip()
            schemes[i].label = label.replace('?',' ')

    elif vars != None:
        for i in range(len(vars)):
            schemes[i].label = vars[i]
    else:
        schemes[0].label = "average"
        schemes[1].label = "run 1"
        schemes[2].label = "run 2"
        schemes[3].label = "run 3"

    fixes.append("POLYMERZ")
    PolymerZSingleCase(prefix=prefix,title=title,pathFinder=pf,seeds=seeds,fixed=fixes,vars=vars,plotSchemes=schemes)

def _Extract_Elements(line:str,schemes):
    parts = line.split(",")
    title = parts[1].strip().rstrip()
    seeds = parts[2].strip().rstrip();
    if seeds.upper() == "ALL":
        seeds = ["seed1", "seed2", "seed3"]
    else:
        seeds = seeds.split()

    fixes = parts[3].split()
    vars = parts[4].split()

    # The user must provide labels
    labels = parts[5].split()
    for i in range(len(labels)):
        label = labels[i].strip().rstrip()
        schemes[i].label = label.replace('?',' ')
    return title,seeds,fixes,vars

def PlotParticleZDistribution(prefix:str,pf:PathFinder,line):
    schemes = []
    lw = 1.4
    m = None
    ms = 0
    schemes.append(PlotSeries(None, None, None, c='purple',  ls='-',lw=lw, m=m, ms=ms)) # for cations
    schemes.append(PlotSeries(None, None, None, c='r',  ls='-', lw=lw, m=m, ms=ms)) # for anions
    schemes.append(PlotSeries(None, None, None, c='k',  ls='-', lw=lw, m=m, ms=ms)) # for first group on polymer
    schemes.append(PlotSeries(None, None, None, c='orange',  ls='-', lw=lw, m=m, ms=ms)) # for 2nd group (AADADMAC only)

    title,seeds,fixes,vars = _Extract_Elements(line,schemes)

    fixes.append("DistributionZ")
    ParticleZDistributionSingleCase(prefix=prefix,title=title,pathFinder=pf,seeds=seeds,fixed=fixes,vars=vars,plotSchemes=schemes)

def PlotRDF(prefix:str,pf:PathFinder,line):
    schemes = []
    lw = 1.4
    m = None
    ms = 0
    schemes.append(PlotSeries(None, None, None, c='k',  ls='-',lw=lw, m=m, ms=ms))
    schemes.append(PlotSeries(None, None, None, c='r',  ls='-', lw=lw, m=m, ms=ms))
    schemes.append(PlotSeries(None, None, None, c='b',  ls='-', lw=lw, m=m, ms=ms))
    schemes.append(PlotSeries(None, None, None, c='orange',  ls='-', lw=lw, m=m, ms=ms))

    title,seeds,fixes,vars = _Extract_Elements(line,schemes)

    fixes.append("RDF")
    RDFSingleCase(prefix=prefix,title=title,pathFinder=pf,seeds=seeds,fixed=fixes,vars=vars,plotSchemes=schemes)

def ChainLengthsSingleCase(prefix:str,title:str,pathFinder:PathFinder,seeds:[str],fixed:[str],vars:[str],plotSchemes:[PlotSeries]):
    # the name of the output jpg file. Split and Rejoin to replace ' ' with '_'
    filename = os.path.join(pathFinder.pwd, "analyze", "{}{}.jpg".format(prefix, "_".join(title.split())))
    plot = Plotter(filename, figsize=(5, 5), dpi=1200, bar_alpha=0.5)
    distributions = {} # all data

    # Add data
    for i, var in enumerate(vars):
        # Read data
        for s in seeds:
            keys = fixed.copy()
            keys.append(var)
            keys.append(s)
            filename = pathFinder.Search(keys)
            distributions[var+s] = DataFileReader(filename)

        # Calculate Seeds Average. Assuming all seeds have identical shape of data
        first_flag = True
        for s in seeds:
            if first_flag:
                first_flag = False
                distributions[var+"seeds_avg"] = distributions[var+s].copy()
            else:
                distributions[var+"seeds_avg"] += distributions[var+s].copy()
        distributions[var+"seeds_avg"] /= len(seeds)

        x = distributions[var+"seeds_avg"][0]
        y = distributions[var+"seeds_avg"][1] # need to smooth it? No?  y = Smooth(y,0.8) ?
        s = plotSchemes[i]
        plot.AddSeries(PlotSeries(x,y,label=var, c=s.color, ls=s.linestyle, lw=s.linewidth), 0)

    ymax = 0.08
    plot.SetAxes(xmin=1.5, xmax=3.0, ymin1=0, ymax1=ymax,
                 xlabel="Length of Chain (nm)",
                 ylabel1="Probability",
                 ymin2=0, ymax2=ymax,
                 ylabel2=None)
    plot.Plot()
    plot.SetLegend(whichaxis=0, loc="upper right")
    plot.SetTitle(title)
def PlotChainLengths(prefix:str,pf:PathFinder,line):
    schemes = []
    lw = 1.4
    m = None
    ms = 0
    schemes.append(PlotSeries(None, None, None, c='k',  ls='-',lw=lw, m=m, ms=ms))
    schemes.append(PlotSeries(None, None, None, c='r',  ls='-', lw=lw, m=m, ms=ms))
    schemes.append(PlotSeries(None, None, None, c='b',  ls='-', lw=lw, m=m, ms=ms))
    schemes.append(PlotSeries(None, None, None, c='orange',  ls='-', lw=lw, m=m, ms=ms))

    title,seeds,fixes,vars = _Extract_Elements(line,schemes)
    fixes.append("ChainLengths")
    ChainLengthsSingleCase(prefix=prefix,title=title,pathFinder=pf,seeds=seeds,fixed=fixes,vars=vars,plotSchemes=schemes)


def Main():
    if(len(sys.argv)>1):
        pwd = sys.argv[1]
    else:
        pwd = "./"
    pf = PathFinder(pwd)

    task_file = None
    task_file_name = os.path.join(pf.pwd,"analyze.txt")
    try:
        task_file = open(task_file_name,"r",encoding='UTF-8')
    except:
        print("Can't open task file [{}], don't know what to do...".format(task_file_name))

    counter = 0;
    for line in task_file.readlines():
        line = line.strip()
        if(line.startswith('#') or len(line)==0):
            continue
        counter+=1
        numbering = "{:03d}.".format(counter)
        if line.upper().startswith("POLYMERZ"):
            PlotAllPolymerZCompareSeeds(numbering,pf,line)
        elif line.upper().startswith("ZDISTRIBUTION"):
            PlotParticleZDistribution(numbering,pf,line)
        elif line.upper().startswith("RDF"):
            PlotRDF(numbering,pf,line)
        elif line.upper().startswith("CHAINLENGTHS"):
            PlotChainLengths(numbering,pf,line)
        elif line.upper().startswith("EXIT"):
            break
        else:
            print("Error reading [{}] at the following line:\n{}".format(task_file_name,line))


if __name__ == '__main__':
    Main()
