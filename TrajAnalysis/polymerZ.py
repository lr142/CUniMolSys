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
        file = open(filename,"r")
    except:
        print("Can't open {}".format(file))
    # 1st line is comment
    file.readline()
    parts = file.readline().strip().split(',')
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
        mpl.rcParams["font.sans-serif"] = ["Helvetica"]
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
                self.ax[i].bar(s.x+width*(j-NSeriesMid),s.y,width=width,alpha=self.bar_alpha,edgecolor=None,facecolor=s.color,fill=True)
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
                 ylabel2="Polymer/Clay Distance (nm)")
    plot.Plot()

    if max_z > 3.0:
        loc = 'center right'
    else:
        loc = 'upper right'
    plot.SetLegend(whichaxis=1,loc=loc)
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

def Main():
    if(len(sys.argv)>1):
        pwd = sys.argv[1]
    else:
        pwd = "./"
    pf = PathFinder(pwd)

    task_file = None
    task_file_name = os.path.join(pf.pwd,"analyze.txt")
    try:
        task_file = open(task_file_name,"r")
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
        else:
            print("Error reading [{}] at the following line:\n{}".format(task_file_name,line))


if __name__ == '__main__':
    Main()
