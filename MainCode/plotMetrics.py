import matplotlib.pyplot as plt
from loadfile import *
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from scipy import ndimage
from calculateMetrics import *

class PlotCommon(object):
    def __init__(self,path,resolution,startSite=0,endSite=0,clmin=0,clmax=50, \
                title="", chr="",other_parameter=0):
        self.path = path
        self.matrix = loadDenseMatrix(path).values
        self.resolution = resolution
        self.startSite = startSite
        self.endSite = endSite
        self.clim = (clmin,clmax)
        self.title = title
        sbin = int(startSite/resolution)
        ebin = int(endSite/resolution)
        self.chr = chr
        self.other_parameter = other_parameter

        if endSite == 0:
            self.matrixRegion = np.nan_to_num(self.matrix.copy())
            ebin = self.matrix.shape[0]
        elif (ebin - sbin) < 5:
            print("The region you choosed is too small")
            exit(1)
        else:
            self.matrixRegion = np.nan_to_num(self.matrix[sbin:ebin+1,sbin:ebin+1])

        self.sbin = sbin
        self.ebin = ebin

        position = self.resolution * np.arange(self.sbin,self.ebin+1,int((self.ebin-self.sbin)/5))
        self.mark = [str(x/1000000)+"M" for x in position]

cmap= LinearSegmentedColormap.from_list("custom1",['#FFFFFF', '#d10a3f'])

class PlotSquare(PlotCommon):
    def draw(self):
        fig = plt.figure(figsize=(5,5))
        plt.imshow(self.matrixRegion,clim= self.clim,cmap=cmap,interpolation="nearest", aspect=1)
        plt.title(self.title,fontsize=20)

        ticks_pos = list(range(0,self.ebin-self.sbin+1,int((self.ebin-self.sbin)/5)))

        plt.yticks(ticks_pos,self.mark)
        plt.xticks([])
        position=fig.add_axes([0.9, 0.73, 0.015, 0.15])
        plt.colorbar(fraction = 0.02,aspect=10,pad=-0.12,cax=position)

class PlotTri(PlotCommon):
    def draw(self):
        tri_matrix = ndimage.rotate(self.matrixRegion,45)
        plt.plot(dpi=300)
        plt.imshow(tri_matrix,clim= self.clim,cmap=cmap,interpolation="nearest", aspect=1)
        plt.title(self.title,fontsize=20)
        tri_shape = tri_matrix.shape[0]
        plt.ylim(tri_shape/2,0-tri_shape/10)
        plt.xlim(0,tri_shape)

        plt.colorbar(fraction = 0.02,aspect=10,pad=-0.12)
        plt.plot([0,tri_shape/2,tri_shape],[tri_shape/2,0,tri_shape/2],"k-",linewidth=0.2)

        ticks_pos = np.arange(0,tri_shape+1,tri_shape/5)

        plt.xticks(ticks_pos,self.mark)
        plt.yticks([])

class PlotBedGraph(PlotTri):

    def draw(self,type,UniqueParameter=None,smoothPC=True,logPC=False,customfile=""):
        if type == 'IS':
            score = InsulationScore(self.path,self.resolution,self.chr).getIS().InsulationScore
            title = "InsulationScore"
        elif type == 'DI':
            score = DirectionalityIndex(self.path,self.resolution,self.chr).getDI().DirectionalityIndex
            title = "DirectionalityIndex"
        elif type == "CI":
            score = ContrastIndex(self.path,self.resolution,self.chr).getCI().ContrastIndex
            title = "ContrastIndex"
        elif type == "TADss":
            score = SeparationScore(self.path,self.resolution,self.chr).getTADss().SeparationScore
            title = "SeparationScore"
        elif type == "DLR":
            score = DistalToLocalRatio(self.path,self.resolution,self.chr).getDLR().DistalToLocalRatio
            title = "DistalToLocalRatio"
        elif type == "PC1":
            score = CompartmentPC1(self.path,self.resolution,self.chr).getPC1(smooth = smoothPC, logOE=logPC).CompartmentPC1
            title = "CompartmentPC1"
        elif type == "intraScore":
            score = intraTADscore(self.path,self.resolution,self.chr).getIntraS().intraTADscore
            title = "intraTADscore"
        elif type == "interScore":
            score = interTADscore(self.path,self.resolution,self.chr).getInterS().interTADscore
            title = "inerTADscore"
        elif type == "custom":
            score = pd.read_csv(customfile,sep="\t",header=None).iloc[:,3]
            title = "customScore"

        scoreRegion = score[self.sbin:self.ebin+1]

        plt.figure(figsize=(10,10))
        #pad + fraction = -0.1, So the bedGraph figure width should be (1+0.1)*width of HiC
        plt.subplot2grid((6,11),(0,0),rowspan=5,colspan=10)
        super().draw()
        plt.subplot2grid((6,11),(5,0),rowspan=1,colspan=11)
        plt.title(title,fontsize=20)
        plt.plot(scoreRegion,c="dodgerblue")
        plt.xlim(self.sbin,self.ebin)
        plt.plot([self.sbin,self.ebin],[score.median(),score.median()],"k--",linewidth=0.4)

        ticks_pos = np.arange(self.sbin,self.ebin+1,(self.ebin-self.sbin)/5)
        plt.xticks(ticks_pos,self.mark)

    def makePDF(self,type,PDFname):
        self.draw(type)
        plt.savefig(PDFname+".pdf")
