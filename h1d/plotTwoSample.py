import matplotlib.pyplot as plt
from .loadfile import *
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from scipy import ndimage
from .calculateMetrics import *
from .calculateTwoSample import *
import os
#from calcuDiffCI import *

cmap= LinearSegmentedColormap.from_list("custom2",['#1310cc', '#FFFFFF', '#d10a3f'])

class DiffDraw(object):
    def __init__(self,path,control_path,resolution,chr="",startSite=0,endSite=0,clmin=-2,clmax=2, \
                title="", startDRF=500000,sizeDRF=2000000,sizeIS=300000,sizeDCI=300000,datatype="matrix",gt=None):

        if datatype == "matrix":
            self.path = path
            self.control_path = control_path
        elif datatype == "rawhic":
            if not gt:
                raise ValueError("Genometable is required for the calculation of IF")
            self.path = hic2matrix(path,resolution,chr,gt)
            self.control_path = hic2matrix(control_path,resolution,chr,gt)
        else: print("Error, please use the matrix or rawhic datatype")
        print(self.path,"diff")
        treat = loadWithNorm(self.path,log = True).values
        control = loadWithNorm(self.control_path,log = True).values
        smooth = 2
        self.treat = ndimage.median_filter(treat,smooth)
        self.control = ndimage.median_filter(control,smooth)
        self.matrix = self.treat-self.control

        self.resolution = resolution
        self.startSite = startSite
        self.endSite = endSite
        self.clim = (clmin,clmax)
        self.title = title
        sbin = int(startSite/resolution)
        ebin = int(endSite/resolution)
        self.chr = chr
        self.sizeDCI = sizeDCI
        self.sizeDRF = sizeDRF
        self.startDRF = startDRF
        self.sizeIS = sizeIS

        if endSite == 0:
            self.plotall = True
            self.matrixRegion = self.matrix.copy()
            self.matrixRegion = np.nan_to_num(self.matrixRegion)
            ebin = self.matrix.shape[0]
        elif (ebin - sbin) < 5:
            print("The region you choosed is too small")
            exit(1)
        else:
            self.matrixRegion = self.matrix[sbin:ebin+1,sbin:ebin+1]

        self.sbin = sbin
        self.ebin = ebin

        position = self.resolution * np.arange(self.sbin,self.ebin+1,int((self.ebin-self.sbin)/5))
        self.mark = [str(x/1000000)+"M" for x in position]

    def draw_square(self):
        fig = plt.figure(figsize=(5,5))
        plt.imshow(self.matrixRegion,clim= self.clim,cmap=cmap,interpolation="nearest", aspect=1)
        plt.title(self.title)

        ticks_pos = list(range(0,self.ebin-self.sbin+1,int((self.ebin-self.sbin)/5)))

        plt.yticks(ticks_pos,self.mark)
        plt.xticks([])
        position=fig.add_axes([0.9, 0.73, 0.015, 0.15])
        plt.colorbar(fraction = 0.02,aspect=10,pad=-0.12,cax=position)

    def draw_tri(self):
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

    def drawTAD(self,squareSize=300000,useNA=True):
        Tad = TADcallIS(self.control_path,self.resolution,self.chr,squareSize,useNA=useNA)
        selectTADbool = np.logical_and(Tad["TADstart"] >= self.startSite,Tad["TADend"] <= self.endSite)
        selectTAD=Tad[selectTADbool]

        originalWidth = self.matrixRegion.shape[0]
        figwidth = ndimage.rotate(self.matrixRegion,45).shape[0]

        xpos=[]
        ypos=[]
        for i in range(selectTAD.shape[0]):
            left = ((selectTAD["TADstart"].iloc[i]-self.startSite)/self.resolution+1) * (figwidth/originalWidth)
            right = ((selectTAD["TADend"].iloc[i]-self.startSite)/self.resolution+1) * (figwidth/originalWidth)
            xpos.append(left)
            xpos.append((left+right)/2)
            xpos.append(right)

            ypos.append(figwidth/2)
            ypos.append(figwidth/2-(right-left)/2)
            ypos.append(figwidth/2)

        self.draw_tri()
        plt.plot(xpos,ypos,"k-",linestyle ="dashed")

    def drawMetric(self,type,customfile="",name="custom"):
        if type == "deltaDLR":
            score = deltaDLR(self.path,self.control_path,self.resolution,self.chr, \
                    sizeDLR=3000000).getDeltaDLR().deltaDLR
            title = "deltaDLR"
        elif type == "ISC":
            score = TADScoreChange(self.path,self.control_path,self.resolution,
                                    self.chr).getChange("IS",300000).InsulationScoreChange
            title = "InsulationScoreChange"
        elif type == "CIC":
            score = TADScoreChange(self.path,self.control_path,self.resolution,
                                    self.chr).getChange("CI",300000).ContrastIndexChange
            title = "ContrastIndexChange"
        elif type == "DRF":
            score = DirectionalRelativeFreq(self.path,self.control_path,self.resolution,self.chr, \
                    start_distance=self.startDRF,end_distance=self.sizeDRF).getDRF().DirectionalRelativeFrequency
            title = "DirectionalRelativeFreq"
        elif type == "intraSC":
            score = intraScoreChange(self.path,self.control_path,self.resolution,self.chr).getIntraSC().IntraSC
            title = "intraScoreChange"
        elif type == "interSC":
            score = interScoreChange(self.path,self.control_path,self.resolution,self.chr).getInterSC().InterSC
            title = "interScoreChange"
        elif type == "custom":
            score = pd.read_csv(customfile,sep="\t",header=None).iloc[:,3]
            title = name

        scoreRegion = score[self.sbin:self.ebin+1]
        plt.figure(figsize=(10,10))
        #pad + fraction = -0.1, So the bedGraph figure width should be (1+0.1)*width of HiC
        plt.subplot2grid((6,11),(0,0),rowspan=5,colspan=10)
        self.draw_tri()
        plt.subplot2grid((6,11),(5,0),rowspan=1,colspan=11)
        plt.title(title,fontsize=20)
        plt.plot(scoreRegion,c='black')
        plt.xlim(self.sbin,self.ebin)
        plt.plot([self.sbin,self.ebin],[score.mean(),score.mean()],"k--",linewidth=0.4)
        if self.plotall: self.ebin = self.ebin-1
        plt.fill_between(np.arange(self.sbin,self.ebin+1,1),scoreRegion,\
                        where = scoreRegion >0,facecolor='#e9a3c9', alpha=0.5)

        plt.fill_between(np.arange(self.sbin,self.ebin+1,1),scoreRegion,\
                        where = scoreRegion <0,facecolor='#a1d76a', alpha=0.5)

        ticks_pos = np.arange(self.sbin,self.ebin+1,(self.ebin-self.sbin)/5)
        plt.xticks(ticks_pos,self.mark)
        os.system("rm -rf MatrixTemp0*")


    def draw_DRF(self):  #已被上面的替代
        score = DirectionalRelativeFreq(self.path,self.control_path,self.resolution,self.chr, \
                start_distance=self.startDRF,end_distance=self.sizeDRF).getDRF().DirectionalRelativeFreq
        title = "DirectionalRelativeFreq"
        scoreRegion = score[self.sbin:self.ebin+1]

        plt.figure(figsize=(10,10))
        #pad + fraction = -0.1, So the bedGraph figure width should be (1+0.1)*width of HiC
        plt.subplot2grid((6,11),(0,0),rowspan=5,colspan=10)
        self.draw_tri()
        plt.subplot2grid((6,11),(5,0),rowspan=1,colspan=11)
        plt.title(title,fontsize=20)
        plt.plot(scoreRegion,c='black')
        plt.xlim(self.sbin,self.ebin)
        plt.plot([self.sbin,self.ebin],[score.mean(),score.mean()],"k--",linewidth=0.4)
        if self.plotall: self.ebin = self.ebin-1
        plt.fill_between(np.arange(self.sbin,self.ebin+1,1),scoreRegion,\
                        where = scoreRegion >0,facecolor='#e9a3c9', alpha=0.99)

        plt.fill_between(np.arange(self.sbin,self.ebin+1,1),scoreRegion,\
                        where = scoreRegion <0,facecolor='#a1d76a', alpha=0.99)

        ticks_pos = np.arange(self.sbin,self.ebin+1,(self.ebin-self.sbin)/5)
        plt.xticks(ticks_pos,self.mark)

    def draw_diffCI(self):
        score = DiffCI(self.path,self.control_path,self.resolution,self.chr,diffCI_size=self.sizeDCI).getDiffCI().diffCI
        title = "Diff CI"
        scoreRegion = score[self.sbin:self.ebin+1]

        plt.figure(figsize=(10,10))
        #pad + fraction = -0.1, So the bedGraph figure width should be (1+0.1)*width of HiC
        plt.subplot2grid((6,11),(0,0),rowspan=5,colspan=10)
        self.draw_tri()
        ciplot= plt.subplot2grid((6,11),(5,0),rowspan=1,colspan=11)
        plt.title(title,fontsize=20)
        plt.plot(scoreRegion,c='black')
        plt.xlim(self.sbin,self.ebin)
        plt.plot([self.sbin,self.ebin],[0,0],"k--",linewidth=0.6)
        if np.mean(scoreRegion)<0:
            plt.fill_between(np.arange(self.sbin,self.ebin+1,1),scoreRegion, 0,\
                            where = scoreRegion <0,facecolor='blue', alpha=0.5)
        else:
            plt.fill_between(np.arange(self.sbin,self.ebin+1,1),scoreRegion, 0,\
                            where = scoreRegion >0,facecolor='red', alpha=0.5)


        ticks_pos = np.arange(self.sbin,self.ebin+1,(self.ebin-self.sbin)/5)
        plt.xticks(ticks_pos,self.mark)
