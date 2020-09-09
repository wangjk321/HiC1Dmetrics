import os
os.chdir("MainCode")

from calculateMetrics import *
from scipy.signal import argrelextrema
from plotMetrics import *
from plotDiff import *
import sys

def TADcallIS(matrixPath,resolution,chromosome,squareSize=150000):
    ISbedgraph = InsulationScore(matrixPath,resolution,chromosome,square_size=squareSize).getIS()
    ISone = ISbedgraph.InsulationScore

    # local minimal
    localMinPos = argrelextrema(np.array(ISone), np.less)
    localMinIS = ISone.iloc[localMinPos]

    # 0< IS <0.8
    localMinIS = localMinIS[localMinIS!=0]
    localMinIS = localMinIS[localMinIS<0.8]

    #Around TAD boundary
    binNum = int(100000/resolution)
    localMinAround = []
    aroundZero=[]
    diffrightleft=[]

    for i in localMinIS.index:
        # local minimal around 100kb
        localMinAround.append(ISone.loc[i-binNum:i+binNum].min())

        #IS_around !=0
        aroundZero.append(ISone.loc[i-binNum*2:i+binNum*2].min())

        #strong boundary
        minusbin = ISone.loc[i-binNum] - ISone.loc[i]
        plusbin = ISone.loc[i+binNum] - ISone.loc[i]
        diffrightleft.append(minusbin+plusbin)

    bool1 = (np.array(localMinIS)-np.array(localMinAround)) <= 0
    bool2 = np.array(aroundZero)>0
    bool3 = np.array(diffrightleft)>0.05
    localMinIS = localMinIS[bool1 * bool2 * bool3]

    # build a table as output
    TADnumber = len(localMinIS)
    chrlist = [chromosome] * (TADnumber-1)
    TADstart = (np.array(localMinIS.index)[:-1])*resolution
    TADend = (np.array(localMinIS.index)[1:])*resolution
    TADout = pd.DataFrame()
    TADout["chr"] = chrlist
    TADout["TADstart"] = TADstart
    TADout["TADend"] = TADend
    #Maximum TAD size 2MB
    TADout = TADout[(TADout["TADend"]-TADout["TADstart"])<2000000]

    return(TADout)

class PlotTAD(PlotTri):

    def drawTAD(self,squareSize=150000):
        Tad = TADcallIS(self.path,self.resolution,self.chr,squareSize)
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

        super().draw()
        plt.plot(xpos,ypos,"k-",linestyle ="dashed")

    def makePDF(self,PDFname,squareSize=150000):
        self.drawTAD(squareSize)
        plt.savefig(PDFname+".pdf")


class DirectionalTAD(DiffDraw):
    def extractRegion(self):
        ControlTad = TADcallIS(self.control_path,self.resolution,self.chr,squareSize= self.sizeIS)
        drf = DirectionalRelativeFreq(self.path,self.control_path,self.resolution,self.chr, \
                start_distance=self.startDRF,end_distance=self.sizeDRF).getDRF().DirectionalRelativeFreq

        MinDRF=[]
        MaxDRF=[]
        for i in range(ControlTad.shape[0]):
            regionLeft = int((ControlTad.iloc[i,1])/self.resolution)
            regionRight = int((ControlTad.iloc[i,2])/self.resolution)
            MinDRF.append(drf.iloc[regionLeft:regionRight].min())
            MaxDRF.append(drf.iloc[regionLeft:regionRight].max())

        leftTAD = ControlTad[np.array(MaxDRF)<0]
        rightTAD = ControlTad[np.array(MinDRF)>0]
        return(leftTAD,rightTAD,drf)

    def plotAlldirec(self,type):
        if type == "left":
            region,_,drf = self.extractRegion()
        elif type == "right":
            _,region,drf = self.extractRegion()
        else:
            sys.exit("Please specific 'left' or 'right' mode")

        #region = region.iloc[0:2,:]
        num = region.shape[0]
        plt.figure(figsize=(7*num, 7+7*(10/6)),dpi=300)

        for i in range(num):
            TADleft = region.iloc[i,1] # directional TAD region
            TADright = region.iloc[i,2]
            start = TADleft - 1000000 # plot region
            end = TADright + 1000000
            sbin = int(start/self.resolution)
            ebin = int(end/self.resolution)
            drfRegion = drf[sbin:ebin+1]

            plotleft = ((TADleft-start)/self.resolution+1)*1.4142 #plot region of TAD
            plotright = ((TADright-start)/self.resolution+1)*1.4142

            #plot1: diffiential contact matrix
            ax1 = plt.subplot2grid((17,11*num),(0,i*11),rowspan=5,colspan=10)
            difftri = DiffDraw(self.path,self.control_path,self.resolution,startSite=start, \
                        endSite=end, clmin=-2,clmax=2,chr=self.chr,title="#"+str(i+1)+" Directional TAD")
            difftri.draw_tri()

            ## show the TAD
            TAD_xpos = [plotleft,(plotleft+plotright)/2,plotright]
            TAD_ypos = [ax1.get_ylim()[0],ax1.get_ylim()[0]-(plotright-plotleft)/2,ax1.get_ylim()[0]]
            ax1.plot(TAD_xpos,TAD_ypos,"k-")

            #plot2: DRF score
            ax2 = plt.subplot2grid((17,11*num),(5,i*11),rowspan=1,colspan=11)
            plt.title("DirectionalRelativeFreq",fontsize=20)
            plt.plot(drfRegion,c='black')
            plt.xlim(sbin,ebin)
            plt.plot([sbin,ebin],[drf.mean(),drf.mean()],"k--",linewidth=0.4)
            plt.fill_between(np.arange(sbin,ebin+1,1),drfRegion,\
                            where = drfRegion >0,facecolor='#e9a3c9', alpha=0.99)

            plt.fill_between(np.arange(sbin,ebin+1,1),drfRegion,\
                            where = drfRegion <0,facecolor='#a1d76a', alpha=0.99)

            ticks_pos = np.arange(sbin,ebin+1,(ebin-sbin)/5)
            plt.xticks(ticks_pos,difftri.mark)

            #plot3: Control matrix
            ax3 = plt.subplot2grid((17,11*num),(7,i*11),rowspan=5,colspan=10)
            #PlotTAD(self.control_path,self.resolution,startSite=start,endSite=end, \
            #    clmin=0,clmax=50,chr=self.chr,title="#"+str(i+1)+" Contact matrix of Control").drawTAD(self.sizeIS)
            PlotTri(self.control_path,self.resolution,startSite=start,endSite=end,
                clmin=0,clmax=50,chr=self.chr,title="#"+str(i+1)+" Contact matrix of Control").draw()
            ax3.plot(TAD_xpos,TAD_ypos,"k-")

            #plot4: Treated matrix
            ax4 = plt.subplot2grid((17,11*num),(12,i*11),rowspan=5,colspan=10)
            PlotTri(self.path,self.resolution,startSite=start,endSite=end,
                clmin=0,clmax=50,chr=self.chr,title="#"+str(i+1)+" Contact matrix of treated").draw()
            ax4.plot(TAD_xpos,TAD_ypos,"k-")

    def makePDF(self,type,outname):
        self.plotAlldirec(type)
        plt.savefig(outname,bbox_inches='tight')

"""
cla = DirectionalTAD("../test_data/Rad21KD_1/observed.KR.chr21.matrix.gz","../test_data/Control_1/observed.KR.chr21.matrix.gz",25000,chr="chr21")
cla.makePDF("left","testleft.pdf")
PlotSquare("../test_data/Control_1/observed.KR.chr21.matrix.gz",25000,startSite=1100*25000,endSite=1300*25000,chr="chr21").draw()
plt.savefig("SquareHiC.pdf")

PlotBedGraph("../test_data/Control_1/observed.KR.chr21.matrix.gz",25000,startSite=1100*25000,endSite=1300*25000,title="One sample",chr="chr21").draw("IS")
plt.savefig("ISsingle.pdf")

diffsq = DiffDraw("../test_data/Rad21KD_1/observed.KR.chr21.matrix.gz","../test_data/Control_1/observed.KR.chr21.matrix.gz",25000,\
                    startSite=1100*25000,endSite=1300*25000, clmin=-3,clmax=3,title="Comparison of two samples",chr="chr21")
diffsq.draw_DRF()
plt.savefig("DRFtwo.pdf")"""

#leftTAD,rightTAD,_ = cla.extractRegion()
#leftTAD
#leftTAD.to_csv("leftTAD" + ".bedGraph", sep="\t", header=False, index=False)
#rightTAD.to_csv("rightTAD" + ".bedGraph", sep="\t", header=False, index=False)

IS = InsulationScore("../test_data/Rad21KD_1/observed.KR.chr21.matrix.gz",25000,"chr21")
IS.getIS()
direcTAD = DirectionalTAD("../test_data/Rad21KD_1/observed.KR.chr21.matrix.gz","../test_data/Control_1/observed.KR.chr21.matrix.gz",25000,chr="chr21",startDRF=500000,sizeDRF=1000000,sizeIS=150000)
leftTAD,rightTAD,_ = direcTAD.extractRegion()
rightTAD
plotDirecTAD = DirectionalTAD("../test_data/Rad21KD_1/observed.KR.chr21.matrix.gz","../test_data/Control_1/observed.KR.chr21.matrix.gz",25000,chr="chr21",clmin=-2,clmax=2,title="Directional TAD on chr21",startDRF=500000,sizeDRF=1000000,sizeIS=150000)
plotDirecTAD.plotAlldirec("right")
