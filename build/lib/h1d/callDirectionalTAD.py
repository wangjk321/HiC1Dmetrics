from .calculateMetrics import *
from scipy.signal import argrelextrema
from .plotMetrics import *
from .plotTwoSample import *
import sys

class PlotTAD(PlotTri):
    def drawTAD(self,squareSize=300000,useNA=True):
        Tad = TADcallIS(self.path,self.resolution,self.chr,squareSize,useNA=useNA)
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

class stripeTAD(BasePara):
    def callStripe(self,squareSize=300000,useNA=True,seg=100000):
        Tad = TADcallIS(self.path,self.resolution,self.chromosome,squareSize,useNA=useNA)
        intraScore = intraTADscore(self.path,self.resolution,self.chromosome).getIntraS().iloc[:,3]
        nonNAIntraScore = intraScore[~intraScore.isnull()]
        bm = nonNAIntraScore.sample(300,random_state=1)

        pl = []
        pr = []
        import statsmodels.stats.weightstats as sw
        from statsmodels.sandbox.stats.multicomp import multipletests
        for i in range(Tad.shape[0]):
            regionLeft = int((Tad.iloc[i,1])/self.resolution)
            regionRight = int((Tad.iloc[i,2])/self.resolution)
            scorei = intraScore.iloc[regionLeft:regionRight]
            segbin = int(seg/self.resolution) #bins of the corner
            l = scorei.iloc[0:segbin]
            r = scorei.iloc[-segbin:]

            pvalue_l = sw.ztest(bm, value=l.mean(), alternative="smaller")[1]
            pvalue_r = sw.ztest(bm, value=r.mean(), alternative="smaller")[1]
            pl.append(pvalue_l)
            pr.append(pvalue_r)

        qvalue_l = multipletests(pl, method='bonferroni')[1]
        qvalue_r = multipletests(pr, method='bonferroni')[1]
        status = []
        for i in range(Tad.shape[0]):
            if qvalue_l[i] < 0.05 and qvalue_r[i] >0.05:
                status.append("leftStripe")
            elif qvalue_r[i] <0.05 and qvalue_l[i] > 0.05:
                status.append("rightStripe")
            elif qvalue_r[i] <0.05 and qvalue_l[i] < 0.05:
                status.append("loopTAD")
            else:
                status.append("otherTAD")
        Tad["TADtype"] = status
        return(Tad)
        #m = scorei.iloc[seg:-seg]
        #ls = l.min() > m.mean() and l.mean() > interScore.median()
        #rs = r.min() > m.mean() and r.mean() > interScore.median()

class DirectionalTAD(DiffDraw):
    def extractRegion(self):
        ControlTad = TADcallIS(self.control_path,self.resolution,self.chr,squareSize= self.sizeIS)
        drf = DirectionalRelativeFreq(self.path,self.control_path,self.resolution,self.chr, \
                start_distance=self.startDRF,end_distance=self.sizeDRF).getDRF().iloc[:,3]

        MinDRF=[]
        MaxDRF=[]
        for i in range(ControlTad.shape[0]):
            regionLeft = int((ControlTad.iloc[i,1])/self.resolution)
            regionRight = int((ControlTad.iloc[i,2])/self.resolution)
            MinDRF.append(drf.iloc[regionLeft:regionRight].min())
            MaxDRF.append(drf.iloc[regionLeft:regionRight].max())

        leftTAD = ControlTad[np.nan_to_num(MaxDRF)<0]
        rightTAD = ControlTad[np.nan_to_num(MinDRF)>0]
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
        plt.figure(figsize=(7*num, 7+7*(10/6)))

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
