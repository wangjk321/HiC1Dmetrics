from calculateMetrics import *
from calculateTwoSample import *
from hmmlearn import hmm
import matplotlib.pyplot as plt
import seaborn as sns
from plotTwoSample import *
from plotMetrics import *
from scipy import stats
from callDirectionalTAD import *

class multiScore:
    def __init__(self,path,res,chr,control_path=""):
        self.path = path
        self.res = res
        self.chr = chr
        self.control_path = control_path

    def obtainOneScore(self,mode,parameter,smoothPC=True,logPC=False,
                        custom_name="InteractionFrequency",normIF=False):
        if mode == "IS":
            score = InsulationScore(self.path,self.res,self.chr,square_size=parameter).getIS()
        elif mode == "CI":
            score = ContrastIndex(self.path,self.res,self.chr,CI_size=parameter).getCI()
        elif mode == "DI":
            score = DirectionalityIndex(self.path,self.res,self.chr,DI_distance=parameter).getDI()
        elif mode == "TADss":
            score = SeparationScore(self.path,self.res,self.chr,TADss_size=parameter).getTADss()
        elif mode == "DLR":
            score = DistalToLocalRatio(self.path,self.res,self.chr,sizeDLR=parameter).getDLR()
        elif mode == "intraS":
            score = intraTADscore(self.path,self.res,self.chr).getIntraS(IS_size = parameter)
        elif mode == "interS":
            score = interTADscore(self.path,self.res,self.chr).getInterS(IS_size = parameter)
        elif mode == "PC1":
            score = CompartmentPC1(self.path,self.res,self.chr).getPC1(signCorr = parameter,smooth = smoothPC, logOE=logPC)
        elif mode == "custom":
            all = pd.read_csv(parameter,sep="\t",header=None)
            score = all[all[0] == self.chr]
            if normIF:
                beforlog = score[3].copy()
                afterlog = np.log1p(beforlog)
                score[3] = afterlog / np.mean(afterlog[afterlog>0])
            score.index = range(score.shape[0])
            score.columns = ["chr","start","end",custom_name]
        elif mode == "stripe":
            score = stripeTAD(self.path,self.res,self.chr).callStripe(seg=parameter)
        else: print("Error: Please use the right mode")
        return(score)

    def allOneScore(self,typelist=["IS","CI","DI","TADss","DLR","intraS","interS","PC1","custom"],
                    parameterlist=[300000,300000,1000000,300000,3000000,300000,300000,"NotSpecified","customPath"],
                    smoothPC=True,logPC=False):
        for i,type in enumerate(typelist):
            if i == 0:
                multiType = self.obtainOneScore(mode=typelist[i],parameter=parameterlist[i],smoothPC=smoothPC,logPC=logPC)
            else:
                next = self.obtainOneScore(mode=typelist[i],parameter=parameterlist[i],smoothPC=smoothPC,logPC=logPC).iloc[:,3]
                multiType = pd.concat([multiType,next],axis=1)

        return(multiType)

    def plotOneScore(self,start,end,res,clmax=100,plotTAD=False,
                    typelist=["IS","CI","DI","TADss","DLR","intraS","interS","PC1","custom"],
                    parameterlist=[300000,300000,1000000,300000,3000000,300000,300000,"NotSpecified","customPath"],
                    smoothPC=True,logPC=False):
        import matplotlib.colors as mcolors
        from callDirectionalTAD import PlotTAD

        cols = list(mcolors.TABLEAU_COLORS.keys())
        scoreMT = self.allOneScore(typelist,parameterlist,smoothPC,logPC)
        nScore = len(typelist)

        plt.figure(figsize=(10,9+nScore*1))
        plt.subplot2grid((5+nScore,11),(0,0),rowspan=5,colspan=10)
        hp = PlotTAD(self.path,res,start,end,clmax=clmax)
        if plotTAD == True:
            hp.drawTAD()
        elif plotTAD == False:
            hp.draw()
        for i in range(nScore):
            scoreRegion = scoreMT.iloc[start//res:end//res,3+i]
            plt.subplot2grid((5+nScore,11),(5+i,0),rowspan=1,colspan=11)
            plt.plot(scoreRegion,c=cols[i],label=scoreRegion.name)
            plt.legend(loc="upper left")
            plt.xlim(start//res,end//res)
            ticks_pos = np.arange(hp.sbin,hp.ebin+1,(hp.ebin-hp.sbin)/5)
            if i < nScore-1:
                plt.xticks(ticks_pos,[])
            else:
                plt.xticks(ticks_pos,hp.mark)

    def obtainTwoScore(self,mode,parameter,smoothPC=True,logPC=False,normIF=True):
        if mode == "ISC":
            score = TADScoreChange(self.path,self.control_path,self.res,self.chr).getChange("IS",parameter)
        elif mode == "CIC":
            score = TADScoreChange(self.path,self.control_path,self.res,self.chr).getChange("CI",parameter)
        elif mode == "DIC":
            score = TADScoreChange(self.path,self.control_path,self.res,self.chr).getChange("DI",parameter)
        elif mode == "TADssC":
            score = TADScoreChange(self.path,self.control_path,self.res,self.chr).getChange("TADss",parameter)
        elif mode == "deltaDLR":
            score = deltaDLR(self.path,self.control_path,self.res,self.chr,sizeDLR=parameter).getDeltaDLR()
        elif mode == "intraSC":
            score = intraScoreChange(self.path,self.control_path,self.res,self.chr,IS_size=parameter).getIntraSC()
        elif mode == "interSC":
            score = interScoreChange(self.path,self.control_path,self.res,self.chr,IS_size=parameter).getInterSC()
        elif mode == "DRF":
            score = DirectionalRelativeFreq(self.path,self.control_path,self.res,self.chr,
                                            start_distance=parameter[0], end_distance=parameter[1]).getDRF()
        elif mode == "CorrD":
            score = CorrelationDifference(self.path,self.control_path,self.res,self.chr,method=parameter).getCorrD()
        elif mode == "PC1C":
            score = PC1change(self.path,self.control_path,self.res,self.chr,corr_file=parameter,smoothPC = smoothPC, logPC=logPC).getPC1change()
        elif mode == "IFC":
            scoreTreat = pd.read_csv(parameter[0],sep="\t",header=None)
            scoreControl= pd.read_csv(parameter[1],sep="\t",header=None)
            treat = scoreTreat[scoreTreat[0]==self.chr]
            control = scoreControl[scoreControl[0]==self.chr]
            treatlog = np.log1p(treat.iloc[:,3])
            controllog = np.log1p(control.iloc[:,3])
            t = treatlog / np.mean(treatlog[treatlog>0])
            c = controllog / np.mean(controllog[controllog>0])
            
            score = pd.DataFrame({"chr":treat[0],"start":treat[1],"end":treat[2],"IFchange":t-c})
            score.index = range(score.shape[0])
        else: print("Error: Please specify the correct mode")

        return(score)

    def allTwoScore(self,typelist=["ISC","CIC","DIC","TADssC","deltaDLR","intraSC","interSC","DRF","CorrD","PC1C"],
                    parameterlist=[300000,300000,1000000,300000,3000000,300000,300000,[200000,5000000],"pearson","NotSpecified"],smoothPC=True,logPC=False):
        for i,type in enumerate(typelist):
            if i == 0:
                multiType = self.obtainTwoScore(mode=typelist[i],parameter=parameterlist[i],smoothPC=smoothPC,logPC=logPC)
            else:
                next = self.obtainTwoScore(mode=typelist[i],parameter=parameterlist[i],smoothPC=smoothPC,logPC=logPC).iloc[:,3]
                multiType = pd.concat([multiType,next],axis=1)

        return(multiType)

    def plotTwoScore(self,start,end,res,clmax=2,plotTAD=False,
                    typelist=["ISC","CIC","DIC","TADssC","deltaDLR","intraSC","interSC","DRF","CorrD","PC1C"],
                    parameterlist=[300000,300000,1000000,300000,3000000,300000,300000,[200000,5000000],"pearson","NotSpecified"],
                    smoothPC=True,logPC=False):
        import matplotlib.colors as mcolors
        from callDirectionalTAD import PlotTAD

        cols = list(mcolors.TABLEAU_COLORS.keys())
        scoreMT = self.allTwoScore(typelist,parameterlist,smoothPC,logPC)
        nScore = len(typelist)

        plt.figure(figsize=(10,9+nScore*1))
        plt.subplot2grid((5+nScore,11),(0,0),rowspan=5,colspan=10)
        hp = DiffDraw(self.path,self.control_path,res,startSite=start,endSite=end,clmax=clmax)
        if plotTAD == True:
            hp.drawTAD()
        elif plotTAD == False:
            hp.draw_tri()
        for i in range(nScore):
            scoreRegion = scoreMT.iloc[start//res:end//res+1,3+i]
            plt.subplot2grid((5+nScore,11),(5+i,0),rowspan=1,colspan=11)
            plt.plot(scoreRegion,c=cols[i],label=scoreRegion.name)
            plt.legend(loc="upper left")
            plt.xlim(start//res,end//res)

            if scoreRegion.name != "CorrD":
                plt.plot([hp.sbin,hp.ebin],[0,0],"k--",linewidth=0.6)
                plt.fill_between(np.arange(hp.sbin,hp.ebin+1,1),scoreRegion, 0,\
                                where = scoreRegion <=0,facecolor='grey', alpha=0.5)
                plt.fill_between(np.arange(hp.sbin,hp.ebin+1,1),scoreRegion, 0,\
                                where = scoreRegion >=0,facecolor=cols[i], alpha=0.4)

            ticks_pos = np.arange(hp.sbin,hp.ebin+1,(hp.ebin-hp.sbin)/5)
            if i < nScore-1:
                plt.xticks(ticks_pos,[])
            else:
                plt.xticks(ticks_pos,hp.mark)

class metricHMM:
    def __init__(self,df,ncluster,nRun=10,covMethod= "spherical",random_state=None,HMMtype="Gaussian"):
        if len(df.shape) == 1: df= pd.DataFrame(df)
        self.rawdf = df
        self.label = df.columns
        self.df = pd.DataFrame(np.nan_to_num(df))
        self.nRun = nRun
        self.ncluster = ncluster
        self.state = ["state" + str(i) for i in range(ncluster)]
        self.index = df.index
        self.covMethod = covMethod
        self.limx = [min(df.index),max(df.index)]

        scorelist = []
        modellist = []
        if HMMtype == "Gaussian":
            hmmTrain = hmm.GaussianHMM
        elif HMMtype == "GMM":
            hmmTrain = hmm.GMMHMM

        for i in range(self.nRun):
            model = hmmTrain(n_components=self.ncluster, n_iter=10000,random_state=random_state,
                                    covariance_type=self.covMethod).fit(self.df)
            scorelist.append(model.score(self.df))
            modellist.append(model)
        self.bestmodel = modellist[np.argmax(scorelist)]

    def oneSampleMultiMetric(self,outtype="predict"):
        predictMT = self.bestmodel.predict(self.df)
        emissionMT = pd.DataFrame(self.bestmodel.means_,columns=self.label,index=self.state)
        transitionMT = pd.DataFrame(self.bestmodel.transmat_,columns=self.state,index=self.state)

        if outtype == "predict":
            return(predictMT)
        elif outtype == "emission":
            return(emissionMT)
        elif outtype == "transition":
            return(transitionMT)

    def plotHMM(self,outtype="predict",norm="local"):
        mt = self.oneSampleMultiMetric(outtype)
        if outtype == "predict":
            plt.scatter(self.index,mt,c=mt,marker="8")
            plt.ylim(-1,self.ncluster)
            plt.yticks(range(self.ncluster),self.state)
        elif outtype == "emission":
            if norm == "local":
                zMT = mt.apply(stats.zscore)
            elif norm == "overall":
                u = self.rawdf.mean(axis=0)
                uMT = pd.DataFrame([u]*self.ncluster,index=mt.index)
                sd = self.rawdf.std(axis=0)
                sdMT = pd.DataFrame([sd]*self.ncluster,index=mt.index)
                zMT = (mt - uMT)/sdMT
            sns.heatmap(zMT.T,cmap="coolwarm",vmax=2,vmin=-2)
        elif outtype == "transition":
            sns.heatmap(mt,cmap="coolwarm")

    def plotHiC(self,mode,path,resolution,startSite,endSite,control_path="",clmax=50):
        mt = self.oneSampleMultiMetric("predict")
        plt.figure(figsize=(10,10))
        plt.subplot2grid((6,11),(0,0),rowspan=5,colspan=10)
        if mode == "single":
            hicplot = PlotTri(path,resolution,startSite,endSite,clmax=clmax)
            hicplot.draw()
        elif mode == "differ":
            hicplot = DiffDraw(path,control_path,resolution,startSite,endSite)
            hicplot.draw_tri()
        plt.subplot2grid((6,11),(5,0),rowspan=1,colspan=11)
        plt.scatter(self.index,mt,c=mt,marker="8")
        plt.yticks(range(self.ncluster),self.state)
        plt.xlim(self.limx[0],self.limx[1])

    def oneSampleOneMetric(self):
        pass
    def multiSampleOneMetric(self):
        pass
    def multiSampleMultiMetric(self):
        pass
