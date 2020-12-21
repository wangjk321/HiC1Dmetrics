from calculateMetrics import *
from calculateTwoSample import *
from hmmlearn import hmm
import matplotlib.pyplot as plt
import seaborn as sns
from plotTwoSample import *
from plotMetrics import *

class multiScore:
    def __init__(self,path,res,chr,control_path=""):
        self.path = path
        self.res = res
        self.chr = chr
        self.control_path = control_path

    def obtainOneScore(self,mode,parameter,smoothPC=True,logPC=False):
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
        return(score)

    def allOneScore(self,typelist=["IS","CI","DI","TADss","DLR","intraS","interS","PC1"],
                    parameterlist=[300000,300000,1000000,300000,3000000,300000,300000,"NotSpecified"],smoothPC=True,logPC=False):
        for i,type in enumerate(typelist):
            if i == 0:
                multiType = self.obtainOneScore(mode=typelist[i],parameter=parameterlist[i],smoothPC=smoothPC,logPC=logPC)
            else:
                next = self.obtainOneScore(mode=typelist[i],parameter=parameterlist[i],smoothPC=smoothPC,logPC=logPC).iloc[:,3]
                multiType = pd.concat([multiType,next],axis=1)

        return(multiType)

    def obtainTwoScore(self,mode,parameter,smoothPC=True,logPC=False):
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

class metricHMM:
    def __init__(self,df,ncluster,nRun=10):
        if len(df.shape) == 1: df= pd.DataFrame(df)
        self.label = df.columns
        self.df = pd.DataFrame(np.nan_to_num(df))
        self.nRun = nRun
        self.ncluster = ncluster
        self.state = ["state" + str(i) for i in range(ncluster)]
        self.index = df.index

    def oneSampleMultiMetric(self,outtype="predict"):
        scorelist = []
        modellist = []
        for i in range(self.nRun):
            model = hmm.GaussianHMM(n_components=self.ncluster, n_iter=10000, covariance_type="full").fit(self.df)
            scorelist.append(model.score(self.df))
            modellist.append(model)
        bestmodel = modellist[np.argmax(scorelist)]
        predictMT = bestmodel.predict(self.df)
        emissionMT = pd.DataFrame(bestmodel.means_,columns=self.label,index=self.state)
        transitionMT = pd.DataFrame(bestmodel.transmat_,columns=self.state,index=self.state)

        if outtype == "predict":
            return(predictMT)
        elif outtype == "emission":
            return(emissionMT)
        elif outtype == "transition":
            return(transitionMT)

    def plotOSMM(self,outtype="predict",plotHiC=False,plotHiC_para=[]):
        mt = self.oneSampleMultiMetric(outtype)
        if outtype == "predict":
            plt.scatter(self.index,mt,c=mt,marker="8")
            plt.yticks(range(self.ncluster),self.state)

        elif outtype == "emission" or outtype == "transition":
            sns.heatmap(mt,cmap="coolwarm")

    def plotHiC(self,mode,path,resolution,startSite,endSite,control_path=""):
        mt = self.oneSampleMultiMetric("predict")
        plt.figure(figsize=(10,10))
        plt.subplot2grid((6,11),(0,0),rowspan=5,colspan=10)
        if mode == "single":
            hicplot = PlotTri(path,resolution,startSite,endSite)
            hicplot.draw()
        elif mode == "differ":
            hicplot = DiffDraw(path,control_path,resolution,startSite,endSite)
            hicplot.draw_tri()
        plt.subplot2grid((6,11),(5,0),rowspan=1,colspan=11)
        plt.scatter(self.index,mt,c=mt,marker="8")
        plt.yticks(range(self.ncluster),self.state)

    def oneSampleOneMetric(self):
        pass
    def multiSampleOneMetric(self):
        pass
    def multiSampleMultiMetric(self):
        pass
