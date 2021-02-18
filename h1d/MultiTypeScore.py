from .calculateMetrics import *
from .calculateTwoSample import *
#from hmmlearn import hmm
import matplotlib.pyplot as plt
import seaborn as sns
from .plotTwoSample import *
#from plotMetrics import *
from scipy import stats
import os
import subprocess
#from callDirectionalTAD import *

class multiScore:
    def __init__(self,path,res,chr,control_path=""):
        self.rawpath = path # store for IF
        self.path = path
        self.res = res
        self.chr = chr
        self.rawcontrol = control_path
        self.control_path = control_path



    def obtainOneScore(self,mode,parameter=None,smoothPC=True,logPC=False,
                        custom_name="InteractionFrequency",normIF=True,gt=None,datatype="matrix"):
        if datatype == "rawhic" and not gt: raise ValueError("rawhic requires Genometable")

        if datatype == "rawhic" and mode != "IF":
            self.path = hic2matrix(self.path,self.res,self.chr,gt)

        if mode == "IS":
            if not parameter: parameter=300000
            score = InsulationScore(self.path,self.res,self.chr,square_size=int(parameter)).getIS()
        elif mode == "CI":
            if not parameter: parameter=300000
            score = ContrastIndex(self.path,self.res,self.chr,CI_size=int(parameter)).getCI()
        elif mode == "DI":
            if not parameter: parameter=1000000
            score = DirectionalityIndex(self.path,self.res,self.chr,DI_distance=int(parameter)).getDI()
        elif mode == "SS":
            if not parameter: parameter=300000
            score = SeparationScore(self.path,self.res,self.chr,TADss_size=int(parameter)).getTADss()
        elif mode == "DLR":
            if not parameter: parameter=3000000
            score = DistalToLocalRatio(self.path,self.res,self.chr,sizeDLR=int(parameter)).getDLR()
        elif mode == "intraS":
            if not parameter: parameter=300000
            score = intraTADscore(self.path,self.res,self.chr).getIntraS(IS_size = int(parameter))
        elif mode == "interS":
            if not parameter: parameter=300000
            score = interTADscore(self.path,self.res,self.chr).getInterS(IS_size = int(parameter))
        elif mode == "PC1":
            if not parameter:
                warnings.warn("The sign of eigenvector is arbitrary unless specify a geneDensity file")
            score = CompartmentPC1(self.path,self.res,self.chr).getPC1(signCorr = parameter,smooth = smoothPC, logOE=logPC)
        elif mode == "IF":
            if not gt:
                raise ValueError("Genometable is required for the calculation of IF")
            codepath = os.path.dirname(os.path.realpath(__file__))
            soft = codepath+"/InteractionFreq.sh"
            juicer = codepath+"/jc/jctool_1.11.04.jar"
            chrnum = self.chr.split("chr")[1]
            os.system("sh '"+soft+"' '"+juicer+"' "+self.path+" "+chrnum+" "+str(self.res)+" "+gt+" "+"IF_"+self.chr) #in case of space
            score = pd.read_csv("IF_"+self.chr+".bedGraph",sep="\t",header=None)
            if normIF:
                beforlog = score[3].copy()
                afterlog = np.log1p(beforlog)
                score[3] = afterlog / np.mean(afterlog[afterlog>0])
            score.index = range(score.shape[0])
            score.columns = ["chr","start","end",custom_name]
            os.system("rm "+"IF_"+self.chr+".bedGraph")

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
                    smoothPC=True,logPC=False,normIF=True,datatype="matrix",gt=None):
        for i,type in enumerate(typelist):
            if i == 0:
                multiType = self.obtainOneScore(mode=typelist[i],parameter=parameterlist[i],smoothPC=smoothPC,logPC=logPC,
                                                datatype=datatype,gt=gt)
            else:
                if datatype=="rawhic" and self.path !=self.rawpath:
                    datatype2="matrix"
                else: datatype2 =datatype
                if typelist[i] == "IF":
                    self.path = self.rawpath; datatype2="rawhic"
                next = self.obtainOneScore(mode=typelist[i],parameter=parameterlist[i],smoothPC=smoothPC,logPC=logPC,
                                        datatype=datatype2,gt=gt).iloc[:,3]
                multiType = pd.concat([multiType,next],axis=1)

        return(multiType)

    def plotOneScore(self,typelist=["IS","CI","DI","TADss","DLR","intraS","interS","PC1","custom"],
                    parameterlist=[300000,300000,1000000,300000,3000000,300000,300000,"NotSpecified","customPath"],
                    start=0,end=0,clmax=100,plotTAD=False,smoothPC=True,logPC=False,datatype="matrix",gt=None):
        import matplotlib.colors as mcolors
        from .callDirectionalTAD import PlotTAD

        cols = list(mcolors.TABLEAU_COLORS.keys())
        scoreMT = self.allOneScore(typelist,parameterlist,smoothPC,logPC,datatype=datatype,gt=gt)
        nScore = len(typelist)

        plt.figure(figsize=(10,9+nScore*1))
        plt.subplot2grid((5+nScore,11),(0,0),rowspan=5,colspan=10)

        if datatype=="rawhic" and typelist[-1] == "IF":
            self.path = hic2matrix(self.path,self.res,self.chr,gt)
        hp = PlotTAD(self.path,self.res,self.chr,start,end,clmax=clmax)
        if plotTAD == True:
            hp.drawTAD()
        elif plotTAD == False:
            hp.draw()
        for i in range(nScore):
            scoreRegion = scoreMT.iloc[start//self.res:end//self.res,3+i]
            plt.subplot2grid((5+nScore,11),(5+i,0),rowspan=1,colspan=11)
            plt.plot(scoreRegion,c=cols[i],label=scoreRegion.name)
            plt.legend(loc="upper left")
            plt.xlim(start//self.res,end//self.res)
            ticks_pos = np.arange(hp.sbin,hp.ebin+1,(hp.ebin-hp.sbin)/5)
            if i < nScore-1:
                plt.xticks(ticks_pos,[])
            else:
                plt.xticks(ticks_pos,hp.mark)

        return(scoreMT)

    def obtainTwoScore(self,mode,parameter,smoothPC=True,logPC=False,normIF=True,gt=None,datatype="matrix"):
        if datatype == "rawhic" and mode != "IFC":
            self.path = hic2matrix(self.path,self.res,self.chr,gt)
            self.control_path = hic2matrix(self.control_path,self.res,self.chr,gt)

        if mode == "ISC":
            if not parameter: parameter=300000
            score = TADScoreChange(self.path,self.control_path,self.res,self.chr).getChange("IS",int(parameter))
        elif mode == "CIC":
            if not parameter: parameter=300000
            score = TADScoreChange(self.path,self.control_path,self.res,self.chr).getChange("CI",int(parameter))
        elif mode == "DIC":
            if not parameter: parameter=1000000
            score = TADScoreChange(self.path,self.control_path,self.res,self.chr).getChange("DI",int(parameter))
        elif mode == "SSC":
            if not parameter: parameter=300000
            score = TADScoreChange(self.path,self.control_path,self.res,self.chr).getChange("TADss",int(parameter))
        elif mode == "deltaDLR":
            if not parameter: parameter=3000000
            score = deltaDLR(self.path,self.control_path,self.res,self.chr,sizeDLR=int(parameter)).getDeltaDLR()
        elif mode == "IASC":
            if not parameter: parameter=300000
            score = intraScoreChange(self.path,self.control_path,self.res,self.chr,IS_size=int(parameter)).getIntraSC()
        elif mode == "IESC":
            if not parameter: parameter=300000
            score = interScoreChange(self.path,self.control_path,self.res,self.chr,IS_size=int(parameter)).getInterSC()
        elif mode == "DRF":
            if not parameter: parameter=[200000,5000000]
            score = DirectionalRelativeFreq(self.path,self.control_path,self.res,self.chr,
                                            start_distance=parameter[0], end_distance=parameter[1]).getDRF()
        elif mode == "CD":
            if not parameter: parameter="pearson"
            score = CorrelationDifference(self.path,self.control_path,self.res,self.chr,method=parameter).getCorrD()
        elif mode == "PC1C":
            score = PC1change(self.path,self.control_path,self.res,self.chr,corr_file=parameter,smoothPC = smoothPC, logPC=logPC).getPC1change()
        elif mode == "IFC":
            if not gt: raise ValueError("Genometable is required for the calculation of IF")
            if datatype=="matrix": raise ValueError("Calculation of IF require rawhic not matrix")
            score = InteractionFrequencyChange(self.path,self.control_path,self.res,self.chr,gt=gt,datatype=datatype,normIF=normIF).getIFC()
        elif mode == "IFCback":
            scoreTreat = pd.read_csv(parameter[0],sep="\t",header=None)
            scoreControl= pd.read_csv(parameter[1],sep="\t",header=None)
            treat = scoreTreat[scoreTreat[0]==self.chr]
            control = scoreControl[scoreControl[0]==self.chr]
            treatlog = np.log1p(treat.iloc[:,3])
            controllog = np.log1p(control.iloc[:,3])
            t = treatlog / np.mean(treatlog[treatlog>0])
            c = controllog / np.mean(controllog[controllog>0])
            score = pd.DataFrame({"chr":treat[0],"start":treat[1],
                                "end":treat[2],"InteractionFrequency Change":t-c})
            score.index = range(score.shape[0])
        else: print("Error: Please specify the correct mode")

        #if datatype == "rawhic" and mode != "IF": os.system("rm -rf MatrixTemp*")
        return(score,self.path,self.control_path)

    def allTwoScore(self,typelist=["ISC","CIC","DIC","SSC","deltaDLR","intraSC","interSC","DRF","CorrD","PC1C"],
                    parameterlist=[300000,300000,1000000,300000,3000000,300000,300000,[200000,5000000],"pearson","NotSpecified"],
                    smoothPC=True,logPC=False,datatype="matrix",gt=None):
        for i,type in enumerate(typelist):
            if i == 0:
                multiType,_,_ = self.obtainTwoScore(mode=typelist[i],parameter=parameterlist[i],smoothPC=smoothPC,logPC=logPC,
                                                    datatype=datatype,gt=gt)
            else:
                if datatype=="rawhic" and self.path !=self.rawpath:
                    datatype2="matrix"
                else: datatype2 =datatype
                if typelist[i] == "IFC":
                    self.path = self.rawpath
                    self.control_path = self.rawcontrol
                    datatype2="rawhic"
                next = self.obtainTwoScore(mode=typelist[i],parameter=parameterlist[i],smoothPC=smoothPC,logPC=logPC,
                                            datatype=datatype2,gt=gt)[0].iloc[:,3]
                multiType = pd.concat([multiType,next],axis=1)

        return(multiType)

    def plotTwoScore(self,typelist=["ISC","CIC","DIC","TADssC","deltaDLR","intraSC","interSC","DRF","CorrD","PC1C"],
                    parameterlist=[300000,300000,1000000,300000,3000000,300000,300000,[200000,5000000],"pearson","NotSpecified"],
                    start=0,end=0,clmax=2,plotTAD=False,smoothPC=True,logPC=False,
                    outname="default",datatype="matrix",gt=None):
        import matplotlib.colors as mcolors
        from .callDirectionalTAD import PlotTAD

        cols = list(mcolors.TABLEAU_COLORS.keys())
        cols.remove("tab:gray")
        cols.remove("tab:pink")
        cols.append("black")

        scoreMT = self.allTwoScore(typelist,parameterlist,smoothPC,logPC,datatype=datatype,gt=gt)
        print(scoreMT)
        nScore = len(typelist)

        plt.figure(figsize=(10,9+nScore*1))
        plt.subplot2grid((5+nScore,11),(0,0),rowspan=5,colspan=10)

        if datatype=="rawhic" and typelist[-1] == "IFC":
            self.path = hic2matrix(self.path,self.res,self.chr,gt)
            self.control_path = hic2matrix(self.control_path,self.res,self.chr,gt)

        hp = DiffDraw(self.path,self.control_path,self.res,startSite=start,endSite=end,clmax=clmax)
        if plotTAD == True:
            hp.drawTAD()
        elif plotTAD == False:
            hp.draw_tri()
        for i in range(nScore):
            scoreRegion = scoreMT.iloc[start//self.res:end//self.res+1,3+i]
            plt.subplot2grid((5+nScore,11),(5+i,0),rowspan=1,colspan=11)
            plt.plot(scoreRegion,c=cols[i],label=scoreRegion.name)
            plt.legend(loc="upper left")
            plt.xlim(start//self.res,end//self.res)

            if scoreRegion.name != "CorrelationDifference":
                if scoreRegion.name == "DirectionalRelativeFrequency":
                    posColor = '#e9a3c9'; negColor = '#a1d76a'; opacity=1
                else: posColor = cols[i]; negColor = 'grey'; opacity=0.5

                plt.plot([hp.sbin,hp.ebin],[0,0],"k--",linewidth=0.6)
                plt.fill_between(np.arange(hp.sbin,hp.ebin+1,1),scoreRegion, 0,\
                                where = scoreRegion <=0,facecolor=negColor, alpha=opacity)
                plt.fill_between(np.arange(hp.sbin,hp.ebin+1,1),scoreRegion, 0,\
                                where = scoreRegion >=0,facecolor=posColor, alpha=opacity)

            ticks_pos = np.arange(hp.sbin,hp.ebin+1,(hp.ebin-hp.sbin)/5)
            if i < nScore-1:
                plt.xticks(ticks_pos,[])
            else:
                plt.xticks(ticks_pos,hp.mark)
        return(scoreMT)

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
