from .calculateMetrics import *
import seaborn as sns
from .loadfile import *
import matplotlib.pyplot as plt

def readIF(parameter,chr,normIF=True,custom_name="InteractionFrequency"):
    all = pd.read_csv(parameter,sep="\t",header=None)
    score = all[all[0] == chr]
    if normIF:
        beforlog = score[3].copy()
        afterlog = np.log1p(beforlog)
        score[3] = afterlog / np.mean(afterlog[afterlog>0])
    score.index = range(score.shape[0])
    score.columns = ["chr","start","end",custom_name]
    return(score)

def getMultiSamplesScore(sampleList, labels, res, chr, mode, UniqueParameter=None,smoothPC=True,logPC=False,
                        datatype="matrix",gt=None):
    if datatype == "rawhic" and not gt: raise ValueError("rawhic requires Genometable")
    if datatype == "rawhic" and mode != "IF":
        newlist = []
        for rawpath in sampleList:
            newpath = hic2matrix(rawpath,res,chr,gt)
            newlist.append(newpath)
        sampleList = newlist

    if mode == 'IS':
        if not UniqueParameter: UniqueParameter=300000
        for i,path in enumerate(sampleList):
            if i==0: metricMT = InsulationScore(path,res,chr,square_size=UniqueParameter).getIS()
            else:
                next = InsulationScore(path,res,chr,square_size=UniqueParameter).getIS().iloc[:,3:4]
                metricMT = pd.concat([metricMT,next],axis=1)

    elif mode == 'raw':
        for i,path in enumerate(sampleList):
            if i==0: metricMT = pd.DataFrame(loadDenseMatrix(path).values.flatten())
            else:
                next = pd.DataFrame(loadDenseMatrix(path).values.flatten())
                metricMT = pd.concat([metricMT,next],axis=1)

    elif mode == "DI":
        if not UniqueParameter: UniqueParameter=1000000
        for i,path in enumerate(sampleList):
            if i==0: metricMT = DirectionalityIndex(path,res,chr,DI_distance=UniqueParameter).getDI()
            else:
                next = DirectionalityIndex(path,res,chr,DI_distance=UniqueParameter).getDI().iloc[:,3:4]
                metricMT = pd.concat([metricMT,next],axis=1)

    elif mode == "CI":
        if not UniqueParameter: UniqueParameter=300000
        for i,path in enumerate(sampleList):
            if i==0: metricMT = ContrastIndex(path,res,chr,CI_size=UniqueParameter).getCI()
            else:
                next = ContrastIndex(path,res,chr,CI_size=UniqueParameter).getCI().iloc[:,3:4]
                metricMT = pd.concat([metricMT,next],axis=1)

    elif mode == "SS":
        if not UniqueParameter: UniqueParameter=300000
        for i,path in enumerate(sampleList):
            if i==0: metricMT = SeparationScore(path,res,chr,TADss_size=UniqueParameter).getTADss()
            else:
                next = SeparationScore(path,res,chr,TADss_size=UniqueParameter).getTADss().iloc[:,3:4]
                metricMT = pd.concat([metricMT,next],axis=1)

    elif mode == "DLR":
        if not UniqueParameter: UniqueParameter=3000000
        for i,path in enumerate(sampleList):
            if i==0: metricMT = DistalToLocalRatio(path,res,chr,sizeDLR=UniqueParameter).getDLR()
            else:
                next = DistalToLocalRatio(path,res,chr,sizeDLR=UniqueParameter).getDLR().iloc[:,3:4]
                metricMT = pd.concat([metricMT,next],axis=1)

    elif mode == "PC1":
        if not UniqueParameter:
            warnings.warn("The sign of eigenvector is arbitrary unless specify a geneDensity file")
        for i,path in enumerate(sampleList):
            if i==0: metricMT = CompartmentPC1(path,res,chr).getPC1(signCorr = UniqueParameter)
            else:
                next = CompartmentPC1(path,res,chr).getPC1(signCorr = UniqueParameter,smooth = smoothPC, logOE=logPC).iloc[:,3:4]
                metricMT = pd.concat([metricMT,next],axis=1)

    elif mode == "IAS":
        if not UniqueParameter: UniqueParameter=300000
        for i,path in enumerate(sampleList):
            if i==0: metricMT = intraTADscore(path,res,chr).getIntraS(IS_size = UniqueParameter)
            else:
                next = intraTADscore(path,res,chr).getIntraS(IS_size = UniqueParameter).iloc[:,3:4]
                metricMT = pd.concat([metricMT,next],axis=1)

    elif mode == "IAS":
        if not UniqueParameter: UniqueParameter=300000
        for i,path in enumerate(sampleList):
            if i==0: metricMT = interTADscore(path,res,chr).getInterS(IS_size = UniqueParameter)
            else:
                next = interTADscore(path,res,chr).getInterS(IS_size = UniqueParameter).iloc[:,3:4]
                metricMT = pd.concat([metricMT,next],axis=1)

    elif mode == "IF":
        for i,path in enumerate(sampleList):
            if i==0: metricMT = InteractionFrequency(path,res,chr,gt=gt,datatype=datatype).getIF()
            else:
                next = InteractionFrequency(path,res,chr,gt=gt,datatype=datatype).getIF().iloc[:,3:4]
                metricMT = pd.concat([metricMT,next],axis=1)

    if mode != "raw":
        #metricMT.index = metricMT.start.tolist()
        #metricMT = metricMT.iloc[:,3:]
        metricMT.columns = ["chr","start","end"]+labels
    elif mode == "raw":
        metricMT.columns = labels
    return metricMT

class repQC:
    def __init__(self,pathlist,namelist,res,chr,mode,UniqueParameter,method="pearson",smoothPC=True,logPC=False,gt=None,datatype="matrix"):
        self.pathlist = pathlist
        self.namelist = namelist
        self.res = res
        self.chr = chr
        self.mode = mode
        self.UniqueParameter = UniqueParameter
        self.nScore = len(namelist)

        self.score = getMultiSamplesScore(self.pathlist,namelist,res=res,chr=chr,mode=mode,
                                        UniqueParameter=UniqueParameter,smoothPC=smoothPC,logPC=logPC,
                                        gt=gt,datatype=datatype)
        self.corrMT = self.score.iloc[:,3:].corr(method=method)

    def corr_plot(self):
        if self.mode!="raw":
            sns.clustermap(self.corrMT ,cmap="RdPu")

    def calcuRepScore(self):
        maxCorr = abs(np.nan_to_num(self.corrMT[self.corrMT<1])).max()
        minCorr = abs(np.array(self.corrMT)).min()
        return(maxCorr-minCorr)

    def heatmap(self,start,end,figs=(10,5),vmin=None,vmax=None):
        sbin = start//self.res
        ebin = end//self.res
        plt.figure(figsize=figs)
        plt.imshow(self.score.iloc[sbin:ebin,:].T,aspect="auto",interpolation='nearest',vmin=vmin,vmax=vmax)
        plt.colorbar()

    def heatmap_tri(self,hic_path,start,end,clmax=100,heatmin=None,plottype="heat"):
        sbin = start//self.res
        ebin = end//self.res
        if ebin ==0: ebin = self.score.shape[0]

        from .callDirectionalTAD import PlotTAD
        plt.figure(figsize=(10,9+self.nScore*1))
        plt.subplot2grid((5+self.nScore,11),(0,0),rowspan=5,colspan=10)
        hp = PlotTAD(hic_path,self.res,self.chr,start,end,clmax=clmax)
        hp.draw()

        plt.subplot2grid((5+self.nScore,11),(5,0),rowspan=self.nScore//2,colspan=11)
        df = self.score.iloc[sbin:ebin,3:].T
        if plottype == "heat":
            plt.imshow(df,aspect="auto",interpolation='none',cmap="Purples",vmin=heatmin)
            plt.yticks(range(len(self.namelist)),self.namelist)
        elif plottype == "line":
            df.columns = range(df.shape[1])
            print(df.T)
            df.T.plot(ax=plt.gca())
            plt.xlim(0,df.shape[1]-1)
            plt.legend()
        plt.xticks(np.arange(0,df.shape[1]+1,(df.shape[1])/5),hp.mark)

    def anova_like(self,start,end):
        sbin = start//self.res
        ebin = end//self.res
        df = self.score.iloc[sbin:ebin,:].T

        nLoci = df.shape[1]
        nScore = df.shape[0]
        arrays = np.zeros(nLoci)*np.NaN
        for i in range(nLoci):
            if (i-2)<0 or (i+2)>=nLoci: continue
            df_i = df.iloc[:,i-2:i+3]
            df_ilist = [df_i.iloc[j] for j in range(nScore)]
            pvalue = stats.f_oneway(*df_ilist)[1]
            arrays[i] = pvalue

        return(arrays)

        #from statsmodels.sandbox.stats.multicomp import multipletests
        #qvalue=multipletests(arrays, method='fdr_bh')
