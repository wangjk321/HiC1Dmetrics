from calculateMetrics import *
import seaborn as sns
from loadfile import *

def getMultiSamplesScore(sampleList, labels, res, chr, mode, UniqueParameter):
    if mode == 'IS':
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
        for i,path in enumerate(sampleList):
            if i==0: metricMT = DirectionalityIndex(path,res,chr,DI_distance=UniqueParameter).getDI()
            else:
                next = DirectionalityIndex(path,res,chr,DI_distance=UniqueParameter).getDI().iloc[:,3:4]
                metricMT = pd.concat([metricMT,next],axis=1)

    elif mode == "CI":
        for i,path in enumerate(sampleList):
            if i==0: metricMT = ContrastIndex(path,res,chr,CI_size=UniqueParameter).getCI()
            else:
                next = ContrastIndex(path,res,chr,CI_size=UniqueParameter).getCI().iloc[:,3:4]
                metricMT = pd.concat([metricMT,next],axis=1)

    elif mode == "TADss":
        for i,path in enumerate(sampleList):
            if i==0: metricMT = SeparationScore(path,res,chr,TADss_size=UniqueParameter).getTADss()
            else:
                next = SeparationScore(path,res,chr,TADss_size=UniqueParameter).getTADss().iloc[:,3:4]
                metricMT = pd.concat([metricMT,next],axis=1)

    elif mode == "DLR":
        for i,path in enumerate(sampleList):
            if i==0: metricMT = DistalToLocalRatio(path,res,chr,sizeDLR=UniqueParameter).getDLR()
            else:
                next = DistalToLocalRatio(path,res,chr,sizeDLR=UniqueParameter).getDLR().iloc[:,3:4]
                metricMT = pd.concat([metricMT,next],axis=1)

    elif mode == "PC1":
        for i,path in enumerate(sampleList):
            if i==0: metricMT = CompartmentPC1(path,res,chr).getPC1(signCorr = UniqueParameter)
            else:
                next = CompartmentPC1(path,res,chr).getPC1(signCorr = UniqueParameter).iloc[:,3:4]
                metricMT = pd.concat([metricMT,next],axis=1)

    elif mode == "intraScore":
        for i,path in enumerate(sampleList):
            if i==0: metricMT = intraTADscore(path,res,chr).getIntraS(IS_size = UniqueParameter)
            else:
                next = intraTADscore(path,res,chr).getIntraS(IS_size = UniqueParameter).iloc[:,3:4]
                metricMT = pd.concat([metricMT,next],axis=1)

    elif mode == "interScore":
        for i,path in enumerate(sampleList):
            if i==0: metricMT = interTADscore(path,res,chr).getInterS(IS_size = UniqueParameter)
            else:
                next = interTADscore(path,res,chr).getInterS(IS_size = UniqueParameter).iloc[:,3:4]
                metricMT = pd.concat([metricMT,next],axis=1)

    if mode != "raw":
        metricMT.index = metricMT.start.tolist()
        metricMT = metricMT.iloc[:,3:]
    metricMT.columns = labels
    return metricMT

class repQC:
    def __init__(self,pathlist,namelist,res,chr,mode,UniqueParameter,method="pearson"):
        self.pathlist = pathlist
        self.namelist = namelist
        self.res = res
        self.chr = chr
        self.mode = mode
        self.UniqueParameter = UniqueParameter

        score = getMultiSamplesScore(self.pathlist,namelist,res=res,chr=chr,mode=mode,UniqueParameter=UniqueParameter)
        self.corrMT = score.corr(method=method)

    def corr_plot(self):
        sns.clustermap(self.corrMT ,cmap="RdPu")

    def calcuRepScore(self):
        maxCorr = abs(np.nan_to_num(self.corrMT[self.corrMT<1])).max()
        minCorr = abs(np.array(self.corrMT)).min()
        return(maxCorr-minCorr)
