from calculateMetrics import *

def getMultiSamplesScore(sampleList, labels, res, chr, mode, UniqueParameter):
    if mode == 'IS':
        for i,path in enumerate(sampleList):
            if i==0: metricMT = InsulationScore(path,res,chr,square_size=UniqueParameter).getIS()
            else:
                next = InsulationScore(path,res,chr,square_size=UniqueParameter).getIS().iloc[:,3:4]
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

    metricMT.index = metricMT.start.tolist()
    metricMT = metricMT.iloc[:,3:]
    metricMT.columns = labels
    return metricMT
