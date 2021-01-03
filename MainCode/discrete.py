from MultiTypeScore import *

def getDiscrete(path,res,chr,mode,parameter):
    ob = multiScore(path,res,chr)
    score = ob.obtainOneScore(mode,parameter)
    state = np.zeros(score.shape[0])*np.NaN
    state[score.iloc[:,3] > 0] = 1
    state[score.iloc[:,3] < 0] = -1
    score.iloc[:,3] =state

    return(score)

def getMultiDiscrete(sampleList, labels, res, chr, mode, UniqueParameter):
    for i,path in enumerate(sampleList):
        if i==0: metricMT = getDiscrete(path,res,chr,mode,UniqueParameter)
        else:
            next = getDiscrete(path,res,chr,mode,UniqueParameter).iloc[:,3:4]
            metricMT = pd.concat([metricMT,next],axis=1)

    if mode != "raw":
        metricMT.index = metricMT.start.tolist()
        metricMT = metricMT.iloc[:,3:]
    metricMT.columns = labels
    return metricMT
