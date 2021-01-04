from MultiTypeScore import *
from calculateMetrics import *

def getDiscrete(path,res,chr,mode,parameter,control_path=""):
    if mode == "PC1":  #Acompartment~1,Bcompartment~-1
        ob = multiScore(path,res,chr)
        score = ob.obtainOneScore(mode,parameter)
        state = np.zeros(score.shape[0])*np.NaN
        state[score.iloc[:,3] > 0] = 1
        state[score.iloc[:,3] < 0] = -1
        score.iloc[:,3] =state

    elif mode in ["deltaDLR","ISC","CIC","intraSC","interSC","DRF"]:
        #positive-decompaction: 1; negative-compaction:-1
        # ISC: positive-more interaction: 1; negative-less interaction:-1
        # CIC: postive 1 stronger boundary. negative -1 weaker boundart
        ob = multiScore(path,res,chr,control_path=control_path)
        score = ob.obtainTwoScore(mode,parameter)
        state = np.array(["none"]*score.shape[0])
        state[score.iloc[:,3] > 0] = "up"+mode
        state[score.iloc[:,3] < 0] = "down"+mode
        score.iloc[:,3] =state

    elif mode == "CorrD":
        ob = multiScore(path,res,chr,control_path=control_path)
        score = ob.obtainTwoScore(mode,parameter)
        thresh = score.CorrD.describe()[5]
        state = np.array(["none"]*score.shape[0])
        state[score.iloc[:,3] >= thresh] = str("highCorr")
        state[score.iloc[:,3] < thresh] = str("lowCorr")
        score.iloc[:,3] =state

    elif mode == "PC1C":
        ob = multiScore(path,res,chr,control_path=control_path)
        score = ob.obtainTwoScore(mode,parameter)
        state = np.array(["none"]*score.shape[0])
        state[score.iloc[:,3] > 0] = "BtoA"
        state[score.iloc[:,3] < 0] = "AtoB"
        score.iloc[:,3] =state

    elif mode == "border":
        tad = TADcallIS(path,res,chr)
        bd = np.concatenate([np.array(tad.TADstart),np.array(tad.TADend)])
        bd = np.unique(bd)
        IS = multiScore(path,res,chr).obtainOneScore("IS",parameter)

        state = np.array(["nonBorder"]*IS.shape[0])
        for i in bd:
            state[IS.start == i] = "border"
            state[IS.start == i-res] = "border"
            state[IS.start == i+res] = "border"
        score = IS
        score.iloc[:,3] =state
        score.columns=["chr","start","end","TADborder"]

    return(score)

class multiTypeDiscrete:
    def __init__(self,path,control_path,res,chr,
                typelist=["PC1","border","deltaDLR","ISC","CIC","intraSC","interSC","DRF","CorrD","PC1C"],
                parameterlist=["NoDefault",300000,3000000,300000,300000,300000,300000,[200000,5000000],"pearson","NoDefault"]):
        self.path = path
        self.res = res
        self.chr = chr
        self.control_path = control_path
        self.typelist = typelist
        self.parameterlist = parameterlist

    def multiDiscrete(self):
        for i,type in enumerate(self.typelist):
            if i == 0:
                mt = getDiscrete(self.path,self.res,self.chr,self.typelist[i],
                                self.parameterlist[i],control_path=self.control_path)
            else:
                next = getDiscrete(self.path,self.res,self.chr,self.typelist[i],
                                self.parameterlist[i],control_path=self.control_path).iloc[:,3]
                mt = pd.concat([mt,next],axis=1)

        return(mt)


class multiSampleDiscrete:
    def __init__(self,pathlist,namelist,res,chr,mode,UniqueParameter):
        self.pathlist = pathlist
        self.namelist = namelist
        self.res = res
        self.chr = chr
        self.mode = mode
        self.UniqueParameter = UniqueParameter
        self.nScore = len(namelist)

    def getMultiDiscrete(self):
        for i,path in enumerate(self.pathlist):
            if i==0: metricMT = getDiscrete(path,self.res,self.chr,self.mode,self.UniqueParameter)
            else:
                next = getDiscrete(path,self.res,self.chr,self.mode,self.UniqueParameter).iloc[:,3:4]
                metricMT = pd.concat([metricMT,next],axis=1)

        metricMT.index = metricMT.start.tolist()
        metricMT = metricMT.iloc[:,3:]
        metricMT.columns = self.namelist
        return(metricMT)

    def plotMultiDiscrete(self,hic_path,start,end,clmax=100,heatmin=None,interpolation='none'):
        sbin = start//self.res
        ebin = end//self.res

        from callDirectionalTAD import PlotTAD
        plt.figure(figsize=(10,9+self.nScore))
        plt.subplot2grid((5+int(self.nScore/1.5),11),(0,0),rowspan=5,colspan=10)
        hp = PlotTAD(hic_path,self.res,start,end,clmax=clmax)
        hp.draw()

        plt.subplot2grid((5+int(self.nScore/1.5),11),(5,0),rowspan=(self.nScore//5)+1,colspan=11)
        df = self.getMultiDiscrete().iloc[sbin:ebin,:].T
        plt.imshow(df,aspect="auto",interpolation=interpolation,vmin=heatmin)
        plt.yticks(range(self.nScore),self.namelist)

def plot_discrete(mt,res,hic_path,start,end,clmax=100,heatmin=None):
    sbin = start//self.res
    ebin = end//self.res

    from callDirectionalTAD import PlotTAD
    plt.figure(figsize=(10,9+self.nScore*1))
    plt.subplot2grid((5+self.nScore,11),(0,0),rowspan=5,colspan=10)
    hp = PlotTAD(hic_path,self.res,start,end,clmax=clmax)
    hp.draw()

    plt.subplot2grid((5+self.nScore,11),(5,0),rowspan=self.nScore//2,colspan=11)
    df = self.score.iloc[sbin:ebin,:].T
    plt.imshow(df,aspect="auto",interpolation='none',cmap="Purples",vmin=heatmin)
