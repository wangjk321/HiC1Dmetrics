from MultiTypeScore import *
from calculateMetrics import *

def getDiscrete(path,res,chr,mode,parameter,control_path=""):
    ob = multiScore(path,res,chr)
    score = ob.obtainOneScore(mode,parameter)
    state = np.zeros(score.shape[0])*np.NaN
    if mode == "PC1":
        state[score.iloc[:,3] > 0] = 1
        state[score.iloc[:,3] < 0] = -1
    elif mode == "border":
        tad = TADcallIS()

    score.iloc[:,3] =state

    return(score)

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

        plt.subplot2grid((5+int(self.nScore/1.5),11),(5,0),rowspan=(self.nScore//4)+1,colspan=11)
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
