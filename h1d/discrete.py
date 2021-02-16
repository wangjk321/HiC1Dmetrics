from MultiTypeScore import *
from calculateMetrics import *
import os

def getDiscrete(path,res,chr,mode,parameter,control_path="",label=True):
    if mode == "PC1":  #Acompartment~1,Bcompartment~-1
        ob = multiScore(path,res,chr)
        score = ob.obtainOneScore(mode,parameter)
        if label:
            state = np.array(["NA"]*score.shape[0],dtype=object)
            state[score.iloc[:,3] > 0] = "CompartA"
            state[score.iloc[:,3] < 0] = "CompartB"
        else:
            state = np.array([np.NaN]*score.shape[0])
            state[score.iloc[:,3] > 0] = 1
            state[score.iloc[:,3] < 0] = 0
        score.iloc[:,3] =state

    elif mode in ["deltaDLR"]:
        #positive-decompaction: 1; negative-compaction:-1
        ob = multiScore(path,res,chr,control_path=control_path)
        score = ob.obtainTwoScore(mode,parameter)
        state = np.array(["NA"]*score.shape[0],dtype=object)
        state[score.iloc[:,3] > 0] = "decompact"
        state[score.iloc[:,3] < 0] = "compact"
        score.iloc[:,3] =state

    elif mode in ["ISC","CIC","intraSC","interSC","DRF"]:
        # ISC: positive-more interaction: 1; negative-less interaction:-1
        # CIC: postive 1 stronger boundary. negative -1 weaker boundart
        ob = multiScore(path,res,chr,control_path=control_path)
        score = ob.obtainTwoScore(mode,parameter)
        state = np.array(["NA"]*score.shape[0],dtype=object)
        state[score.iloc[:,3] > 0] = "up"+mode
        state[score.iloc[:,3] < 0] = "down"+mode
        score.iloc[:,3] =state

    elif mode == "CorrD":
        ob = multiScore(path,res,chr,control_path=control_path)
        score = ob.obtainTwoScore(mode,parameter)
        thresh = score.CorrD.describe()[5]
        state = np.array(["NA"]*score.shape[0],dtype=object)
        state[score.iloc[:,3] >= thresh] = str("highCorr")
        state[score.iloc[:,3] < thresh] = str("lowCorr")
        score.iloc[:,3] =state

    elif mode == "PC1C":
        ob = multiScore(path,res,chr,control_path=control_path)
        score = ob.obtainTwoScore(mode,parameter)
        state = np.array(["none"]*score.shape[0],dtype=object)
        state[score.iloc[:,3] > 0] = "BtoA"
        state[score.iloc[:,3] < 0] = "AtoB"
        score.iloc[:,3] =state

    elif mode == "border":
        tad = TADcallIS(path,res,chr)
        bd = np.concatenate([np.array(tad.TADstart),np.array(tad.TADend)])
        bd = np.unique(bd)
        IS = multiScore(path,res,chr).obtainOneScore("IS",parameter)

        state = np.array(["nonBorder"]*IS.shape[0],dtype=object)
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
            print(i,type)
            if i == 0:
                mt = getDiscrete(self.path,self.res,self.chr,self.typelist[i],
                                self.parameterlist[i],control_path=self.control_path)
            else:
                next = getDiscrete(self.path,self.res,self.chr,self.typelist[i],
                                self.parameterlist[i],control_path=self.control_path).iloc[:,3]
                mt = pd.concat([mt,next],axis=1)

        return(mt)

    def makecsv(self,outname="multiDiscrete.txt",s=450,e=650):
        self.multiDiscrete().iloc[s:e,3:].to_csv(outname,sep="\t",index=None)

    def useHMM(self):
        self.makecsv()
        import rpy2.robjects as robjects
        r = robjects.r
        if not os.path.exists("seqHMM"):
            os.mkdir("seqHMM")

        r("library(seqHMM)")
        r("mt <- read.csv('multiDiscrete.txt',sep='\t')")
        r("mt[mt=='NA'] <- NA")

        r("mtlist <- list()")
        r("for(i in seq(1,dim(mt)[2])){mtlist[[i]] <- seqdef(t(as.matrix(mt[,i])))}")

        r("pdf('seqHMM/observedstate.pdf',width=10,height = 8)")
        r("ssplot(mtlist, title = 'Observed states plots')")
        r("dev.off()")

        r("mod <- build_hmm(observations = mtlist, n_states = 5)")
        r("hmm <- fit_model(mod,global_step = TRUE, local_step = TRUE,threads=10,control_em = list(restart=list(times = 100)))")

        r("tMT <- hmm$model$transition_probs")
        r("write.table(tMT,'seqHMM/transitionMatrix.txt',sep='\t')")

        r("hMT <- hidden_paths(hmm$model)")
        r("write.table(hMT,'seqHMM/hiddenState.txt',sep='\t')")

        r("emp<- hmm$model$emission_probs")
        r("eMT <- matrix(data = NA,dim(emp[[1]])[1],length(emp))")
        r("label <- c()")
        r("for(i in seq(1,length(emp))){eMT[,i] <- emp[[i]][,1];label[i] <- colnames(emp[[i]])[1]}")
        r("colnames(eMT) <- label")
        r("row.names(eMT) <- row.names(emp[[1]])")
        r("write.table(eMT,'seqHMM/emissionMatrix.txt',sep='\t')")

        hMT = pd.read_csv("seqHMM/hiddenState.txt",sep="\t")
        eMT = pd.read_csv("seqHMM/emissionMatrix.txt",sep="\t")
        tMT = pd.read_csv("seqHMM/transitionMatrix.txt",sep="\t")

        return(hMT,eMT,tMT)

    def plotHMM(self,type="hidden",plotHiC="",s=None,e=None,clmax=100):
        hMT,eMT,tMT = self.useHMM()
        if type == "hidden":
            from sklearn.preprocessing import LabelEncoder
            hmt = hMT.iloc[0,:]
            le = LabelEncoder().fit(hmt)
            hmtN = le.transform(hmt)
            if plotHiC in ["single","differ"]:
                plt.figure(figsize=(10,10))
                plt.subplot2grid((6,11),(0,0),rowspan=5,colspan=10)
                if plotHiC == "single":
                    hicplot = PlotTri(self.path,self.res,startSite=s,endSite=e,clmax=clmax)
                    hicplot.draw()
                elif plotHiC == "differ":
                    hicplot = DiffDraw(self.path,self.control_path,self.res,startSite=s,endSite=e)
                    hicplot.draw_tri()
                plt.subplot2grid((6,11),(5,0),rowspan=1,colspan=11)
            plt.scatter(range(len(hmt)),hmtN,c=hmtN,marker="8")
            plt.yticks(range(5),list(le.classes_))
        elif type == "transition":
            sns.heatmap(tMT,cmap="coolwarm")
        elif type == "emission":
            sns.heatmap(eMT.T,cmap="coolwarm")


class multiSampleDiscrete:
    def __init__(self,pathlist,namelist,res,chr,mode,UniqueParameter):
        self.pathlist = pathlist
        self.namelist = namelist
        self.res = res
        self.chr = chr
        self.mode = mode
        self.UniqueParameter = UniqueParameter
        self.nScore = len(namelist)

    def getMultiDiscrete(self,label=True):
        for i,path in enumerate(self.pathlist):
            if i==0: metricMT = getDiscrete(path,self.res,self.chr,self.mode,self.UniqueParameter,label=label)
            else:
                next = getDiscrete(path,self.res,self.chr,self.mode,self.UniqueParameter,label=label).iloc[:,3:4]
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
        plt.subplot2grid((5+int(self.nScore/1.2),11),(0,0),rowspan=5,colspan=10)
        hp = PlotTAD(hic_path,self.res,self.chr,start,end,clmax=clmax)
        hp.draw()

        plt.subplot2grid((5+int(self.nScore/1.2),11),(5,0),rowspan=(self.nScore//7)+1,colspan=11)
        df = self.getMultiDiscrete(label=False).iloc[sbin:ebin,:].T
        plt.imshow(df,aspect="auto",interpolation=interpolation,vmin=heatmin)
        plt.yticks(range(self.nScore),self.namelist)
        plt.xticks(np.arange(0,df.shape[1]+1,(df.shape[1])/5),hp.mark)
