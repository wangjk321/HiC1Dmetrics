from multiprocessing import Pool
from callDirectionalTAD import *
from MultiTypeScore import *
#from multiprocessing.pool import ThreadPool as Pool #多线程

# calculate each dTAD of chromosome separately
class paralfunc(object):
    def __init__(self,myfun,num_processer):
        self.myfun = myfun
        self.num_processer = num_processer

        #chrlist= list(range(1,23)); chrlist.append("X")
        #chrolist = ["chr"+str(a) for a in chrlist]
        #self.chrolist = chrolist

    def run(self,pathName,pathControl,resolution,maxchr=22):
        chrlist= list(range(1,maxchr+1)); chrlist.append("X")
        chrolist = ["chr"+str(a) for a in chrlist]
        self.chrolist = chrolist

        resultlist=[]
        p = Pool(self.num_processer)
        for r in self.chrolist:
            result = p.apply_async(self.myfun, args=(r,pathName,pathControl,resolution,))
            resultlist.append(result)
        p.close()
        p.join()

        output_left = pd.DataFrame()
        output_right = pd.DataFrame()
        for i in resultlist:
            output_single = i.get().copy()
            output_left = output_left.append(output_single["leftTAD"])
            output_right = output_right.append(output_single["rightTAD"])
        return(output_left,output_right)

def call_dTAD(chrom,pathName,pathControl,resolution):
    filename = pathName+"/observed.KR."+chrom+".matrix.gz"
    controlfilename = pathControl + "/observed.KR."+chrom+".matrix.gz"

    dTAD = DirectionalTAD(filename,controlfilename, \
                      resolution,chr=chrom,startDRF=500000,sizeDRF=2000000,\
                      sizeIS=300000)
    leftTAD,rightTAD,_ = dTAD.extractRegion()
    return({"leftTAD":leftTAD,"rightTAD":rightTAD})

def alldTAD(pathName,pathControl,resolution,outname="twoSample",maxchr=22):
    ltad,rtad = paralfunc(call_dTAD,30).run(pathName,pathControl,resolution,maxchr)
    ltad.to_csv(outname+"_leftTAD.csv",sep="\t",index=False)
    rtad.to_csv(outname+"_rightTAD.csv",sep="\t",index=False)

#alldTAD("/Users/wangjiankang/figureServer/Nov2020/Rad21KD1_HiCmatrix",
#        "/Users/wangjiankang/figureServer/Nov2020/Control1_HiCmatrix",50000,maxchr=23)

#=====================================#
#simply call all TAD

class paralfuncOneSample(object):
    def __init__(self,myfun,num_processer):
        self.myfun = myfun
        self.num_processer = num_processer

    def run(self,pathName,resolution,maxchr=22,type=None,parameter=None):
        chrlist= list(range(1,maxchr+1)); chrlist.append("X")
        chrolist = ["chr"+str(a) for a in chrlist]
        self.chrolist = chrolist

        resultlist=[]
        p = Pool(self.num_processer)
        for r in self.chrolist:
            result = p.apply_async(self.myfun, args=(r,pathName,resolution,type,parameter))
            resultlist.append(result)
        p.close()
        p.join()

        output_all = pd.DataFrame()
        for i in resultlist:
            output_single = i.get().copy()
            output_all = output_all.append(output_single)
        return(output_all)

def TAD1sample(chrom,pathName,resolution,type,parameter):
    filename = pathName+"/observed.KR."+chrom+".matrix.gz"
    tad = TADcallIS(filename,resolution,chrom,squareSize=300000)
    return(tad)

def runTAD1sample(pathName,resolution,outname="allTAD",maxchr=22):
    tad = paralfuncOneSample(TAD1sample,30).run(pathName,resolution,maxchr)
    tad.to_csv(outname+"_TAD.csv",sep="\t",index=False)

def oneScoreSinglechr(chrom,pathName,resolution,type,parameter):
    filename = pathName+"/observed.KR."+chrom+".matrix.gz"
    score = multiScore(filename,resolution,chrom).obtainOneScore(type,parameter)
    return(score)

def oneScoreAllchr(pathName,resolution,type,parameter,outname="OneScore",maxchr=22):
    allscore = paralfuncOneSample(oneScoreSinglechr,30).run(pathName,resolution,maxchr,type,parameter)
    return(allscore)
    #allscore.to_csv(outname+"_"+type+".csv",sep="\t",index=False)

def multiScoreSinglechr(chrom,pathName,resolution,typelist,parameterlist):
    filename = pathName+"/observed.KR."+chrom+".matrix.gz"
    multiType = multiScore(filename,resolution,chrom).allOneScore(typelist,parameterlist)
    return(multiType)

def multiScoreAllchr(pathName,resolution,typelist,parameterlist,outname="multiScore",maxchr=22):
    allmultiscore = paralfuncOneSample(multiScoreSinglechr,30).run(pathName,resolution,maxchr,typelist,parameterlist)
    return(allmultiscore)

#runTAD1sample("/Users/wangjiankang/figureServer/Nov2020/Rad21KD1_HiCmatrix",50000,maxchr=5)
