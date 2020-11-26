from Multi_samples_metrics import *
import seaborn as sns

path=["/Users/wangjiankang/figureServer/Nov2020/Rad21KD1_HiCmatrix/observed.KR.chr21.matrix.gz",
        "/Users/wangjiankang/figureServer/Nov2020/NIPBLKD1_HiCmatrix/observed.KR.chr21.matrix.gz",
        "/Users/wangjiankang/figureServer/Nov2020/Control1_HiCmatrix/observed.KR.chr21.matrix.gz",
        "/Users/wangjiankang/figureServer/Nov2020/Control2_HiCmatrix/observed.KR.chr21.matrix.gz"]
label=["Rad21KD","NIPBLKD","Control1","Control2"]

IS=getMultiSamplesScore(path,label,res=50000,chr="chr21",mode="IS",UniqueParameter=200000)
IS.corr()
sns.clustermap(IS.corr(),cmap="RdPu")
sns.clustermap(IS,row_cluster=False)
