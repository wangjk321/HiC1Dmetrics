# HiC1Dmetrics
This repository contain code and basic tutorial for "HiC1Dmetrics"

![figure1](https://github.com/wangjk321/HiC1Dmetrics/blob/master/IMG/Figure1.png)

# Introduction

Basically, HiC1Dmetrics mainly provide three kinds of functions:

1. Calculate and visualize multiple 1-D metrics for one Hi-C samples.
    - Directional Index (DI) ([PMID: 22495300](https://pubmed.ncbi.nlm.nih.gov/22495300/))
    - Insulation Score (IS) ([PMID: 26030525](https://pubmed.ncbi.nlm.nih.gov/26030525/))
    - Contrast Index (CI) ([PMID: 24981874](https://pubmed.ncbi.nlm.nih.gov/24981874/))
    - TAD separation score (TADsep) ([PMID: 29335486](https://pubmed.ncbi.nlm.nih.gov/29335486/))
    - Distal-to-Local Ratio (DLR)  ([PMID: 30146161](https://pubmed.ncbi.nlm.nih.gov/30146161/))
    - More in building ...
2. Calculate and visulize 1-D metrics for comparing two Hi-C samples 
    - **Directional Ratio Frequency**, DRF (Original metric)
    - Differential DLR ([PMID: 30146161](https://pubmed.ncbi.nlm.nih.gov/30146161/))
    - Insulation score changes ([PMID: 31495782](https://pubmed.ncbi.nlm.nih.gov/31495782/))
    - More in building ...
3. Extract and visualize all "directional TAD" sites, which are defined by DRF metrics.

## Requirement

HiC1Dmetrics was based on `python 3.6`, and it require `numpy`, `pandas`, `scipy` and `matplotlib` library.

# Quick Start

The key function of "HiC1Dmetrics" is annotating and plotting "directional TAD". Here the example shows the right directional TAD on chromosome 21. The analysis is performed by comparing si-Rad21 with control RPE cells.

```python
from callDirectionalTAD import *
DirectionalTAD("../test_data/Rad21KD_1/observed.KR.chr21.matrix.gz","../test_data/Control_1/observed.KR.chr21.matrix.gz", \
               25000,chr="chr21").plotAlldirec("right")
```

![figure2](https://github.com/wangjk321/HiC1Dmetrics/blob/master/IMG/Figure2.png)

# Input file

We use zipped dense matrix as the input file. The format is like:

|       |  0   | 25000 | 50000 | 75000 | ...  |
| :---: | :--: | :---: | :---: | :---: | ---- |
|   0   |  0   |   0   |   0   |   0   | ...  |
| 25000 |  0   |   8   |   3   |   5   | ...  |
| 50000 |  0   |   3   |   8   |   4   | ...  |
| 75000 |  0   |   5   |   4   |   0   | ...  |
|  ...  | ...  |  ...  |  ...  |  ...  | ...  |

The easy way to generate input file is using JuicerTools:

```shell
java -jar juicer_tools.jar dump observed KR RPE.hic chr21 chr21 BP 25000 output.txt
```



# Usage

HiC1Dmetrics could be used in python and shell (will be added in version 0.2)

#### 1. Calculate 1-D metrics of ONE sample  (The output format is bedGraph): 

```python
from calculateMetrics import *
IS = InsulationScore("./Rad21KD_1/observed.KR.chr21.matrix.gz",25000,"chr21",\
                     out_name="InsulationScore",square_size=150000)
IS.getIS() #get IS
IS.getCSV() #export bedGraph file
```

|      |   chr |    start |      end | InsulationScore |
| :--- | ----: | -------: | -------: | --------------- |
| 0    | chr21 |        0 |    25000 | 0.0             |
| 1    | chr21 |    25000 |    50000 | 0.0             |
| 2    | chr21 |    50000 |    75000 | 0.0             |
| ...  |   ... |      ... |      ... | ...             |
| 1866 | chr21 | 46650000 | 46675000 | 0.0             |
| 1867 | chr21 | 46675000 | 46700000 | 0.0             |
| 1868 | chr21 | 46700000 | 46725000 | 0.0             |

`ContrastIndex`, `DirectionalityIndex`, `SeparationScore`, `DistalToLocal` are responsible for calculating CI,DI,TADsep,DLR, respectively.



#### 2. Plot 1-D metrics of one sample:

```python
from plotMetrics import *
plotm = PlotBedGraph("../test_data/Control_1/observed.KR.chr21.matrix.gz",25000,\
                     chr="chr21",startSite=1100*25000,endSite=1300*25000,\
                     title="One sample",clmin=0,clmax=50,other_parameter=0)
plotm.draw("IS") #plot scores
plotm.makePDF("IS","outputname") #export pdf
```

<img src="https://github.com/wangjk321/HiC1Dmetrics/blob/master/IMG/Figure3.png" width = "50%">

It is also possible to calculate 1-D metrics of multiple sample simultaneously:

```python
from Multi_samples_metrics import *
samples = ["../test_data/Rad21KD_1/observed.KR.chr20.matrix.gz",\
           "../test_data/Control_1/observed.KR.chr20.matrix.gz"]
labels = ["Rad21KD","Control"]
IS = getMultiSamplesScore(samples,labels,25000,"chr20","IS",150000)
CI = getMultiSamplesScore(samples,labels,25000,"chr20","CI",150000)
```



#### 3. Calculate 1-D metrics of two paired samples:

```python
from calculateMetrics import *
drf = DirectionalRelativeFreq("../test_data/Rad21KD_1/observed.KR.chr21.matrix.gz",\
                              "../test_data/Control_1/observed.KR.chr21.matrix.gz",\
                              25000,"chr21",out_name="DRF",
                start_distance=500000, end_distance=1000000)
drf.makeCSV() # export bedGraph file of DRF
```

#### 4. Plot 1-D metrics of two samples:

```python
from plotDiff import *
dplot = DiffDraw("../test_data/Rad21KD_1/observed.KR.chr21.matrix.gz",\
                 "../test_data/Control_1/observed.KR.chr21.matrix.gz",\
                 25000,startSite=1100*25000,endSite=1300*25000, clmin=-3,clmax=3,\
                 title="Comparison of two samples",chr="chr21",\
                 startDRF=500000,sizeDRF=1000000)
dplot.draw_square() # plot differential square matrix
dplot.draw_tri() #plot differential triangle matrix
dplot.draw_DRF() #plot DRF score
```

<img src="https://github.com/wangjk321/HiC1Dmetrics/blob/master/IMG/Figure4.png" width = "50%">

#### 5. Extract regions of "directional TAD" from differential contact matrix.

```python
from callDirectionalTAD import *
direcTAD = DirectionalTAD("../test_data/Rad21KD_1/observed.KR.chr21.matrix.gz",\
                          "../test_data/Control_1/observed.KR.chr21.matrix.gz",\
                          25000,chr="chr21",startDRF=500000,sizeDRF=1000000,\
                          sizeIS=150000)
leftTAD,rightTAD,_ = direcTAD.extractRegion() # The result is [leftTAD,rightTAD,DRFscore]
```

`rightTAD` is the region of "right directional TAD"

|  chr  | TADstart |  TADend  |
| :---: | :------: | :------: |
| chr21 | 15575000 | 15875000 |
| chr21 | 19000000 | 19450000 |
| chr21 | 22400000 | 22800000 |
| chr21 | 25375000 | 25575000 |
| chr21 | 28325000 | 28750000 |
| chr21 | 36100000 | 36375000 |

#### 6. Plot "directional TAD"

``` python
from callDirectionalTAD import *
plotDirecTAD = DirectionalTAD("../test_data/Rad21KD_1/observed.KR.chr21.matrix.gz",\
                              "../test_data/Control_1/observed.KR.chr21.matrix.gz",\
                              25000,chr="chr21",clmin=-2,clmax=2,title="Directional TAD on chr21",\
                              startDRF=500000,sizeDRF=1000000,sizeIS=150000)
plotDirecTAD.plotAlldirec("right") #plot rightTAD
plotDirecTAD.makePDF("right","rightTAD.pdf") #export to .pdf 
```

![figure6](https://github.com/wangjk321/HiC1Dmetrics/blob/master/IMG/Figure6.png)

#### 7. In addition. HiC1Dmetrics also provide simply function to visualize Hi-C data.

```python
from plotMetrics import *
#1 plot contact matrix (square)
PlotSquare("../test_data/Control_1/observed.KR.chr21.matrix.gz",25000,,chr="chr21",\
           startSite=1100*25000,endSite=1300*25000,title="Square matrix",\
           clmin=0,clmax=50).draw()

#2 plot contact matrix (triangle)
PlotTri("../test_data/Control_1/observed.KR.chr21.matrix.gz",25000,chr="chr21",\
        startSite=1100*25000,endSite=1300*25000,title="Tri matrix",\
        clmin=0,clmax=50).draw()

#3 plot TAD (TAD regions are calculated within HiC1Dmetrics)
PlotTAD("../test_data/Control_1/observed.KR.chr21.matrix.gz",25000,chr="chr21",\
        startSite=1100*25000,endSite=1300*25000,title="Plot TAD",\
        clmin=0,clmax=50).drawTAD()
```

![figure7](https://github.com/wangjk321/HiC1Dmetrics/blob/master/IMG/Figure7.png)

