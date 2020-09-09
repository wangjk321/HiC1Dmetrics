# HiC1Dmetrics
This repository contain code and basic tutorial for "HiC1Dmetrics"

![figure1](https://github.com/wangjk321/HiC1Dmetrics/blob/master/IMG/Figure1.png)

# Introduction

Basically, HiC1Dmetrics mainly provide three types of function:

1. Calculate and visulize multiple 1-D metrics for one Hi-C samples (Already published idea).
    - Directional Index (DI)
    - Insulation Score (IS)
    - Contrast Index (CI)
    - TAD separation score (TADsep)
    - Distal-to-Local Ratio (DLR)
    - More in building..

2. Calculate and visulize 1-D metrics for comparing two Hi-C samples 
    - **Directional Frequency Ratio**, DFR (Original metric)
    - Differential DLR
    - More in building ..

3. Extract and visulzie all "directional TAD" sites, which are defined by DFR metrics.

# Quick Start

```python
DirectionalTAD("../test_data/Rad21KD_1/observed.KR.chr21.matrix.gz","../test_data/Control_1/observed.KR.chr21.matrix.gz",25000,chr="chr21").plotAlldirec("right")
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

The esaist way to generate input file is using JuicerTools:

```shell
java -jar juicer_tools.jar dump observed KR RPE.hic chr21 chr21 BP 25000 output.txt
```



# Usage
1. Calculate 1-D metrics of ONE sample  (The output format is bedGraph): 

```python
from calculateMetrics import *
IS = InsulationScore("./Rad21KD_1/observed.KR.chr21.matrix.gz",25000,"chr21",out_name="InsulationScore",square_size=150000)
IS.getCSV()
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

`ContrastIndex`, `DirectionalityIndex`, `SeparationScore`, `DistalToLocal` is responsible for calculate CI,DI,TADsep,DLR, respectively.



2. Plot 1-D metrics of one sample:

```python
plotm = PlotBedGraph("../test_data/Control_1/observed.KR.chr21.matrix.gz",25000,startSite=1100*25000,endSite=1300*25000,title="One sample",chr="chr21",clmin=0,clmax=50,other_parameter=0)
plotm.draw("IS")
plotm.makePDF("IS","outputname")
```

![figure3](https://github.com/wangjk321/HiC1Dmetrics/blob/master/IMG/Figure3.png)

3. Calculate 1-D metrics of two samples:

```python
drf = DirectionalRelativeFreq("../test_data/Rad21KD_1/observed.KR.chr21.matrix.gz","../test_data/Control_1/observed.KR.chr21.matrix.gz",25000,"chr21",out_name="DRF",
                start_distance=500000, end_distance=1000000)
drf.makeCSV()
```

4. Plot 1-D metrics of two samples:

```python
dplot = DiffDraw("../test_data/Rad21KD_1/observed.KR.chr21.matrix.gz","../test_data/Control_1/observed.KR.chr21.matrix.gz",25000,startSite=1100*25000,endSite=1300*25000, clmin=-3,clmax=3,title="Comparison of two samples",chr="chr21",startDRF=500000,sizeDRF=1000000)
dplot.draw_square()
dplot.draw_tri()
dplot.draw_DRF()
```

![figure4](https://github.com/wangjk321/HiC1Dmetrics/blob/master/IMG/Figure4.png){:height="50%" width="50%"}

5. Extract regions of "directional TAD" from differential contact matrix.

```python
direcTAD = DirectionalTAD("../test_data/Rad21KD_1/observed.KR.chr21.matrix.gz","../test_data/Control_1/observed.KR.chr21.matrix.gz",25000,chr="chr21",startDRF=500000,sizeDRF=1000000,sizeIS=150000)
leftTAD,rightTAD,_ = direcTAD.extractRegion()
rightTAD
```

|  chr  | TADstart |  TADend  |
| :---: | :------: | :------: |
| chr21 | 15575000 | 15875000 |
| chr21 | 19000000 | 19450000 |
| chr21 | 22400000 | 22800000 |
| chr21 | 25375000 | 25575000 |
| chr21 | 28325000 | 28750000 |
| chr21 | 36100000 | 36375000 |

6. Plot "directional TAD"

``` python
plotDirecTAD = DirectionalTAD("../test_data/Rad21KD_1/observed.KR.chr21.matrix.gz","../test_data/Control_1/observed.KR.chr21.matrix.gz",25000,chr="chr21",clmin=-2,clmax=2,title="Directional TAD on chr21",startDRF=500000,sizeDRF=1000000,sizeIS=150000)
plotDirecTAD.plotAlldirec("right")
plotDirecTAD.makePDF("right","rightTAD.pdf")
```

![figure6](https://github.com/wangjk321/HiC1Dmetrics/blob/master/IMG/Figure6.png)

7. Others. HiC1Dmetrics all provide function to simply visualize Hi-C data.

```python
#1 plot contact matrix (square)
PlotSquare("../test_data/Control_1/observed.KR.chr21.matrix.gz",25000,startSite=1100*25000,endSite=1300*25000,title="Square matrix",chr="chr21",clmin=0,clmax=50).draw()
#2 plot contact matrix (triangle)
PlotTri("../test_data/Control_1/observed.KR.chr21.matrix.gz",25000,startSite=1100*25000,endSite=1300*25000,title="Tri matrix",chr="chr21",clmin=0,clmax=50).draw()
#3 plot TAD (TAD regions are calculated within HiC1Dmetrics)
PlotTAD("../test_data/Control_1/observed.KR.chr21.matrix.gz",25000,startSite=1100*25000,endSite=1300*25000,title="Plot TAD",chr="chr21",clmin=0,clmax=50).drawTAD()
```

![figure7](https://github.com/wangjk321/HiC1Dmetrics/blob/master/IMG/Figure7.png)

