import pandas as pd
import seaborn as sns
import numpy as np
from scipy.stats import pearsonr
import matplotlib.pyplot as plt


e1 = pd.read_csv("../test_data/emissions_8.txt",sep="\t",index_col=0)
e2 = pd.read_csv("../test_data/emissions_9.txt",sep="\t",index_col=0)
e2.columns

e1=e1[e2.columns]
e1.index = ["e1_state"+str(i) for i in range(1,9)]
e2.index = ["e2_state"+str(i) for i in range(1,10)]
e1

corMT = pd.DataFrame(index=["e1_state"+str(i) for i in range(1,9)],
                    columns=["e2_state"+str(i) for i in range(1,10)])
for i in range(8):
    for j in range(9):
        corMT.iloc[i,j] = pearsonr(e1.iloc[i,:],e2.iloc[j,:])[0]
df = corMT.astype(float)

plt.figure(figsize=(9,8))
sns.heatmap(df,cmap="coolwarm")
