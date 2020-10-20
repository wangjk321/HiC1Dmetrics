import time
import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
#import altair as alt
from PIL import Image
#import datetime
from plotMetrics import *
from plotDiff import *
#import SessionState
import base64
from callDirectionalTAD import *
st.set_option('deprecation.showPyplotGlobalUse', False)
st.set_option('deprecation.showfileUploaderEncoding', False)

#Title
st.title("HiC1Dmetrics Web App (beta version)")
st.subheader("Here is a web-application for calculating, visualizing and download the results of 'HiC1Dmetrics'")
#image = Image.open('Figure1.png')
#st.image(image,use_column_width=True)

st.markdown('***')
st.markdown("# 🛠 For one sample")
st.markdown("### STEP1. Please choose your file from `juicer dump`; Or select test data from the sidebar.")

#Upload file
file = st.file_uploader("- Select your file:")
#or selcet example file
st.sidebar.markdown("# Test data:")
st.sidebar.markdown("## One sample 1D-metrics")
examplefile = st.sidebar.selectbox("- select the example data",["--Select--","Control.chr21.txt","Rad21KD.chr21.txt"])

if file is not None:
    with st.spinner("File loading..."):
        data = pd.read_csv(file,sep="\t",index_col=0)
    st.success("Successful loading from your data!")
elif examplefile is not "--Select--":
    with st.spinner("File loading..."):
        data = pd.read_csv(examplefile,sep="\t",index_col=0)
    st.success("Successful loading from test data!")

#Parameters
st.markdown("### STEP2. Set parameters")

resolution = st.number_input("- Resolution of your input", 25000,step = 5000)
chromosome = st.selectbox("- Chromosome:", ["chr1","chr2","chr3","chr4","chr5","chr6",\
                        "chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16",\
                        "chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"])
modeName = st.radio('- Choose the 1D metrics you wanted',('Insulation_Score', 'Contrast_Index',
                    'Directionality_Index','TAD_Separation_Score',"Distal_to_Local_Ratio"))
addtionPara = st.number_input("- Set the parameter for the selected metrics", 150000,step = 50000)
st.markdown("*(For insulation score, it is square_size; For contrast index, it is CI_size ...)*")

#Calculating
st.markdown("### STEP3. Run")
st.markdown("#### 👉 Obtain 1-D metrics")
whetherRun = st.button("Click to run",key="a")
bar = st.progress(0)
if whetherRun:
    if examplefile is "--Select--" and file is None:
        st.error("No input: Please choose your data,or select example data from sidebar.")
    else:
        if modeName == "Insulation_Score":
            metricClass = InsulationScore(data,resolution,chromosome,square_size=addtionPara)
        elif modeName == "Contrast_Index":
            metricClass = ContrastIndex(data,resolution,chromosome,CI_size=addtionPara)
        else:
            st.error("Under construction...")
        bar.progress(33)

        metric = metricClass.getScore()
        bar.progress(66)

        output = metric.to_csv(sep="\t", header=False, index=False)
        b64 = base64.b64encode(output.encode()).decode()
        outName = modeName+"_"+chromosome+".bedGraph"
        href = f'<a href="data:file/csv;base64,{b64}" download={outName}>Download bedGraph File</a> (right-click and save as &lt;some_name&gt;.bedGraph)'
        bar.progress(100)
        st.success("Done!")

        st.write("📎 Results of " + modeName+":",metric)
        st.write("📎 Download your results:")
        st.markdown(href, unsafe_allow_html=True)

#Visualization
st.markdown("#### 👉 Visualization")
whetherPlot = st.button("Click to plot selected region",key="b")

if examplefile is "--Select--" and file is None:
    if whetherPlot:
        st.error("No input: Please choose your data,or select example data from sidebar.")
else:
    plotRange = st.slider('Select the chromasom region:',0, resolution * data.shape[0],
                        (0, resolution*data.shape[0]),step=100000)
    if whetherPlot:
        if modeName == "Insulation_Score":
            type = "IS"
        elif modeName == "Contrast_Index":
            type = 'CI'

        plotClass = PlotBedGraph(data,resolution,chr=chromosome,other_parameter=addtionPara,
                                startSite=plotRange[0],endSite=plotRange[1])
        fig = plotClass.draw(type)
        st.pyplot(fig)

############################################################################################################
st.markdown('***')
# Two sample
st.markdown("# 🛠 For comparation of two sample")

@st.cache
def loadWithNorm(filename,method= "RPM",log = False):
    data = pd.read_csv(filename, delimiter='\t', index_col=0)
    if method == "RPM":
        data = (10000000 * data) / np.nansum(data)
    if log:
        return np.log1p(data)
    else:
        return data

#Upload file
st.markdown("### STEP1. Please upload (or select from the sidebar) two data for comparison (sample1 vs sample2)")
file1 = st.file_uploader("- Select your sample1 (treated):")
file2 = st.file_uploader("- Select your sample2 (control):")
#selce example file
st.sidebar.markdown("## Two sample 1D-metrics")
examplefile1 = st.sidebar.selectbox("- select the example data1 (treated)",["--Select--","Rad21KD.chr21.txt"])
examplefile2 = st.sidebar.selectbox("- select the example data2 (control)",["--Select--","Control.chr21.txt"])


if file1 is not None and file2 is not None:
    with st.spinner("File loading..."):
        data1 = loadWithNorm(file1,log = True)
        data2 = loadWithNorm(file2,log = True)
        dataForTreat= pd.read_csv(file1,sep="\t",index_col=0)
        dataForControl= pd.read_csv(file2,sep="\t",index_col=0)
    st.success("Successful loading from your data!")
elif examplefile1 is not "--Select--" and examplefile2 is not "--Select--":
    with st.spinner("File loading..."):
        data1 = loadWithNorm(examplefile1,log = True)
        data2 = loadWithNorm(examplefile2,log = True)
        dataForTreat= pd.read_csv(examplefile1,sep="\t",index_col=0)
        dataForControl= pd.read_csv(examplefile2,sep="\t",index_col=0)
    st.success("Successful loading from test data!")

#set parameter
st.markdown("### STEP2. Set parameters")
resolution = st.number_input("- Resolution of two sample", 25000,step = 5000)
chromosome = st.selectbox("- Chromosome of two sample:", ["chr1","chr2","chr3","chr4","chr5","chr6",\
                        "chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16",\
                        "chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"])
modeName = st.radio('- Choose the 1D metrics you wanted',('Directional_Ratio_of_Frequency', 'delta-DLR'))
DRFpara1 = st.number_input("- Parameter1 for DRF", min_value=0,value=500000,step = 50000)
DRFpara2 = st.number_input("- Parameter2 for DRF", min_value=0,value=1000000,step = 50000)

#Calculating
st.markdown("### STEP3. Run")
st.markdown("#### 👉 Obtain 1-D metrics for comparison of two sample")

whetherRunDRF = st.button("Click to run",key="drf")
bar = st.progress(0)
if whetherRunDRF:
    if (file1 is None or file2 is None) and (examplefile1 is "--Select--" or examplefile2 is "--Select--"):
        st.error("No input: Please choose your data,or select example data from sidebar.")
    else:
        if modeName == "Directional_Ratio_of_Frequency":
            metricClass = DirectionalRelativeFreq(data1,data2,resolution,chromosome,
                                        start_distance=DRFpara1, end_distance=DRFpara2)
        elif modeName == "delta-DLR":
            pass
        bar.progress(33)
        metric = metricClass.getDRF()
        bar.progress(66)

        output = metric.to_csv(sep="\t", header=False, index=False)
        b64 = base64.b64encode(output.encode()).decode()
        outName = modeName+"_"+chromosome+".bedGraph"
        href = f'<a href="data:file/csv;base64,{b64}" download={outName}>Download bedGraph File</a> (right-click and save as &lt;some_name&gt;.bedGraph)'
        bar.progress(100)
        st.success("Done!")

        st.write("📎 Results:",metric)
        st.write("📎 Download your results")
        st.markdown(href, unsafe_allow_html=True)

#Visualization
st.markdown("#### 👉Visualization")
whetherPlotDRF = st.button("Click to plot selected region",key="drf")
if  (file1 is None or file2 is None) and (examplefile1 is "--Select--" or examplefile2 is "--Select--"):
    if whetherPlotDRF:
        st.error("No input: Please choose your data,or select example data from sidebar.")
else:
    plotRangeDRF = st.slider('Select the chromosome range of plots: ()',0, resolution * data1.shape[0],
                    (0, resolution * data1.shape[0]),step=100000)

    if whetherPlotDRF:
        if modeName == "Directional_Ratio_of_Frequency":
            plotDRFClass = DiffDraw(data1,data2,resolution,chr=chromosome,startDRF=DRFpara1, sizeDRF=DRFpara2,
                                    startSite=plotRangeDRF[0],endSite=plotRangeDRF[1],clmin=-3,clmax=3)
            try:
                fig = plotDRFClass.draw_DRF()
                st.pyplot(fig)
            except:
                st.error("Please change the chromosome region")
        elif modeName == "delta-DLR":
            pass

#Call all directional TAD
st.markdown("#### 👉Extract directional TAD")
whetherExtractDRF = st.button("Click to extract directional region",key="extract")
bar = st.progress(0)
if whetherExtractDRF:
    if modeName == "Directional_Ratio_of_Frequency":
        extractClass = DirectionalTAD(data1,data2,25000,chr="chr21",startDRF=DRFpara1,
                                    sizeDRF=DRFpara2,sizeIS=150000)
    elif modeName == "delta-DLR":
        pass
    bar.progress(20)

    leftTAD,rightTAD,_ = extractClass.extractRegion(dataForControl)
    bar.progress(40)

    outleft = leftTAD.to_csv(sep="\t", header=False, index=False)
    b64left = base64.b64encode(outleft.encode()).decode()
    outNameleft = "leftTAD"+chromosome+".bedGraph"
    hrefleft = f'<a href="data:file/csv;base64,{b64left}" download={outNameleft}>Download bedGraph File</a> (right-click and save as &lt;some_name&gt;.bedGraph)'
    st.write("📎 Results of left-dTAD:",leftTAD)
    st.write("📎 Download the left-dTAD")
    st.markdown(hrefleft, unsafe_allow_html=True)
    bar.progress(60)

    outright = rightTAD.to_csv(sep="\t", header=False, index=False)
    b64right = base64.b64encode(outright.encode()).decode()
    outNameright = "rightTAD"+chromosome+".bedGraph"
    hrefright = f'<a href="data:file/csv;base64,{b64right}" download={outNameright}>Download bedGraph File</a> (right-click and save as &lt;some_name&gt;.bedGraph)'
    st.write("📎 Results of right-dTAD:",rightTAD)
    st.write("📎 Download the right-dTAD")
    st.markdown(hrefright, unsafe_allow_html=True)
    bar.progress(80)

    #plotleft
    st.write("📎 Show the left-dTAD")
    figleft = extractClass.plotAlldirec("left",dataForControl,dataForTreat)
    st.pyplot(figleft)
    bar.progress(90)

    st.write("📎 Show the right-dTAD")
    figright = extractClass.plotAlldirec("right",dataForControl,dataForTreat)
    st.pyplot(figright)
    bar.progress(100)
    st.success("Finish!")


#防止点击按钮后刷新，之后添加此功能
#def main():
#    st.subheader("new")
#
#    session_state = SessionState.get(name="", button_sent=False)
#
#    session_state.name = st.text_input("Enter your name")
#    button_sent = st.button("Send")
#
#    if button_sent:
#        session_state.button_sent = True
#
#    if session_state.button_sent:
#        st.write(session_state.name)
#
#        session_state.bye = st.checkbox("bye")
#        session_state.welcome = st.checkbox("welcome")
#
#        if session_state.bye:
#            st.write("I see")
#        if session_state.welcome:
#            st.write("you see")

#main()
