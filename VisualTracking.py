import matplotlib as mpl
import matplotlib.pyplot as plt
from utils import *
import numpy as np
from numpy import genfromtxt
import datetime
from scipy import stats
import pandas as pd

from main import *
from main_stat import *
from main_plot import *
from main_create_plot_cvs import *
from main_stat_correlation import *
from main_pixelheatmap import *
from main_statCal import *
import time


start = time.time()

if __name__ == '__main__':
    print(datetime.datetime.now().strftime("%d.%b %Y %H:%M:%S"))
    
    fast_test=False
    
    strModality='Full' # Select this for 5 seconds analysis
    #strModality='FirstS' # Select this for 1st second analysis
    #strModality='Timeshift'
    nTimeInterval=1000
    nShiftInterval=100

    strEye='Right' # Select this for Right eye data
    #strEye='Left' # Select this for Left eye data
    #strVectorType='Vector'
    #strVectorType='DifferencesVector'
    #strVectorType='DistancesVector' # Select this for R/S and DFA Anlysis  
    strVectorType='2DPlot' # Select this for FD Anlysis 
    fLimitFDWindow=False
    bHurstIndex=False # Change TRUE for R/S Anlysis
    bDFA=False # for DFA both bHurstIndex and bDFA have to be true
    showPlot=True   #False True

    main(strModality=strModality, nTimeInterval=nTimeInterval, nShiftInterval=nShiftInterval,
              strEye=strEye,strVectorType=strVectorType,fLimitFDWindow=fLimitFDWindow,bHurstIndex=bHurstIndex,bDFA=bDFA,showPlot=showPlot,fast_test=fast_test)

    end = time.time()

