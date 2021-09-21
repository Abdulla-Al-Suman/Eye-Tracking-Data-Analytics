import matplotlib as mpl
import matplotlib.pyplot as plt
from utils import *
import numpy as np
from numpy import genfromtxt
import datetime
from scipy import stats
import pandas as pd

def main(strModality, nTimeInterval, nShiftInterval, 
                strEye,strVectorType,fLimitFDWindow=False,bHurstIndex=False,bDFA=False,showPlot=False,fast_test=False):

    fileList = []
    with open('filelist_Deidentified_Extracted.txt') as f
        for line in f:
            fileList.append(line[:-1])

    print('Total files: ' + str(len(fileList)))
    if bDFA==True:
        bHurstIndex=True

    allFDs = []
    bFirstFD=True
    i=0
    for file in fileList:
        print(datetime.datetime.now().strftime("%d.%b %Y %H:%M:%S"))
        print('Reading File: ' + file)
        i+=1
        if (fast_test==True):
            vt_data = genfromtxt('test' + str(i) + '.txt', delimiter=',', skip_header=4, max_rows=8000, dtype=np.object)
        else:
            vt_data = genfromtxt("Deidentified _Extracted _Data/" + file, delimiter='\t', skip_header=0, dtype=np.object)
        vt_Trials=np.unique(vt_data[:,2])
        FDs = []
        FDs.append(('','',file[:-4]))
        for trial in vt_Trials[1:]:
            condition=vt_data[:,2]==trial
            condition&=vt_data[:,9]==b'Information'
            title_data=vt_data[condition,:]
            
            #Start to analyze just the Right Eye, so in this moment I comment the left eye information conditions
            condition=vt_data[:,2]==trial
            condition&=vt_data[:,9]!=b'Information'
            #condition&=((vt_data[:,11]==b'Fixation') | (vt_data[:,11]==b'Saccade'))     #Category Right
            ##condition&=((vt_data[:,12]==b'Fixation') | (vt_data[:,12]==b'Saccade'))    #Category Left
            #condition&=vt_data[:,15]!=b'0.0' #Pupil Diameter Right
            ##condition&=vt_data[:,16]!=b'0.0' #Pupil Diameter Left
            #condition&=vt_data[:,33]!=b'0.0' #Pupil Pos Right X
            #condition&=vt_data[:,34]!=b'0.0' #Pupil Pos Right Y
            ##condition&=vt_data[:,35]!=b'0.0' #Pupil Pos Left X
            ##condition&=vt_data[:,36]!=b'0.0' #Pupil Pos Left Y
            filtered_data=vt_data[condition,:]
            TimeData=filtered_data[:,[0]].astype(np.float)
            StimulusData=title_data[:,[3]]
            StimulusName=StimulusData[0][0].decode('utf-8')
            print(datetime.datetime.now().strftime("%d.%b %Y %H:%M:%S"))
            print('File: ' + file + ' - Current Stimulus:' + StimulusName)
            if (len(filtered_data)==0):
                print('No useful data for ' + StimulusName)
                m=0
                FDs.append((trial,StimulusName,m))
            else:                
                #shifting FD values
                if (strModality=='Timeshift'):
                    #split TimeData in nTimeInterval millisec intervals, shifting nShiftInterval ms
                    print('Time init: ' + str(TimeData[0,0]) + '; Time final: ' + str(TimeData[-1,0]))
                    StartingTime=TimeData[0,0]
                    EndingTime=StartingTime+nTimeInterval
                    if EndingTime>TimeData[-1,0]:
                        EndingTime=TimeData[-1,0]
                    while EndingTime<=TimeData[-1,0]:
                        print('Staring time: ' + str(StartingTime) + '; Ending time: ' + str(EndingTime))
                        StartingTime=StartingTime+nShiftInterval
                        EndingTime=StartingTime+nTimeInterval
                        if (EndingTime>TimeData[-1,0]) and (TimeData[-1,0]-StartingTime>(nTimeInterval-nShiftInterval)):
                            EndingTime=TimeData[-1,0]
                        
                        conditionTime=filtered_data[:,0].astype(np.float)>=StartingTime
                        conditionTime&=filtered_data[:,0].astype(np.float)<=EndingTime
                        conditionTime&=((filtered_data[:,11]==b'Fixation') | (filtered_data[:,11]==b'Saccade'))     #Category Right
                        conditionTime&=filtered_data[:,13]!=b'0.0' #Pupil Diameter Right
                        conditionTime&=filtered_data[:,22]!=b'0.0' #Pupil Pos Right X
                        conditionTime&=filtered_data[:,23]!=b'0.0' #Pupil Pos Right Y
                        SerieTime=filtered_data[conditionTime,[0]].astype(np.float)
                        PointofRegardRightX=filtered_data[conditionTime,[14]].astype(np.float)
                        PointofRegardRightY=filtered_data[conditionTime,[15]].astype(np.float)
                        # PointofRegardLeftX=filtered_data[conditionTime,[19]].astype(np.float)
                        # PointofRegardLeftY=filtered_data[conditionTime,[20]].astype(np.float)
                        EvaluateFDAndPlot(FDs=FDs, SerieX=PointofRegardRightX, SerieY=PointofRegardRightY, SerieTime=SerieTime, strVectorType=strVectorType, fLimitFDWindow=fLimitFDWindow, trial=trial, StimulusName=StimulusName, showPlot=showPlot, bHurstIndex=bHurstIndex, bDFA=bDFA)
                else:
                    # # for selecting 1st seconds
                    #conditionTime=((filtered_data[:,10]==b'Fixation') | (filtered_data[:,10]==b'Saccade'))
                    # StartingTime=TimeData[0,0]
                    # EndingTime=StartingTime+1000 
                    # if EndingTime>TimeData[-1,0]:
                    #     EndingTime=TimeData[-1,0]
                    # conditionTime=filtered_data[:,0].astype(np.float)>=StartingTime
                    # conditionTime&=filtered_data[:,0].astype(np.float)<=EndingTime 
                    # #conditionTime&=(filtered_data[:,10]==b'Fixation') # only Fixation points 
                    # conditionTime&=((filtered_data[:,10]==b'Fixation') | (filtered_data[:,10]==b'Saccade'))

                    # for selecting 5 seconds
                    conditionTime=((filtered_data[:,10]==b'Fixation') | (filtered_data[:,10]==b'Saccade')) # for 1st second, this line need deactivate and activate upper lines
                    conditionTime&=filtered_data[:,12]!=b'0.0' #Pupil Diameter Right
                    conditionTime&=filtered_data[:,21]!=b'0.0' #Pupil Pos Right X 
                    conditionTime&=filtered_data[:,22]!=b'0.0' #Pupil Pos Right Y
                    SerieTime=filtered_data[conditionTime,[0]].astype(np.float)
                    PointofRegardRightX=filtered_data[conditionTime,[13]].astype(np.float)
                    PointofRegardRightY=filtered_data[conditionTime,[14]].astype(np.float)
   
                    EvaluateFDAndPlot(FDs=FDs, SerieX=PointofRegardRightX, SerieY=PointofRegardRightY, SerieTime=SerieTime, strVectorType=strVectorType, fLimitFDWindow=fLimitFDWindow, trial=trial, StimulusName=StimulusName, showPlot=showPlot, bHurstIndex=bHurstIndex, bDFA=bDFA) # suman what returns


        FDs = np.asarray(FDs)
        #np.savetxt('FDs_' + strEye + '_' + strVectorType + '_' + strModality + '_' + file, FDs, delimiter=";", fmt="%s")
        if (bFirstFD==True):
            allFDs=FDs
            bFirstFD=False
        else:
            fd_Trials=np.unique(allFDs[:,0])
            bFirstRow=True
            for trial in fd_Trials:
                cond1=allFDs[:,0]==trial
                cond2=FDs[:,0]==trial
                numData=max(len(allFDs[cond1]),len(FDs[cond2]))
                if (bFirstRow==True):
                    allFDPadded = np.asarray(allFDs[cond1,:])
                else:
                    allFDPadded=np.concatenate((allFDPadded,allFDs[cond1,:]))
                if ((numData-len(allFDs[cond1,:]))>0):
                    tmpRow=np.c_[[[trial,allFDs[cond1,1][0]]],[np.zeros(len(allFDPadded[0,:])-2)]]
                    allFDPadded=np.concatenate((allFDPadded,np.tile(tmpRow,((numData-len(allFDs[cond1,:])),1))))
                
                if (bFirstRow==True):
                    FDPadded=np.asarray(FDs[cond2,:])
                else:
                    FDPadded=np.concatenate((FDPadded,FDs[cond2,:]))

                bFirstRow=False
                if ((numData-len(FDs[cond2,:]))>0):
                    FDPadded=np.concatenate((FDPadded,[[trial,FDs[cond2,1][0],0]]*(numData-len(FDs[cond2,:]))))

            allFDs=np.c_[allFDPadded,FDPadded[:,2]]

    if fLimitFDWindow==True:
        strLimit='_LimFDWin'
    else:
        strLimit=''
    if (strModality=='Timeshift'):
        strShiftParams='_' + str(nTimeInterval) + '_' + str(nShiftInterval)
    else:
        strShiftParams=''
    print(datetime.datetime.now().strftime("%d.%b %Y %H:%M:%S"))
    if bHurstIndex==True:
        if bDFA==True:
            np.savetxt('allHI_DFAs_BothFixationsSaccades63participants_' + strEye + '_' + strVectorType + '_' + strModality + strShiftParams + strLimit + '.csv', allFDs, delimiter=",", fmt="%s")
        else:
            np.savetxt('allHIs_BothFixationsSaccades63participants_' + strEye + '_' + strVectorType + '_' + strModality + strShiftParams + strLimit + '.csv', allFDs, delimiter=",", fmt="%s")
    else:
        np.savetxt('allFDs_BothFixationsSaccades63participants_' + strEye + '_' + strVectorType + '_' + strModality + strShiftParams + strLimit + '.csv', allFDs, delimiter=",", fmt="%s")
    print(datetime.datetime.now().strftime("%d.%b %Y %H:%M:%S"))

