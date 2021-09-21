import matplotlib.pyplot as plt
from numpy.polynomial.polynomial import polyfit
import math
import numpy as np
from hurst import compute_Hc
from dfa import dfa

def getFD(TimeSeriesCoords, fLimitFDWindow=False, showPlot=False):
    NumElem=len(TimeSeriesCoords[:,0])
    if (NumElem<=2):
        return -1 # suman 
    stepIncreseRatio=1.1
    
    #dinamic NumSteps based on series lenght
    NumSteps=math.log(NumElem)/math.log(stepIncreseRatio)
    #eventually limit NumSteps to avoid multifractal behaviour
    #(to apply in 2DPlot, Vector and DistanceVector, maybe not in the DifferenceVector)
    if (fLimitFDWindow==True):
        NumSteps/=2

    if (NumSteps<=1):
        return -1
    ArrayX=[]
    ArrayY=[]
    step=1.0
    FromValue=0
    ToValue=len(TimeSeriesCoords[:,0])
    for i in range(int(NumSteps)):
        numRules=GetSerieExactNumRulers(step,FromValue,ToValue,TimeSeriesCoords)
        if (numRules>0):    #CR: fix to the Nan issue
            ArrayX.append(step)
            ArrayY.append(numRules)
            step*=stepIncreseRatio

    ArrayX = np.asarray(ArrayX)
    ArrayX = 1/ ArrayX
    logX = np.log(ArrayX);
    logY = np.log(ArrayY);
    # Fit with polyfit
    b, m = polyfit(logX,logY, 1)
    ArrayY2 =  np.asarray(logX) * m + b
    if (showPlot):
        plt.plot(logX,logY,'ro')
        plt.title("FD = " + str(m))
        plt.plot(logX, ArrayY2, '-')
        plt.show()    
    return m


def __to_inc(x):
    incs = x[1:] - x[:-1]
    return incs

def __to_pct(x):
    pcts = x[1:] / x[:-1] - 1.
    return pcts


def getHurstIndex(TimeSeriesCoords, fLimitFDWindow, kind, showPlot=False, bDFA=False):
    ##commented for this particular case of Visual Tracking, see below*
    ##NumElem=len(TimeSeriesCoords[:,0])
    ##if (NumElem<=2):
    ##    return -1
    ##step=(TimeSeriesCoords[-1,0]-TimeSeriesCoords[0,0])/NumElem

    ##*in this particular case, I know the sampling rate is 4 ms, so I force the step to be 4 regardless af the points in the serie
    #step=4
    #NumElem=((TimeSeriesCoords[-1,0]-TimeSeriesCoords[0,0])/step)+1
    #series=[]
    #currVal=0
    #for i in range(int(NumElem)-1):
    #    #fill missing data
    #    #currVal=bisect_left(TimeSeriesCoords[1:,0], (i+1)*step)
    #    currVal=currVal+next(x[0] for x in enumerate(TimeSeriesCoords[currVal:,0]) if x[1] >= ((i+1)*step+TimeSeriesCoords[0,0]))
    #    currDist=((i+1)*step+TimeSeriesCoords[0,0])-TimeSeriesCoords[currVal-1,0]
    #    currY=((TimeSeriesCoords[currVal,1]-TimeSeriesCoords[currVal-1,1])/(TimeSeriesCoords[currVal,0]+TimeSeriesCoords[currVal-1,0]))*currDist+TimeSeriesCoords[currVal-1,1]
    #    series.append(currY)

    #series.append(TimeSeriesCoords[-1,1])
    #if (len(series)<100):
    #    #cumpute_HC does not work for series inferior to 100 elements
    #    return 0
    if (len(TimeSeriesCoords[:,1])<100) and (bDFA!=True):
        #cumpute_HC does not work for series inferior to 100 elements
        return 0
    # Evaluate Hurst equation
    try:
        if bDFA==True:
            NumElem=len(TimeSeriesCoords[:,0])
            NumSteps=int(math.log(NumElem,2))
            if kind == 'random_walk':
                x = __to_inc(TimeSeriesCoords[:,1])
            elif kind == 'price':
                x = __to_pct(TimeSeriesCoords[:,1])
            elif kind == 'change':
                x = TimeSeriesCoords[:,1]
            
            if fLimitFDWindow==True:       
                scales, fluct, alpha = dfa(x,scale_lim=[4,NumSteps], show=showPlot)
            else:
                scales, fluct, alpha = dfa(x,scale_lim=[2,NumSteps], show=showPlot)
            H = alpha
            data=scales
        else:
            if fLimitFDWindow==True:
                #H, c, data = compute_Hc(series, kind=kind, min_window=10, max_window=None, simplified=True)
                H, c, data = compute_Hc(TimeSeriesCoords[:,1], kind=kind, min_window=10, max_window=None, simplified=True)
            else:
                H, c, data = compute_Hc(TimeSeriesCoords[:,1], kind=kind, min_window=4, max_window=None, simplified=True)
    except:
        H = 0
        c = 0

    # Plot
    if showPlot==True and bDFA==False:
        f, ax = plt.subplots()
        ax.plot(data[0], c*data[0]**H, color="deepskyblue")
        ax.scatter(data[0], data[1], color="purple")
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel('Time interval')
        ax.set_ylabel('R/S ratio')
        ax.grid(True)
        plt.show()

    if bDFA==False:
        print("H={:.4f}, c={:.4f}".format(H,c))
    else:
        print("H={:.4f}".format(H))
    return H

   
#this is only approximative (works on time series and in 2DPlots, but is only an approximation)
def GetSerieLength(step,fromVal,toVal,TimeSeriesCoords):
    NumElem=toVal-fromVal
    if ((NumElem<=1) | (step<1)):
        return -1

    l=0.0;
    if (toVal>NumElem):
        toVal=NumElem
    for i in np.arange(fromVal+step,toVal,step):
        l+=math.sqrt(pow(TimeSeriesCoords[int(i),1]-TimeSeriesCoords[int(i-step),1],2)+pow((TimeSeriesCoords[int(i),0]-TimeSeriesCoords[int(i-step),0])/1000,2))
    return l

#this is only approximative (works quite good on time series, nut not very good in 2DPlots!)
def GetSerieNumRulers(step,fromVal,toVal,TimeSeriesCoords):
    NumElem=toVal-fromVal
    if ((NumElem<=1) | (step<1)):
        return -1

    l=0.0;
    if (toVal>NumElem):
        toVal=NumElem
    for i in np.arange(fromVal+step,toVal,step):
        l+=math.sqrt(pow(TimeSeriesCoords[int(i),1]-TimeSeriesCoords[int(i-step),1],2)+pow((TimeSeriesCoords[int(i),0]-TimeSeriesCoords[int(i-step),0])/1000,2))
    return l/step

#the right number of Rulers needed to cover the timeserie
def GetSerieExactNumRulers(step,fromVal,toVal,TimeSeriesCoords):
    NumElem=toVal-fromVal
    if ((NumElem<=1) | (step<1)):
        return -1

    N=0
    if (toVal>NumElem):
        toVal=NumElem
    currX=TimeSeriesCoords[fromVal,0] #suman paper of compass method
    currY=TimeSeriesCoords[fromVal,1]
    currVal=fromVal+1
    while currVal<toVal:
        #find the first farthest point from the circle centered in currX,currY with radius step
        bFound=False
        EPS=0.1 #use EPS to correct precision issues (if less that ESP, do not check circle with line intersection)
        for i in np.arange(currVal,toVal):
            currDist=math.sqrt(pow(TimeSeriesCoords[i,1]-currY,2)+pow(TimeSeriesCoords[i,0]-currX,2))# Suman
            if (abs(currDist-step)<=EPS):
                bFound=True
                N+=1
                currVal=i
                currX=TimeSeriesCoords[currVal,0]
                currY=TimeSeriesCoords[currVal,1]
            else:
                if (currDist>=step):
                    bFound=True
                    #optimization: if the next point is on the same segment of the current, I do not need to do the intersection
                    if (currVal==i):
                        N+=int(currDist/step)
                        currX+=(TimeSeriesCoords[currVal,0]-currX)*int(currDist/step)/(currDist/step)
                        currY+=(TimeSeriesCoords[currVal,1]-currY)*int(currDist/step)/(currDist/step)
                    else:
                        N+=1
                        currVal=i                    
                        #find intersection between circle and segment
                        
                        inp=LineIntersectCircle((currX,currY,step),(TimeSeriesCoords[currVal-1,0],TimeSeriesCoords[currVal-1,1]),(TimeSeriesCoords[currVal,0],TimeSeriesCoords[currVal,1]))
                        if len(inp)==1:
                            currX=inp[0][0]
                            currY=inp[0][1]
                        else:
                            if len(inp)==0:
                                #if this occurs, it means I got a resolution issue, so I get no intersection but I should take the new point as current
                                #N.B: I should never get there if EPS is big enought to avoid precision issue
                                currX=TimeSeriesCoords[currVal,0]
                                currY=TimeSeriesCoords[currVal,1]
                            else:
                                if ((pow(TimeSeriesCoords[currVal,1]-inp[0][1],2)+pow(TimeSeriesCoords[currVal,0]-inp[0][0],2))<((pow(TimeSeriesCoords[currVal,1]-inp[1][1],2)+pow(TimeSeriesCoords[currVal,0]-inp[1][0],2)))):
                                    currX=inp[0][0]
                                    currY=inp[0][1]
                                else:
                                    currX=inp[1][0]
                                    currY=inp[1][1]
                        
                    break
        if bFound==False:
            currVal=toVal
    return N

def LineIntersectCircle(p,lsp,lep):
    # p is the circle parameter, lsp and lep is the two end of the line
    x0,y0,r0 = p
    x1,y1 = lsp
    x2,y2 = lep
    if x1 == x2:
        if abs(r0) >= abs(x1 - x0):
            p1 = x1, y0 - math.sqrt(r0**2 - (x1-x0)**2)
            p2 = x1, y0 + math.sqrt(r0**2 - (x1-x0)**2)
            inp = [p1,p2]
            # select the points lie on the line segment
            inp = [p for p in inp if p[1]>=min(y1,y2) and p[1]<=max(y1,y2)]
        else:
            inp = []
    else:
        k = (y1 - y2)/(x1 - x2)
        b0 = y1 - k*x1
        a = k**2 + 1
        b = 2*k*(b0 - y0) - 2*x0
        c = (b0 - y0)**2 + x0**2 - r0**2
        delta = b**2 - 4*a*c
        if delta >= 0:
            p1x = (-b - math.sqrt(delta))/(2*a)
            p2x = (-b + math.sqrt(delta))/(2*a)
            p1y = k*p1x + b0
            p2y = k*p2x + b0
            inp = [[p1x,p1y],[p2x,p2y]]
            # select the points lie on the line segment
            inp = [p for p in inp if p[0]>=min(x1,x2) and p[0]<=max(x1,x2)]
        else:
            inp = []
    return inp

def EvaluateFDAndPlot(FDs, SerieX, SerieY, SerieTime, strVectorType, fLimitFDWindow, trial, StimulusName, showPlot, bHurstIndex, bDFA):
    #different modalities: 
    #one consider THE time, so add X and Y coords togheter, evaluate the variation and compare VS time 
    Vector=np.sqrt(np.square(SerieX)+np.square(SerieY))
    VectorSerie=np.c_[SerieTime, Vector]
    VectorDiff=[]
    VectorDistances=[]
    for i in range(len(Vector)-1):
        value = Vector[i+1] - Vector[i]
        VectorDiff.append(value)
        value = np.sqrt(np.square(SerieX[i+1]-SerieX[i])+np.square(SerieY[i+1]-SerieY[i]))
        VectorDistances.append(value)

    VectorSerieDiff=np.c_[SerieTime[1:], VectorDiff] 
    VectorSerieDistances=np.c_[SerieTime[1:], VectorDistances]

    if (strVectorType=='Vector'):
        if (showPlot):
            plt.plot(VectorSerie[:,0],VectorSerie[:,1])
            #plt.title("Vector")
            plt.xlabel('Time')
            plt.ylabel('Vector')
            plt.show()                
        if bHurstIndex==True:
            m = getHurstIndex(TimeSeriesCoords=VectorSerie, fLimitFDWindow=fLimitFDWindow, kind="random_walk", showPlot=showPlot, bDFA=bDFA)
        else:
            m = getFD(TimeSeriesCoords=VectorSerie, fLimitFDWindow=fLimitFDWindow, showPlot=showPlot)
        FDs.append((trial,StimulusName,m))
    if (strVectorType=='DifferencesVector'):
        if (showPlot):
            plt.plot(VectorSerieDiff[:,0],VectorSerieDiff[:,1])
            plt.title("Differences vector")
            plt.show()
        if bHurstIndex==True:
            m = getHurstIndex(TimeSeriesCoords=VectorSerieDiff, fLimitFDWindow=fLimitFDWindow, kind="change", showPlot=showPlot, bDFA=bDFA)
        else:
            m = getFD(TimeSeriesCoords=VectorSerieDiff, fLimitFDWindow=fLimitFDWindow, showPlot=showPlot)
        FDs.append((trial,StimulusName,m))

    if (strVectorType=='DistancesVector'):
        if (showPlot):
            plt.plot(VectorSerieDistances[:,0],VectorSerieDistances[:,1])
            #plt.title("Distances Vector")
            plt.xlabel('Time')
            plt.ylabel('Distance Vector')
            plt.show()
        if bHurstIndex==True:
            m = getHurstIndex(TimeSeriesCoords=VectorSerieDistances, fLimitFDWindow=fLimitFDWindow, kind="change", showPlot=showPlot, bDFA=bDFA)
        else:
            m = getFD(TimeSeriesCoords=VectorSerieDistances, fLimitFDWindow=fLimitFDWindow, showPlot=showPlot)
        FDs.append((trial,StimulusName,m))
                
    if (strVectorType=='2DPlot'):        
        #the second consider just the 2D plot of X and Y coords
        # if (showPlot):
        #     plt.plot(SerieX,SerieY)
        #     plt.xlim([0,1919])
        #     plt.ylim([1079,0])
        #     #plt.title("Brain-1")
        #     plt.title("FD = " + str(m))
        #     plt.xlabel('X')
        #     plt.ylabel('Y')
        #     plt.show()
        if bHurstIndex==True:
            m = getHurstIndex(TimeSeriesCoords=np.c_[SerieX,SerieY], fLimitFDWindow=fLimitFDWindow, kind="random_walk", showPlot=showPlot, bDFA=bDFA)
        else:
            m = getFD(TimeSeriesCoords=np.c_[SerieX,SerieY], fLimitFDWindow=fLimitFDWindow, showPlot=showPlot)
        FDs.append((trial,StimulusName,m))
        if (showPlot):
            plt.plot(SerieX,SerieY)
            plt.plot(SerieX,SerieY)
            plt.xlim([0,1919])
            plt.ylim([1079,0])
            #plt.title("Brain-1")
            plt.title("FD = " + str(np.abs(m)))
            plt.xlabel('X')
            plt.ylabel('Y')
            plt.show()

    if (showPlot):
        #plt.plot(SerieX[:,1].astype(np.float),SerieY[:,1].astype(np.float), color='C4', linestyle='--')
        plt.plot(SerieX,SerieY, color='r', linewidth=1)
        plt.xlim([0,1919])
        plt.ylim([1079,0])
        #plt.title(StimulusName)
        plt.title("Noise With Cross")
        im = plt.imread("Experiment Stimuli/" + StimulusName)
        implot = plt.imshow(im)
        plt.show()
