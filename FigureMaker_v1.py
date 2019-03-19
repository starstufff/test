#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 22 15:30:21 2019

@author: oleg
"""


import sys
sys.path.append('/Users/StarStuff/Dropbox/pythonfiles/Python files 280119/')

import pandas as pd
import numpy as np
import matplotlib as mpl
mpl.use('GTKAgg') # to use GTK UI
import matplotlib.pyplot as plt
import random
from mpl_toolkits.mplot3d import Axes3D
from sklearn.decomposition import PCA 
from scipy.special import binom 
from scipy.special import erf 


from CellTrackModule import CellTrackClass
from CellTrackModule import CellTrackArrayClass
from CellTrackModule import DoContourWritheCalc
from CellTrackModule import CellTrackMatchClass
from CellTrackModule import ContourLength
#import importlib
#import FitToolModule #import the module here, so that it can be reloaded.
#from FitToolModule import FitTool
#import FitToolModule.FitTool

#%%
ExcelFilePath = '/Users/StarStuff/Dropbox/pythonfiles/intravital measurements_013019.xlsx'
trks = CellTrackArrayClass()
trks.DoLoadDataFromExcelFile(ExcelFilePath)
trks.DoSortTracksByParam('dt_s') 
HSC = trks.DoGetCellTrackSubarrayByParamValue({'cellType' : 'hsc'})
MPHG = trks.DoGetCellTrackSubarrayByParamValue({'cellType' : 'macrophage'})
#%% exclude post Rx data
print(len(HSC.tracks))
HSCnodrug= HSC.DoExcludeCellTracksBySubstring({'metaname' : 'PostRx'})
HSCnodrug= HSCnodrug.DoExcludeCellTracksBySubstring({'metaname' : 'postRx'})
HSCnodrug= HSCnodrug.DoExcludeCellTracksBySubstring({'metaname' : 'Post'})
print(len(HSCnodrug.tracks))
#%% exclude pre Rx data
print(len(HSC.tracks))
HSCdrug= HSC.DoExcludeCellTracksBySubstring({'metaname' : 'PreRx'})
HSCdrug= HSCdrug.DoExcludeCellTracksBySubstring({'metaname' : 'preRx'})
HSCdrug= HSCdrug.DoExcludeCellTracksBySubstring({'metaname' : 'Pre'})
print(len(HSCdrug.tracks))
#%%
HSCnodrug.DoCalculateMSD()
MPHG.DoCalculateMSD()
HSCdrug.DoCalculateMSD()
#%% Fig. 1 track examples
plt.cla()
MPHG5801 =  MPHG.DoGetCellTrackSubarrayByParamValue({'animalID': 5801})
MPHG5801.DoPlotXYtrack(startXYloc = [0, 0], markStartEnd = False, picker = 5, label = None, color = "lightgrey")
#%%
HSC5272 =  HSCnodrug.DoGetCellTrackSubarrayByParamValue({'animalID': 5272}) 
HSC5272.DoPlotXYtrack(label = None, startXYloc = [0, 0], markStartEnd = False)

#%%
HSC3509 =  HSCnodrug.DoGetCellTrackSubarrayByParamValue({'animalID': 3509})
HSC3509.DoPlotXYtrack(label = None, startXYloc = [0, 0], markStartEnd = False)

#%% Fig. 1 track - combine tracks by trk number
for trk in HSCnodrug.tracks[:]:
    #if trk.maxDist >= 20 : #filter tracks with maxdist >= 20um
        trk.DoPlotXYtrack(label = None, startXYloc = [0, 0], markStartEnd = False, picker = 5, alpha = 0.5)

#%%
plt.cla()
for trk in MPHG.tracks[30:40]:
    trk.DoPlotXYtrack(label = None, startXYloc = [0, 0], markStartEnd = False, picker = 5)
    
#%% Fig. 1 population - track correlation
#% load population data
D = pd.read_excel('/Users/StarStuff/Dropbox/pythonfiles/tomato+ cell counts_intravital biopsiesNumeric names.xlsx')

#%% first check maximal displacement histograms
plt.cla()
maxDist = [trk.maxDist for trk in HSCnodrug.tracks]
maxDistMPHG = [trk.maxDist for trk in MPHG.tracks]
H, bins = np.histogram(maxDist, bins = np.arange(0, 100, 10))
fig, ax1 = plt.subplots()
plt.cla()


plt.plot((bins[:-1]+bins[1:])/2, H, '-o', label = 'Tomato+', color = 'tab:red')
plt.xlabel('Maximal distance on track ($\mu m$)')
plt.ylabel('Tomato+ cell No')
plt.ylim((0, 28))
#plt.legend()
ax2 = plt.gca().twinx()  # instantiate a second axes that shares the same x-axis

Hmphg, bins = np.histogram(maxDistMPHG, bins = bins)
plt.plot((bins[:-1]+bins[1:])/2, Hmphg, '-o', label = 'macrophages')
plt.xlabel('Maximal distance on track ($\mu m$)')
plt.ylabel('Macrophage cell No')
YLim = np.array([0, 470])
plt.ylim(YLim)
#plt.legend()

thresh = 20
plt.plot([thresh, thresh], YLim, '--', color = 'k')
plt.show()

#%% Fig. 1 population vs long track correlation (exclude calvaria tracks)
HSCtibia = HSCnodrug.DoExcludeCellTracksBySubstring({'CompleteMetaName' : '98'}) 
HSCtibia = HSCtibia.DoExcludeCellTracksBySubstring({'CompleteMetaName' : '1303'})
HSCtibia = HSCtibia.DoExcludeCellTracksBySubstring({'CompleteMetaName' : '3509'})


CellArray = HSCtibia #HSCnodrug #HSCtibia

maxDist_thresh = 20
meanVar = []  
thresh = 0.88
fracAboveThresh = [] 
HSC = [] 
STHSC = []
MPP2 = []
MPP34 = []
MyP = []
MEP = []
ScaPcKitN = []
ScaNcKitN = []
Lin = []
AnimalID = []
meanRg = []
NoTracks = []
meanMaxDist = []

expl_var = np.array([trk.MSD.PCA_explained_var for trk in CellArray.tracks])
animalID = np.array([trk.animalID for trk in CellArray.tracks])
Rg = np.array([trk.Rg_um for trk in CellArray.tracks])
maxDist = np.array([trk.maxDist for trk in CellArray.tracks])

animalIDunique= set(animalID)

varByID = dict()
RgByID = dict()
maxDistByID = dict()
for aID in animalIDunique:
    varByID[aID] = expl_var[animalID == aID]   
    RgByID[aID] = Rg[animalID == aID]
    maxDistByID[aID] = maxDist[animalID == aID]
    
for aID in varByID:
    NoTracks.append(varByID[aID].size)
    AnimalID.append(aID)
    meanVar.append(np.mean(varByID[aID]))
    meanRg.append(np.mean(RgByID[aID]))
    meanMaxDist.append(np.mean(maxDistByID[aID]))
    fracAboveThresh.append(np.sum(maxDistByID[aID] > maxDist_thresh)/maxDistByID[aID].size)
    HSC.append(D.loc[D['Animal ID'] == aID]['HSC'].values[0])
    STHSC.append(D.loc[D['Animal ID'] == aID]['ST-HSC'].values[0])
    MPP2.append(D.loc[D['Animal ID'] == aID]['MPP2'].values[0])
    MPP34.append(D.loc[D['Animal ID'] == aID]['MPP3/4'].values[0])
    MyP.append(D.loc[D['Animal ID'] == aID]['MyP'].values[0])
    MEP.append(D.loc[D['Animal ID'] == aID]['MEP'].values[0])
    ScaPcKitN.append(D.loc[D['Animal ID'] == aID]['Sca+ cKit-'].values[0])
    ScaNcKitN.append(D.loc[D['Animal ID'] == aID]['Sca- cKit-'].values[0])
    Lin.append(D.loc[D['Animal ID'] == aID]['Lin+'].values[0])

NoTracks =np.array(NoTracks)
meanVar = np.array(meanVar)  
fracAboveThresh = np.array(fracAboveThresh) 
HSC = np.array(HSC)
STHSC = np.array(STHSC)
MPP2 = np.array(MPP2)
MPP34 = np.array(MPP34)
MyP = np.array(MyP)
MEP = np.array(MEP)
ScaPcKitN = np.array(ScaPcKitN)
ScaNcKitN = np.array(ScaNcKitN)
Lin = np.array(Lin)
AnimalID = np.array(AnimalID)
meanMaxDist = np.array(meanMaxDist)

AllPops = HSC + STHSC + MPP2 + MPP34 + MyP + MEP + ScaPcKitN + ScaNcKitN + Lin
HSCfrac = HSC/AllPops

HSCfracTrack =  np.array([HSCfrac[AnimalID == trk.animalID][0] 
                        for trk in CellArray.tracks])


fig, ax1 = plt.subplots()
print('\n HSC correlation = \n', np.corrcoef(HSCfrac, fracAboveThresh))
plt.cla()
plt.plot(HSCfrac, fracAboveThresh,  'o', color = 'r', picker = 5)
plt.xlabel('HSC fraction')
plt.ylabel('Fraction of tracks with max dist > 20um')
#plt.title('All data')
plt.show()


STHSCfrac = STHSC/AllPops
print('\n ST-HSC correlation = \n', np.corrcoef(STHSCfrac, fracAboveThresh))
MPP2frac = MPP2/AllPops
print(np.corrcoef(MPP2frac, fracAboveThresh))
#MPP34frac = MPP34/AllPops
#print(np.corrcoef(MPP34frac, meanRg))
MyPfrac = MyP/AllPops
print(np.corrcoef(MyPfrac, fracAboveThresh))
MEPfrac = MEP/AllPops
print(np.corrcoef(MEPfrac, fracAboveThresh))
ScaPcKitNfrac = ScaPcKitN/AllPops
print(np.corrcoef(ScaPcKitNfrac, fracAboveThresh))
ScaNcKitNfrac = ScaNcKitN/AllPops
print(np.corrcoef(ScaNcKitNfrac, fracAboveThresh))

#%% mean track velocities plot
plt.cla()
CellArray = HSCnodrug
animalIDs = np.array([trk.animalID for trk in CellArray.tracks])
Vmean = np.array([trk.Vmean_um_s for trk in CellArray.tracks])
ticks = []
labels = []
for idx, aID in enumerate([4113,2747,4107,3809,3111,3114,3215,5009,1296]): #enumerate(set(animalIDs)):
    Vanml = Vmean[animalIDs == aID]*60
    Vmphg = np.array([trk.Vmean_um_s for trk in MPHG.tracks if trk.animalID == aID])
    plt.plot(idx*np.ones(Vanml.size), Vanml, 'o', color = 'red')
    plt.plot(idx*np.ones(Vmphg.size), Vmphg, 'o', color = 'blue')
    ticks.append(idx)
    labels.append(str(aID))

plt.xticks(ticks, labels) #print animalID labels
#plt.xticks(ticks) #prints random number
plt.ylabel('Mean track velocity ($\mu$m/min)')
plt.show()
#print list of animal IDs
#print(set(animalIDs))

#%% XY only: mean track velocities plot
plt.cla()
CellArray = HSCnodrug
animalIDs = np.array([trk.animalID for trk in CellArray.tracks])
Vxymean = np.array([trk.Vxymean_um_s for trk in CellArray.tracks])
ticks = []
labels = []
for idx, aID in enumerate([4113,2747,4107,3809,3111,3114,3215,5009,1296]): #enumerate(set(animalIDs)):
    Vanml = Vxymean[animalIDs == aID]*60
    Vmphg = np.array([trk.Vxymean_um_s for trk in MPHG.tracks if trk.animalID == aID])
    plt.plot(idx*np.ones(Vanml.size), Vanml, 'o', color = 'red')
    plt.plot(idx*np.ones(Vmphg.size), Vmphg, 'o', color = 'blue')
    ticks.append(idx)
    labels.append(str(aID))
    
for idx, aID in enumerate([4113,2747,4107,3809,3111,3114,3215,5009,1296]): #enumerate(set(animalIDs)):
    Vanml = Vxymean[animalIDs == aID]*60
    Vmphg = np.array([trk.Vxymean_um_s for trk in MPHG.tracks if trk.animalID == aID])
    
    print (aID)
    print ("ID above, red first, blue second")
    print (Vanml)
    print ("blue")
    print (Vmphg)
    

plt.xticks(ticks, labels) #print animalID labels
#plt.xticks(ticks) #prints random number
plt.ylabel('Mean track velocity XY ($\mu$m/min)')
plt.show()
#print list of animal IDs
#print(set(animalIDs))
#%%
plt.cla()
CellArray = HSCnodrug
animalIDs = np.array([trk.animalID for trk in CellArray.tracks])
Dist = np.array([trk.distance_um for trk in CellArray.tracks])
ticks = []
labels = []
for idx, aID in enumerate([2747,4113,4107,3114,3111,3809,5008,633,1296,5009,5272,2454,2545]): #enumerate(set(animalIDs)):
    DistHSC = Dist[animalIDs == aID]
    Distmphg = np.array([trk.distance_um for trk in MPHG.tracks if trk.animalID == aID])
    plt.plot(idx*np.ones(DistHSC.size), DistHSC, 'o', color = 'red')
    plt.plot(idx*np.ones(Distmphg.size), Distmphg, 'o', color = 'blue')
    ticks.append(idx)
    labels.append(str(aID))

plt.xticks(ticks, labels) #print animalID labels
#plt.xticks(ticks) #prints random number
plt.ylabel('Track Length ($\mu$m/)')
plt.show()
#print list of animal IDs
#print(set(animalIDs))

#%% XY plots before and after Rx
plt.cla()
HSC2454 =  HSCnodrug.DoGetCellTrackSubarrayByParamValue({'animalID': 2454})
HSC2454drug = HSC.DoGetCellTrackSubarrayBySubstring({'metaname' : 'post'}).DoGetCellTrackSubarrayByParamValue({'animalID': 2454})

HSC2454.DoPlotXYtrack(label = None, color = 'r')
HSC2454drug.DoPlotXYtrack(label = None, color = 'g')

#%% stretch
trk = HSCnodrug.tracks[53]
print(ContourLength(trk.X_um, trk.Y_um, Rerror = np.sqrt(0.34**2 + 0.25**2)))

#%% stretch
plt.cla()
HSC3809 =  HSCnodrug.DoGetCellTrackSubarrayByParamValue({'animalID': 3809})
HSC3809.DoPlotXYtrack(label = None)

plt.cla()
for trk in HSC3809.tracks:
    trk.DoShowTrackStretches(clearPlot = False)
plt.show()
