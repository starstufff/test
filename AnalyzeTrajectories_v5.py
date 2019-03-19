#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 11 16:00:38 2018

@author: oleg
"""

import sys
sys.path.append('/Users/StarStuff/Desktop/Python files 220119/')

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
ExcelFilePath = '/Users/StarStuff/Dropbox/pythonfiles/intravital measurements_011619.xlsx'
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

#%%
plt.cla()
MPHG1303 =  MPHG.DoGetCellTrackSubarrayByParamValue({'animalID': 1303})
MPHG1303.DoPlotXYtrack()

#%% remove 4113
print(len(MPHG.tracks))
MPHG= MPHG.DoExcludeCellTracksBySubstring({'CompleteMetaName' : '1303'})
print(len(MPHG.tracks))

#%% try to get sense of location errors on the macrophages
# what are the sampling times?

dt = np.array([trk.dt_s for trk in MPHG.tracks])
print(dt.min())
print(dt.max())

aID = [trk.animalID for trk in MPHG.tracks if trk.dt_s < 60]
print(aID)

#%% 3809 seems to be made with the shortest time delay
plt.cla()
MPHG3809 =  MPHG.DoGetCellTrackSubarrayByParamValue({'animalID': 3809})
MPHG3809.DoPlotXYtrack(label = None)
#%% look at the distribution of displacements
trkNos = [0, 1, 2, 3]
dxAll = np.empty((0,))
dyAll = np.empty((0,))
dx2All = np.empty((0,))
dy2All = np.empty((0,))
for trkNo in range(9):
    dx = MPHG3809.tracks[trkNo].X_um[1:] - MPHG3809.tracks[trkNo].X_um[:-1]
    dy = MPHG3809.tracks[trkNo].Y_um[1:] - MPHG3809.tracks[trkNo].Y_um[:-1]
    dxAll = np.concatenate((dxAll, dx))
    dyAll = np.concatenate((dyAll, dy))
    dx = MPHG3809.tracks[trkNo].X_um[2:] - MPHG3809.tracks[trkNo].X_um[:-2]
    dy = MPHG3809.tracks[trkNo].Y_um[2:] - MPHG3809.tracks[trkNo].Y_um[:-2]
    dx2All = np.concatenate((dx2All, dx))
    dy2All = np.concatenate((dy2All, dy))

Hx, bins = np.histogram(dxAll, 50)
Hy, bins = np.histogram(dyAll, bins = bins)
Hx2, bins = np.histogram(dx2All, bins = bins)
Hy2, bins = np.histogram(dy2All, bins = bins)
plt.cla()
plt.plot((bins[:-1]+bins[1:])/2, Hx, '-o', label = 'Hx')
plt.plot((bins[:-1]+bins[1:])/2, Hy, '-o', label = 'Hy')
plt.plot((bins[:-1]+bins[1:])/2, Hx2, '-o', label = 'Hx2')
plt.plot((bins[:-1]+bins[1:])/2, Hy2, '-o', label = 'Hy2')
plt.legend()
plt.show()

#%% standard deviations
print(np.std(dxAll[np.abs(dxAll)<1]))
print(np.std(dyAll[np.abs(dyAll)<1]))
print(np.std(dx2All[np.abs(dx2All)<1]))
print(np.std(dy2All[np.abs(dy2All)<1]))

# let's say in one direction the std is 0.33, then in 2D this gives 0.33*1.4
# and since the vector is defined by two points we need to take at least 2 sigma
res = 0.31*np.sqrt(2)*2
print(res)
# i.e about 1um

#%% testing and test yet lower resolution
HSCnodrug.tracks[0].DoGetContour(coarseRes = [1, 2], XYonly =True)
#%%
for cont in HSCnodrug.tracks[0].contour:
    print(cont.L)
#%%
plt.cla()
for cont in HSCnodrug.tracks[0].contour:
    plt.plot(cont.XYZ_um[0, :],cont.XYZ_um[1, :], '-o')
plt.show()

#%%
#%% testing and test yet lower resolution
HSCnodrug.tracks[0].DoGetContour(coarseRes = [1, 2], XYonly =True)
#%%
HSCnodrug.DoTrackCalculations()

#%% simulate many tracks
NoTimePnts = 400 #20
Ntracks = 100
SimTrks = CellTrackArrayClass()

for i in range(Ntracks):
    trk = CellTrackClass()
    trk.t_s = np.arange(0, NoTimePnts)
    dX = np.random.normal(size = NoTimePnts)
    trk.X_um = np.cumsum(dX) + np.mod(i, 3)*20
    dY = np.random.normal(size = NoTimePnts)
    trk.Y_um = np.cumsum(dY) + np.floor(i/3)*20
    dZ = np.random.normal(size = NoTimePnts)
    trk.Z_um = np.cumsum(dZ)
    SimTrks.tracks.append(trk)
   
SimTrks.DoTrackCalculations()
#%%\

def running_mean(x, N):
    cumsum = np.cumsum(np.insert(x, 0, 0)) 
    return (cumsum[N:] - cumsum[:-N]) / float(N)

def running_var(x, N):
    return running_mean(x**2, N) - running_mean(x, N)**2

def running_contour_length(x, y, N):
    ds = np.sqrt(np.diff(x)**2 + np.diff(y)**2)
    L = np.cumsum(np.insert(ds, 0, 0)) 
    return (L[(N-1):] - L[:-(N-1)])

def running_contour_Rg_ratio(x, y, N):
    Rg = np.sqrt(running_var(x, N) + running_var(y, N))
    L = running_contour_length(x, y, N)
    LRg = L/Rg
    meanderInd = L/np.sqrt((x[(N-1):] - x[:-(N-1)])**2 + (y[(N-1):] - y[:-(N-1)])**2)
    return LRg, L, Rg, meanderInd
    

#%%

WindSize = 19
trkNo = 9

trk = SimTrks.tracks[trkNo]

#trk.DoCalculateVelocities()
#Vangles = running_mean(trk.cosVXYangles, WindSize)

plt.subplot(1, 2, 1)
plt.cla()
trk.DoPlotXYtrack(marker = 'o', label = None, fillstyle = 'none')
plt.plot(trk.contour[1].XYZ_um[0, :], trk.contour[1].XYZ_um[1, :])
plt.subplot(1, 2, 2)
plt.cla()

#cont = trk.contour[0]
#LRg, L, Rg, meanderInd = running_contour_Rg_ratio(cont.XYZ_um[0, :], cont.XYZ_um[0, :], WindSize)
#plt.plot(LRg)

#cont = trk.contour[1]
#LRg, L, Rg, meanderInd = running_contour_Rg_ratio(cont.XYZ_um[0, :], cont.XYZ_um[0, :], WindSize)
#plt.plot(LRg)

#cont = trk.contour[2]
#LRg, L, Rg, meanderInd = running_contour_Rg_ratio(cont.XYZ_um[0, :], cont.XYZ_um[0, :], WindSize)
#plt.plot(LRg)

#plt.plot(running_mean(LRg, WindSize))
#plt.plot(3*meanderInd)
#plt.plot(meanderInd/LRg)
#plt.plot(Vangles)

#NoTimePnts = trk.X_um.size
#trk = CellTrackClass()
#trk.t_s = np.arange(0, NoTimePnts)
#dX = np.random.normal(size = NoTimePnts)
#trk.X_um = np.cumsum(dX)
#dY = np.random.normal(size = NoTimePnts)
#trk.Y_um = np.cumsum(dY)
#dZ = np.random.normal(size = NoTimePnts)
#trk.Z_um = np.cumsum(dZ)
#LRg, L, Rg, meanderInd = running_contour_Rg_ratio(trk.X_um, trk.Y_um, WindSize)
#trk.DoCalculateVelocities()
Vangles = running_mean(trk.contour[1].cosTangent, WindSize)
#Vangles = running_mean(trk.cosVXYangles, WindSize)


#plt.plot(LRg)
#plt.plot(running_mean(LRg, WindSize))
#plt.plot(3*meanderInd)
#plt.plot(meanderInd/LRg)
plt.plot(Vangles)
plt.show()

#%% are there HSCs in 3809: shortest time delay
plt.cla()
HSC3809 =  HSCnodrug.DoGetCellTrackSubarrayByParamValue({'animalID': 3809})
HSC3809.DoPlotXYtrack(marker = 'o')

#%%
def contour_length(x, y, coarseGrain = 1):
    diffX = x[coarseGrain:] - x[:-coarseGrain] 
    diffY = y[coarseGrain:] - y[:-coarseGrain] 
    ds = np.sqrt(diffX**2 + diffY**2)
    L = ds.sum()/coarseGrain
    ds2 = np.mean(ds**2)
    ds2a
    return L, ds.mean(), ds2

#%%
trk = HSC3809.tracks[1]
N = 1
L, ds = contour_length(trk.X_um, trk.Y_um, N)
print(L, ds, ds/N)
N = 2
L, ds = contour_length(trk.X_um, trk.Y_um, N)
print(L, ds, ds/N)
N = 3
L, ds = contour_length(trk.X_um, trk.Y_um, N)
print(L, ds, ds/N)
N = 4
L, ds = contour_length(trk.X_um, trk.Y_um, N)
print(L, ds, ds/N)


#%%
WindSize = 9
trkNo = 2

trk = HSC3809.tracks[trkNo]

#trk.DoCalculateVelocities()
#Vangles = running_mean(trk.cosVXYangles, WindSize)

plt.subplot(1, 2, 1)
plt.cla()
trk.DoPlotXYtrack(marker = 'o', label = None, fillstyle = 'none')
plt.plot(trk.contour[1].XYZ_um[0, :], trk.contour[1].XYZ_um[1, :])
plt.subplot(1, 2, 2)
plt.cla()

#cont = trk.contour[0]
#LRg, L, Rg, meanderInd = running_contour_Rg_ratio(cont.XYZ_um[0, :], cont.XYZ_um[0, :], WindSize)
#plt.plot(LRg)

#cont = trk.contour[1]
#LRg, L, Rg, meanderInd = running_contour_Rg_ratio(cont.XYZ_um[0, :], cont.XYZ_um[0, :], WindSize)
#plt.plot(LRg)

#cont = trk.contour[2]
#LRg, L, Rg, meanderInd = running_contour_Rg_ratio(cont.XYZ_um[0, :], cont.XYZ_um[0, :], WindSize)
#plt.plot(LRg)

#plt.plot(running_mean(LRg, WindSize))
#plt.plot(3*meanderInd)
#plt.plot(meanderInd/LRg)
#plt.plot(Vangles)

#NoTimePnts = trk.X_um.size
#trk = CellTrackClass()
#trk.t_s = np.arange(0, NoTimePnts)
#dX = np.random.normal(size = NoTimePnts)
#trk.X_um = np.cumsum(dX)
#dY = np.random.normal(size = NoTimePnts)
#trk.Y_um = np.cumsum(dY)
#dZ = np.random.normal(size = NoTimePnts)
#trk.Z_um = np.cumsum(dZ)
#LRg, L, Rg, meanderInd = running_contour_Rg_ratio(trk.X_um, trk.Y_um, WindSize)
#trk.DoCalculateVelocities()
Vangles = running_mean(trk.contour[1].cosTangent, WindSize)
#Vangles = running_mean(trk.cosVXYangles, WindSize)


#plt.plot(LRg)
#plt.plot(running_mean(LRg, WindSize))
#plt.plot(3*meanderInd)
#plt.plot(meanderInd/LRg)
plt.plot(Vangles)
plt.show()

#%% 
dt = np.array([trk.dt_s for trk in HSCnodrug.tracks])
aID = [trk.animalID for trk in HSCnodrug.tracks]

#%% 1303 has 20s dt
plt.cla()
HSC1303 =  HSCnodrug.DoGetCellTrackSubarrayByParamValue({'animalID': 1303})
HSC1303.DoPlotXYtrack(marker = 'o')

#%%
trk = HSC1303.tracks[0]
N = 1
L, ds = contour_length(trk.X_um, trk.Y_um, N)
print(L, ds, ds/N)
N = 2
L, ds = contour_length(trk.X_um, trk.Y_um, N)
print(L, ds, ds/N)
N = 3
L, ds = contour_length(trk.X_um, trk.Y_um, N)
print(L, ds, ds/N)
N = 6
L, ds = contour_length(trk.X_um, trk.Y_um, N)
print(L, ds, ds/N)

#%% 5009 has has one at 46s and another at 237
plt.cla()
HSC5009 =  HSCnodrug.DoGetCellTrackSubarrayByParamValue({'animalID': 5009})
HSC5009.DoPlotXYtrack(marker = 'o')

#%%
trk = HSC5009.tracks[1]
print(trk.dt_s)
N = 1
L, ds = contour_length(trk.X_um, trk.Y_um, N)
print(L, ds, ds/N)
N = 2
L, ds = contour_length(trk.X_um, trk.Y_um, N)
print(L, ds, ds/N)
N = 3
L, ds = contour_length(trk.X_um, trk.Y_um, N)
print(L, ds, ds/N)
N = 4
L, ds = contour_length(trk.X_um, trk.Y_um, N)
print(L, ds, ds/N)

#%% look at the other side of the spectrum about 15 min dts
plt.cla()
HSC3111 =  HSCnodrug.DoGetCellTrackSubarrayByParamValue({'animalID': 3111})
HSC3111.DoPlotXYtrack(marker = 'o')

#%%
trk = HSC3111.tracks[3]
print(trk.dt_s)
N = 1
L, ds = contour_length(trk.X_um, trk.Y_um, N)
print(L, ds, ds/N)
N = 2
L, ds = contour_length(trk.X_um, trk.Y_um, N)
print(L, ds, ds/N)
N = 3
L, ds = contour_length(trk.X_um, trk.Y_um, N)
print(L, ds, ds/N)
N = 4
L, ds = contour_length(trk.X_um, trk.Y_um, N)
print(L, ds, ds/N)

#%% create mesh of distances
N = 1
L, ds, ds2 = contour_length(trk.X_um, trk.Y_um, N)
X1, X2 = np.meshgrid(trk.X_um, trk.X_um)
Y1, Y2 = np.meshgrid(trk.Y_um, trk.Y_um)
Rsq = (X1 - X2)**2 + (Y1 - Y2)**2
I, J = np.meshgrid(np.arange(trk.X_um.size), np.arange(trk.X_um.size))
rat = Rsq/(ds2*np.abs(I-J))
logP = -rat + np.log(trk.X_um.size-np.abs(I-J))
plt.cla()
plt.imshow(-logP)
plt.show()

#%%
plt.cla()
trk.DoPlotXYtrack(marker = 'o')
plt.plot(trk.X_um[16:20], trk.Y_um[16:20])
plt.plot(trk.X_um[20:26], trk.Y_um[20:26])
plt.plot(trk.X_um[4:13], trk.Y_um[4:13])
plt.show()

#%%
from scipy import ndimage as ndi
from skimage.feature import peak_local_max
from skimage import data, img_as_float

P[np.isnan(P)] = 0
coordinates = peak_local_max(P, min_distance=1)
print(coordinates)

#%%
plt.cla()
HSC5009.DoPlotXYtrack(label = None)
plt.show()

#%%
plt.cla()
HSC4107 =  HSCnodrug.DoGetCellTrackSubarrayByParamValue({'animalID': 4107})
HSC4107.DoPlotXYtrack(label = None)

#%%
plt.cla()
HSC4107.DoPlotMSD(MSDtype = 'XYsq_um2', label = None)

#%%
plt.cla()
HSCnodrug.DoPlotMSD(MSDtype = 'XYsq_um2', label = None)

#%%
plt.cla()
HSCnodrug.DoFitMSD(fitfunc = 'BallisticDiffusion2DFit', ShowXYtrack = True)

#%%
plt.cla()
MPHG.DoFitMSD(fitfunc = 'BallisticDiffusion2DFit', ShowXYtrack = True)

#%% simulate many tracks for escape probability 
NoTimePnts = 100 #20
Ntracks = 10000
SimTrks = CellTrackArrayClass()

for i in range(Ntracks):
    trk = CellTrackClass()
    trk.t_s = np.arange(0, NoTimePnts)
    dX = np.random.normal(size = NoTimePnts)
    trk.X_um = np.cumsum(dX)
    dY = np.random.normal(size = NoTimePnts)
    trk.Y_um = np.cumsum(dY)
    dZ = np.random.normal(size = NoTimePnts)
    trk.Z_um = np.cumsum(dZ)
    SimTrks.tracks.append(trk)
    
#%%
Rsq = np.array([(trk.X_um[-1]**2 + trk.Y_um[-1]**2) for trk in SimTrks.tracks])
Rsqm = [(trk.X_um[-2]**2 + trk.Y_um[-2]**2) for trk in SimTrks.tracks]
Rsqm = np.array( Rsqm + [((trk.X_um[-1] - trk.X_um[0])**2 + (trk.Y_um[-1] - trk.Y_um[0])**2) for trk in SimTrks.tracks])
#%%

R = np.arange(0, np.sqrt(Rsq.max()), np.sqrt(Rsq.max())/30)
Pp = np.array([(Rsq > r**2).sum()/Rsq.size for r in R])
Ppm = np.array([(Rsqm > r**2).sum()/Rsqm.size for r in R])
plt.cla()
plt.semilogy(R**2, Pp, 'o')
plt.semilogy(R**2, Ppm, 'o')
plt.semilogy(R**2, np.exp(-R**2/(NoTimePnts +1)/2))
plt.show()

#%% to what degree consecutive segments are independent statistics
NoTimePnts1 = 50
Rsqm = []
for NN  in range(NoTimePnts-NoTimePnts1):
    Rsqm = Rsqm + [((trk.X_um[NN +NoTimePnts1] - trk.X_um[NN])**2 + 
                    (trk.Y_um[NN + NoTimePnts1] - trk.Y_um[NN])**2) for trk in SimTrks.tracks]

Rsqm = np.array(Rsqm)   

#%%
R = np.arange(0, np.sqrt(Rsqm.max()), np.sqrt(Rsqm.max())/30)
Pp = np.array([(Rsqm > r**2).sum()/Rsqm.size for r in R])
plt.cla()
plt.plot(R, Pp, 'o')
plt.plot(R, np.exp(-R**2/(NoTimePnts1 +1)/2))
plt.show()

#%% practice getting sum of ds2 between each step
x = trk.X_um
y = trk.Y_um
#x = np.arange(5)
#y = np.arange(5)
diffX = np.diff(x)
diffY = np.diff(y)
ds2 = diffX**2 + diffY**2
cumsumDs2 = np.cumsum(np.insert(ds2, 0, 0))
CS1, CS2 = np.meshgrid(cumsumDs2, cumsumDs2)
s2 = np.abs(CS1 - CS2)
plt.cla()
plt.imshow(s2)
plt.show()


#%%
trk = HSC3111.tracks[3]
X1, X2 = np.meshgrid(x, x)
Y1, Y2 = np.meshgrid(y, y)
Rsq = (X1 - X2)**2 + (Y1 - Y2)**2
I, J = np.meshgrid(np.arange(x.size), np.arange(x.size))
plt.subplot(1, 2, 1)
rat = Rsq/s2
logP = -rat + np.log(trk.X_um.size-np.abs(I-J))
plt.cla()
plt.imshow(-logP)
plt.subplot(1, 2, 2)
rat = Rsq/(ds2.mean()*np.abs(I-J))
logP = -rat + np.log(trk.X_um.size-np.abs(I-J))
plt.cla()
plt.imshow(-logP)
plt.show()

#%% extract local minima
logP[np.isnan(logP)] = logP.max()
temp = -np.tril(logP)
coord = peak_local_max(temp, threshold_abs = 0)
#print(logP[coord[:, 0], coord[:, 1]])

plt.cla()
plt.imshow(temp)
plt.show()

ii = np.argsort(-temp[coord[:, 0], coord[:, 1]]) #- sign to have a descending sort
print(temp[coord[ii, 0], coord[ii, 1]])

#%%
trk = HSC3111.tracks[2]
#trk = HSC5009.tracks[0]
Llim = 2
x = trk.X_um
y = trk.Y_um
X1, X2 = np.meshgrid(x, x)
Y1, Y2 = np.meshgrid(y, y)
Rsq = (X1 - X2)**2 + (Y1 - Y2)**2
I, J = np.meshgrid(np.arange(x.size), np.arange(x.size))

diffX = np.diff(x)
diffY = np.diff(y)
ds2 = diffX**2 + diffY**2
cumsumDs2 = np.cumsum(np.insert(ds2, 0, 0))
CS1, CS2 = np.meshgrid(cumsumDs2, cumsumDs2)
s2 = np.abs(CS1 - CS2)

plt.subplot(1, 2, 1)
rat = Rsq/s2
rat = Rsq/(ds2.mean()*np.abs(I-J))
logP = -rat + np.log(trk.X_um.size-np.abs(I-J))
plt.cla()
plt.imshow(-logP)

plt.subplot(1, 2, 2)
logP[np.isnan(logP)] = logP.max()
temp = -logP
coord = peak_local_max(temp, threshold_abs = -2, exclude_border = False)
#print(logP[coord[:, 0], coord[:, 1]])

ii = np.argsort(-temp[coord[:, 0], coord[:, 1]]) #- sign to have a descending sort
ii = ii[::2]
print(temp[coord[ii, 0], coord[ii, 1]])
coord = coord[ii, :]
coord.sort(axis = 1)
#dismiss short stretches
coord = coord[(np.diff(coord, axis =1) > Llim).T[0]]
#%% recalibrate probabilites by orthogonal deviation
c = coord[6]
xx = x[c[0]:(c[1]+1)]
yy = y[c[0]:(c[1]+1)]
print(xx.size)
JJ = [xx[-1] - xx[0], yy[-1] - yy[0]]
JJ = JJ/np.linalg.norm(JJ)
JJP = [-JJ[1], JJ[0]]
rrP= np.array([xx - xx[0], yy - yy[0]]).T@JJP
rrPmax = np.square(rrP).max()
rrPth = (xx.size - 1)/8*ds2.mean()

PP = erf(np.sqrt(rrPmax/rrPth/2))
print(PP)





#%%
c = coord[0]
isIn = [c in cc for cc in coord[:]]
print(isIn)
#%%
def intervalOverlap(interval1, intervalArray):
    center = interval1.mean()
    rad = np.abs(np.diff(interval1)[0])/2
    centers = intervalArray.mean(axis =1)
    rads = np.abs(np.diff(intervalArray, axis = 1)[:, 0])/2
    return np.abs(centers - center) < (rad + rads)
#remove overlapping
coord2 = coord
coord1 = []
for c in coord:
    if c in coord2:
        print(c)
        isIn = intervalOverlap(c, coord2[:]) 
        coord1.append(coord2[isIn][0])
        notIsIn = [not s for s in isIn]
        coord2 = coord2[notIsIn]
coord1 = np.array(coord1)  
    
#%%
plt.cla()
trk.DoPlotXYtrack(marker = 'o')
for c in coord1:
    plt.plot(trk.X_um[c[0]:(c[1]+1)], trk.Y_um[c[0]:(c[1]+1)], marker = 'o')
plt.show()

#%% probabilities in terms of angles
trk.DoCalculateVelocities()
phi = np.arccos(trk.cosVXYangles)


#%% summarize automated annotation of the tracks
#trk = HSC3111.tracks[2]
#trk = HSC5009.tracks[0]
trk = SimTrks.tracks[7]
maxProbability = 0.1
dr2_error = 0.34**2 + 0.25**2
#

def intervalOverlap(interval1, intervalArray):
    center = interval1.mean()
    rad = np.abs(np.diff(interval1)[0])/2
    centers = intervalArray.mean(axis =1)
    rads = np.abs(np.diff(intervalArray, axis = 1)[:, 0])/2
    return np.abs(centers - center) < (rad + rads)


def FindDirectedMotion(trk, maxProbability = 0.1, minSteps = 2, Rerror = np.sqrt(0.34**2 + 0.25**2),
                       showPlots = True, RefineProb = True):
    DirectStretchList = []
    Llim = minSteps
    dr2_error = Rerror**2
    x = trk.X_um
    y = trk.Y_um
    X1, X2 = np.meshgrid(x, x)
    Y1, Y2 = np.meshgrid(y, y)
    Rsq = (X1 - X2)**2 + (Y1 - Y2)**2
    I, J = np.meshgrid(np.arange(x.size), np.arange(x.size))
    
    diffX = np.diff(x)
    diffY = np.diff(y)
    ds2 = diffX**2 + diffY**2
    # taking the average of ds on the stretch itself
    #cumsumDs2 = np.cumsum(np.insert(ds2, 0, 0))
    #CS1, CS2 = np.meshgrid(cumsumDs2, cumsumDs2)
    #s2 = np.abs(CS1 - CS2)
    #rat = Rsq/s2
        
    rat = Rsq/((ds2.mean()-dr2_error)*np.abs(I-J))
    logP = -rat + np.log((trk.X_um.size-1)/(np.abs(I-J)-1)) # - 4*np.abs(I-J)/trk.X_um.size
    
    if showPlots:
        plt.cla()
        plt.subplot(1, 2, 1)
        plt.imshow(-logP)
        plt.subplot(1, 2, 2)
        plt.cla()
        trk.DoPlotXYtrack(marker = 'o')

        
        
    logP[np.isnan(logP)] = logP.max()
    temp = -logP
    coord = peak_local_max(temp, threshold_abs = -np.log(3*maxProbability), exclude_border = False)
    if coord.size == 0:
        return DirectStretchList
    #print(logP[coord[:, 0], coord[:, 1]])
    #print(coord.size)
    ii = np.argsort(-temp[coord[:, 0], coord[:, 1]]) #- sign to have a descending sort
    ii = ii[::2] # remove same peaks from the other matrix triangle
#    print(coord)
    coord = coord[ii, :]
#    print(coord)
    coord.sort(axis = 1) #make sure that step indexes in the move are arrange in ascending order
#    print(coord)
    #dismiss short stretches
    coord = coord[(np.diff(coord, axis =1) >= Llim).T[0]]
    if coord.size == 0:
        plt.show()
        return DirectStretchList
    #print(np.exp(logP[coord[:, 0], coord[:, 1]]))
     
    # remove overlaping intervals 
    coord2 = coord
    #print(np.exp(logP[coord2[:, 0], coord2[:, 1]]))
    coord1 = []
    for c in coord:
        ind = np.where((coord2[:, 0]  == c[0]) * (coord2[:, 1]  == c[1]))
        #if c in coord2:
        if ind[0].size > 0:
            #print(c)
            isIn = intervalOverlap(c, coord2[:]) 
            coord1.append(coord2[isIn][0])
            notIsIn = [not s for s in isIn]
            coord2 = coord2[notIsIn]   
    coord1 = np.array(coord1)  
    if coord1.size == 0:
        plt.show()
        return DirectStretchList

    Probabls = np.exp(logP[coord1[:, 0], coord1[:, 1]])
    # recalibrate probabilites by orthogonal deviation
    if RefineProb:
        PP = []
        for c in coord1:
            xx = x[c[0]:(c[1]+1)]
            yy = y[c[0]:(c[1]+1)]
    #        print(xx.size)
            JJ = [xx[-1] - xx[0], yy[-1] - yy[0]]
            JJ = JJ/np.linalg.norm(JJ)
            JJP = [-JJ[1], JJ[0]]
            rrP= np.array([xx - xx[0], yy - yy[0]]).T@JJP
            rrPmax = np.square(rrP).max()
            rrPth = (xx.size - 1)/8*(ds2.mean()-dr2_error)
            
            PP.append(erf(np.sqrt(rrPmax/rrPth/2)))  
        Probabls = np.array(PP)*Probabls
    
        
    coord1 = coord1[Probabls <= maxProbability]
    Probabls = Probabls[Probabls <= maxProbability]
    if coord1.size == 0:
        plt.show()
        return DirectStretchList
    
#    DirectStretchList.append({'Stretch Indices': coord1, 'Probabilities': Probabls, 'Stretch Lengths': np.diff(coord1, axis = 1)[0]})
    # plotting
    if showPlots:
        plt.cla()
        trk.DoPlotXYtrack(marker = 'o')
        for c in coord1:
            plt.plot(trk.X_um[c[0]:(c[1]+1)], trk.Y_um[c[0]:(c[1]+1)], marker = 'o')
        plt.show()
    
    for c, P in zip(coord1, Probabls):
        DirectStretchList.append({'Stretch Indices': c, 
                                  'Probabilities': P, 
                                  'Stretch Lengths': np.diff(c)[0],
                                  'X': trk.X_um[c[0]:(c[1]+1)],
                                  'Y': trk.Y_um[c[0]:(c[1]+1)]})
        
    return DirectStretchList


#%%
AllStretches = []
#%%

trk = SimTrks.tracks[3]   
DS = FindDirectedMotion(trk, maxProbability = 0.1, minSteps = 4, Rerror = np.sqrt(0.34**2 + 0.25**2),
                       showPlots = False)
AllStretches = AllStretches + DS
print(DS)

#%%
AllStretches = []
for ind, trk in enumerate(SimTrks.tracks[:1000]):
    DS = FindDirectedMotion(trk, maxProbability = 0.01, minSteps = 4, Rerror = 0,
                       showPlots = False, RefineProb = False)
    for ds in DS:
        if ds['Stretch Lengths'] == 99:
            print(ind)
    AllStretches = AllStretches + DS


LL = []
PP = []
for Stretch in AllStretches:
    LL.append(Stretch['Stretch Lengths'])
    PP.append(Stretch['Probabilities'])
LL = np.array(LL)
PP = np.array(PP)

H, bins = np.histogram(LL, bins = np.arange(1, 100))
plt.subplot(1, 1, 1)
plt.cla()
plt.plot(bins[:-1], H, '-o')
plt.show()

#%%
plt.cla()
trk = HSCnodrug.tracks[10]
DS = FindDirectedMotion(trk, maxProbability = 0.01, minSteps = 4, 
                    Rerror = np.sqrt(0.34**2 + 0.25**2),
                       showPlots = True, RefineProb = False)

for ds in DS:
    print(ds['Stretch Lengths'])
    print(ds['Probabilities'])


#%%
NoTimePnts = 20 #trk.X_um.size #20
Ntracks = 10000
SimTrks = CellTrackArrayClass()

for i in range(Ntracks):
    trk = CellTrackClass()
    trk.t_s = np.arange(0, NoTimePnts)
    dX = np.random.normal(size = NoTimePnts)
    trk.X_um = np.cumsum(dX) + np.mod(i, 3)*20
    dY = np.random.normal(size = NoTimePnts)
    trk.Y_um = np.cumsum(dY) + np.floor(i/3)*20
    dZ = np.random.normal(size = NoTimePnts)
    trk.Z_um = np.cumsum(dZ)
    SimTrks.tracks.append(trk)
 
#%%
AllStretches = []
for ind, trk in enumerate(SimTrks.tracks):
    DS = FindDirectedMotion(trk, maxProbability = 0.04, minSteps = 3, Rerror = 0,
                       showPlots = False, RefineProb = False)
    for ds in DS:
        ds['trackInd'] = ind
    AllStretches = AllStretches + DS


LL = []
PP = []
Ind = []
for Stretch in AllStretches:
    LL.append(Stretch['Stretch Lengths'])
    PP.append(Stretch['Probabilities'])
    Ind.append(Stretch['trackInd'])
    
LL = np.array(LL)
PP = np.array(PP)
Ind = np.array(Ind)

H, bins = np.histogram(LL, bins = np.arange(1, LL.max()))
plt.subplot(1, 1, 1)
plt.cla()
plt.plot(bins[:-1], H, '-o')
plt.show()
#%%
Hist_L200_P01 = dict()
Hist_L200_P01['bins'] = bins
Hist_L200_P01['H'] = H
Hist_L200_P01['L'] = 200
Hist_L200_P01['maxP'] = 0.1
#%% 
plt.cla()
plt.semilogy(bins[:-1], H, '-o')
plt.show()
#%%
plt.cla()
HH = Hist_L100_P01
plt.semilogy(HH['bins'][:-1], HH['H']/HH['H'].sum()*np.exp(4*HH['bins'][:-1]/HH['L']))
print(HH['H'].sum())
HH = Hist_L100_P003
plt.semilogy(HH['bins'][:-1], HH['H']/HH['H'].sum())
print(HH['H'].sum())
HH = Hist_L100_P001
plt.semilogy(HH['bins'][:-1], HH['H']/HH['H'].sum())
print(HH['H'].sum())
HH = Hist_L200_P001
plt.semilogy(HH['bins'][:-1]/2, 2*HH['H']/HH['H'].sum())
print(HH['H'].sum())
HH = Hist_L200_P01
plt.semilogy(HH['bins'][:-1]/2, 2*HH['H']/HH['H'].sum())
print(HH['H'].sum())
plt.semilogy(bins[:-1], H/H.sum(), '-o')
plt.show()


#%%
plt.cla()
trk = SimTrks.tracks[17]
DS = FindDirectedMotion(trk, maxProbability = 0.036, minSteps = 4, Rerror = 0,
                showPlots = True, RefineProb = False)


#%%
coord2 = coord
coord1 = []
for c in coord[:1]:
    if c in coord2:
        print(c)
        isIn = intervalOverlap(c, coord2[:]) 
        coord1.append(coord2[isIn][0])
        notIsIn = [not s for s in isIn]
        coord2 = coord2[notIsIn]   
coord1 = np.array(coord1)  

#%%
for ds in DS:
    print(ds['Stretch Indices'])
#%%   
tt0 = np.where(coord2[:, 0] == 114) 
tt1 = np.where(coord2[:, 1] == 402)
tt = np.intersect1d(tt0, tt1)
print(tt)
tt2 = np.where(coord2[:, 0] == 18) 
tt3 = np.where(coord2[:, 1] == 37)
ttt = np.intersect1d(tt2, tt3)
print(ttt)
#%%
ind = np.where((coord2[:, 0]  == 114) * (coord2[:, 1]  == 402))

#%%
np.where((coord2[:, 0]  == c[0]) * (coord2[:, 1]  == c[1]))

#%% what are typical track lengths
trkLen = np.array([trk.X_um.size for trk in HSCnodrug.tracks])
H, bins = np.histogram(trkLen)
plt.cla()
plt.plot(bins[:-1], H, '-o')
plt.show()

#%%
trkNo = 5
HSCnodrug.tracks[trkNo].DoFindDirectedMotion()
HSCnodrug.tracks[trkNo].DoShowTrackStretches()

#%%
HSCnodrug.DoFindDirectedMotion(maxProbability = 0.05, minSteps = 3, 
                    Rerror = np.sqrt(0.34**2 + 0.25**2),
                       showPlots = False, RefineProb = False)

#%%
HSCnodrug.DoShowTrackStretches()

#%%
trk = HSCnodrug.tracks[53]
print(ContourLength(trk.X_um, trk.Y_um, Rerror = np.sqrt(0.34**2 + 0.25**2)))

#%%
Ttotal = []
T =[]
for trk in HSCnodrug.tracks:
    if trk.maxDist >=20 :
        Ttotal.append(trk.t_s[-1])
        for stretch in trk.DirectStretchList:
            T.append(stretch['Stretch Time (s)'])
        
T = np.array(T)
Ttotal = np.array(Ttotal)
print(T.mean()/60)
print(T.sum()/Ttotal.sum())    


#%% 3809 seems to be made with the shortest time delay
plt.cla()
HSC633 =  HSCnodrug.DoGetCellTrackSubarrayByParamValue({'animalID': 633})
HSC633.DoPlotXYtrack(label = None)

plt.cla()
for trk in HSC633.tracks:
    trk.DoShowTrackStretches(clearPlot = False)
plt.show()