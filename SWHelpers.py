import os
import glob
import imp
from copy import copy,deepcopy

import numpy as np
from scipy import ndimage
import wx
import Image
import matplotlib.pyplot as plt

from cmpGen import cmpGen
from InstantRead import read
from openCSV import openCSV
import FilenameSort as FS
import ExcelHelper
import ImageContour
from ValueReceived import imshow_vr
from dplot import dplot
import GifTiffLoader as GTL
import SeedWaterSegmenter as SWS

from MyPath import cprint
#from pylab import *

def norm(a,framesExcludeAfter=None):
    if framesExcludeAfter==None:
        return a/np.mean(a)
    else:
        return a/np.mean(a[:framesExcludeAfter])
def deletecases(l,valList):
    if not hasattr(valList,'__iter__'):
        valList=[valList]
    l=copy(l) # work on a copy
    for v in valList:
        for i in range(l.count(v)):
            l.remove(v)
    return l

def GetWoutline(watershed,dilate):
    woutline = np.zeros(watershed.shape,np.int)
    woutline[:-1,:-1] = (np.diff(watershed,axis=0)[:,:-1]!=0) + (np.diff(watershed,axis=1)[:-1]!=0)
    
    if dilate>0:
        woutline = ndimage.binary_dilation(woutline,iterations=dilate)
    
    return woutline

def GetSeedArray(watershed,SeedsPyFile,frameNumber):
    # Load the SeedsPyFile:
    fid = open(SeedsPyFile)
    try:
        exec(fid.read().replace('\r','')) # Loads seedList and seedVals from the file...
    except SyntaxError:
        print 'Invalid Syntax in Seeds.py'        
    fid.close()
    
    # Create a seedArray, Fill using seeds and vals
    seedArray = np.zeros(watershed.shape,np.int)
    
    if seedList[frameNumber]!=None:
        SWS.PointsToArray(seedList[frameNumber], seedVals[frameNumber],seedArray)
    val = 2 #(2 if self.walgorithm=='cv' else 1) # I could do this, but what's the need?
    seedArray[:val,:]  = 1
    seedArray[-val:,:] = 1
    seedArray[:,:val]  = 1
    seedArray[:,-val:] = 1
    
    return seedArray

def CreateOverlay(original,watershed,SeedsPyFile,frameNumber,dilate=0,returnSeparateArrays=False,overlaySolid=False):
    # Generate the outline file from the watershed:
    len(original.shape)
    
    woutline=GetWoutline(watershed,dilate)
    seedArray=GetSeedArray(watershed,SeedsPyFile,frameNumber)
    
    seed = (seedArray!=0).astype(np.uint8)
    outl = woutline.astype(np.uint8)
    seedORoutl = np.array(np.array(seed+outl)>0,np.uint8)
    
    overlayIm = GTL.DivideConvertType(original,8)
    if returnSeparateArrays:
        return overlayIm,outl,seed
    else:
        overlayIm = np.array([overlayIm,overlayIm,overlayIm]).transpose(1,2,0)
        
        overlayIm[:,:,0][np.where(outl==True)] = 255 # outl points turn red
        overlayIm[:,:,1][np.where(seed==True)] = 255 # seeds point turn green
        
        if overlaySolid:
            overlayIm[:,:,1][np.where(outl==True)] = 0 # outl points turn red
            overlayIm[:,:,2][np.where(outl==True)] = 0 # outl points turn red
            overlayIm[:,:,0][np.where(seed==True)] = 0 # seeds point turn green
            overlayIm[:,:,2][np.where(seed==True)] = 0 # seeds point turn green
        
        return overlayIm
def CreateMapPlot(watershedI):
    np.random.seed(0)
    mapPlotRandomArray=np.array([np.random.random(10000)*236+20,
                                 np.random.random(10000)*236+20,
                                 np.random.random(10000)*236+20], dtype=np.uint8)
    mapPlotRandomArray[:,0]=255
    mapPlotRandomArray[:,1]=255
    
    f1=np.vectorize(lambda x: mapPlotRandomArray[0][x])
    f2=np.vectorize(lambda x: mapPlotRandomArray[1][x])
    f3=np.vectorize(lambda x: mapPlotRandomArray[2][x])
    w = watershedI
    rgbMi = np.array([f1(w),f2(w),f3(w)]).transpose(1,2,0)
    return rgbMi
def CreateMapPlots(watershed):
    np.random.seed(0)
    mapPlotRandomArray=np.array([np.random.random(10000)*236+20,
                                 np.random.random(10000)*236+20,
                                 np.random.random(10000)*236+20], dtype=np.uint8)
    mapPlotRandomArray[:,0]=255
    mapPlotRandomArray[:,1]=255
    
    f1=np.vectorize(lambda x: mapPlotRandomArray[0][x])
    f2=np.vectorize(lambda x: mapPlotRandomArray[1][x])
    f3=np.vectorize(lambda x: mapPlotRandomArray[2][x])
    rgbM=[]
    for w in watershed:
        rgbM.append( np.array([f1(w),f2(w),f3(w)]).transpose(1,2,0) )
    return rgbM
def CreateDiffImage(arr,targetArr):
    return ((arr - targetArr)!=0).astype(np.uint8)
def GetTimeAxisAdj():
    sList = []
    timeAdj = 0
    count=0
    for i in FS.getSortedListOfFiles(di,'*.tif'):
        count+=1
        d,h,m,s=map(int,os.path.split(i)[1].split('_')[2:6])
        h += d*24
        m += h*60
        s += m*60 + timeAdj
        if sList!=[]:
            if s - sList[-1]>80: # Pretend any large gaps in time are only 1min
                print count
                timeAdj += sList[-1] - s + 80
                s = sList[-1] + 80
            pass
        sList.append(s)
    return sList

def SaveOverlayImages(tifStackData,SeedsPyFile,SegmentsDirectory,red=True):
    # This function is more-or-less just a template right now, has MAJOR issues...
    r = GTL.LoadFileSequence(di, wildcard='*.tif')
    rgbM = CreateMapPlots(r)
    
    for i in range(len(Files)):
        tifDir,nameTif = os.path.split(Files[i])
        mapDir = os.path.join(tifDir,'Maps')
        overlayDir = os.path.join(tifDir,'Overlay')
        for j in [mapDir,overlayDir]:
            if not os.path.exists(j):
                os.mkdir(j)
        name = os.path.splitext(nameTif)[0]
        nameSeedsPy = '_'.join(name.split('_')[:-2])+'_Seeds.py'
        ### GENERATE MAP PLOT
        namePng = name+'.png'
        rgbMi = CreateMapPlot(r[i])
        im = Image.fromarray(rgbMi)
        im.save(os.path.join(mapDir,namePng))
        
        ### NOW GENERATE OVERLAY IMAGE
        # Grab the frame number from the image name...
        frameNumber = int(name.split('_')[-1])
        # Get the watershed file from the saved tif:
        watershed = GTL.LoadSingle(Files[i])
        # Find and load the correct original array based on the frame number:
        originalDir = '.../2010FEB09_Edge5/TiffStack' # Put in the right directory here...
        origFile = FS.getSortedListOfFiles(originalDir)[frameNumber]
        original = GTL.LoadSingle(origFile)
        
        SeedsPyFile = os.path.join(tifDir,nameSeedsPy)
        overlayIm = CreateOverlay(original,watershed,SeedsPyFile,frameNumber)
        
        # Now, save the overlay as a png
        nameSeedsOutl = name+'SO.png'
        im = Image.fromarray(overlayIm)
        im.save(os.path.join(overlayDir,nameSeedsOutl))
        cprint(str(i))    
def SaveDiffImages(logDir):
    diffDir = os.path.join(logDir,'Diff')
    if not os.path.exists(diffDir):
        os.mkdir(diffDir)
    Files = FS.getSortedListOfFiles(logDir,'*.tif')
    frameUpdateList = [] # [time,fn,newVal]
    for i,f in enumerate(Files):
        name = os.path.splitext(os.path.split(f)[1])[0]
        frameNumber = int(name.split('_')[-1])
        cur = GTL.LoadSingle(f)
        final = GTL.LoadSingle(FrameFiles[frameNumber])
        d = CreateDiffImage(cur,final)
        frameUpdateList.append([timeAxis[i],frameNumber,d.astype(int).sum()])
        im = Image.fromarray(d*255)
        im.save(os.path.join(diffDir,name+'.png'))

# operates on the list in place
def MergeWoundVals(neigh,wvl):
    firstWV = wvl[0]

    if neigh[firstWV-2]==None:
        neigh[firstWV-2]=[]
    # Merge all the wound value's neighbor lists together
    for wv in sorted(wvl[1:])[::-1]:
        if neigh[wv-2]!=None:
            neigh[firstWV-2] += neigh[wv-2]
        neigh[wv-2] = None
    # and remove any duplicate values this created (and wound self-references)
    neigh[firstWV-2] = list(set(deletecases(neigh[firstWV-2],wvl)))
    # Now, go through the whole list and convert references for other wv's to firstWV, removing duplicates
    for i in range(len(neigh)):
        if neigh[i] != None:
            neigh[i] = sorted(list(set([(firstWV if (j in wvl) else j) for j in neigh[i]])))

def GetNeighborPairsByFrame(neighbors,woundVals): # pass Neighbors.neighbors
    neighborPairsByFrame = []
    for frame in range(len(neighbors)):
        MergeWoundVals(neighbors[frame],woundVals)
        neighborPairs = set()
        for i in range(len(neighbors[frame])):
            if neighbors[frame][i]!=None:
                if len(neighbors[frame][i])>0:
                    neighborPairs.update( [tuple(sorted([i+2,j])) for j in neighbors[frame][i] ] )
        neighborPairsByFrame.append(sorted(list(neighborPairs)))
    return neighborPairsByFrame

def GetAllValsByFrame(neighborPairsByFrame):
    allValsByFrame = []
    for neighborPairs in neighborPairsByFrame:
        allValsByFrame.append(sorted(list(set([j for i in neighborPairs for j in i]))))
    return allValsByFrame

def GetWoundNeighborSetsByFrame(neighborPairsByFrame,wv):
    woundNeighborsByFrame = []
    woundNeighbors2ByFrame = []
    secondNeighborsByFrame = []
    woundNeighborOuterPairsByFrame = []
    for i in range(len(neighborPairsByFrame)):
        nPairs = np.array(neighborPairsByFrame[i])
        wh=list(np.where(nPairs==wv))
        woundNeighbors = nPairs[(wh[0],(wh[1]+1)%2)].tolist() # Get the wound value's friend...
        woundNeighbors2 = []
        for n in nPairs:
            if (n[0] in woundNeighbors) and (n[1] in woundNeighbors):
                woundNeighbors2.append(tuple(sorted(n)))
        secondNeighbors = [] # Anything touching a wound neighbor
        woundNeighborOuterPairs = [] # This essentially constitutes everything that would make up the perimeter around the 1st ring of wound neighbors
        for n in nPairs:
            if (n[0] in woundNeighbors) and not (n[1] in (woundNeighbors+[wv])):
                secondNeighbors.append(n[1])
                woundNeighborOuterPairs.append(tuple(sorted(n)))
            elif (n[1] in woundNeighbors) and not (n[0] in (woundNeighbors+[wv])):
                secondNeighbors.append(n[0])
                woundNeighborOuterPairs.append(tuple(sorted(n)))
        secondNeighbors = sorted(list(set(secondNeighbors)))
        
        woundNeighborsByFrame.append(woundNeighbors)
        woundNeighbors2ByFrame.append(woundNeighbors2)
        secondNeighborsByFrame.append(secondNeighbors)
        woundNeighborOuterPairsByFrame.append(woundNeighborOuterPairs)
    return woundNeighborsByFrame,woundNeighbors2ByFrame,secondNeighborsByFrame,woundNeighborOuterPairsByFrame

def GetContourValuesLengthsAndSubContoursByFrame(watershed,allValsByFrame):
    cVLSByFrame = []
    for frame in range(len(watershed)):
        cVLS = [] # for each subcountour: [value, length, subcontour points]
        for v in allValsByFrame[frame]:
            boundingRect=ImageContour.GetBoundingRect(watershed[frame],v)
            # No longer needed: #contour,turns,vals = ImageContour.GetContour(watershed[0],v,boundingRect=boundingRect,byNeighbor=True)
            perimeterVals,perimeterList,subContours = ImageContour.GetPerimeterByNeighborVal(watershed[frame],v,boundingRect=boundingRect,getSubContours=True)
            subContoursAdj = [(np.array(sc)+[boundingRect[0][0],boundingRect[1][0]]).tolist() for sc in subContours] # Will need to - 0.5 to line up on an overlay
            if len(perimeterList)>0:
                cVLS += [ [sorted([v,perimeterVals[i]]),perimeterList[i],subContoursAdj[i]] for i in range(len(perimeterVals))]
        cVLS.sort(cmpGen(lambda x: x[0]))
        for i in range(len(cVLS)-1,0,-1):
            if cVLS[i-1][0]==cVLS[i][0]:                    # if 2 subcoutours are the same,
                cVLS[i-1][1] = min(cVLS[i-1][1],cVLS[i][1]) # keep only the one with the minimum length computation
                del(cVLS[i])
        cVLSByFrame.append(cVLS)
    return cVLSByFrame

def ContourPlotFromImage(im,neighborPairs):
    _=imshow_vr(im,interpolation='nearest',cmap=cm.gray)
    for i,nPair in enumerate(neighborPairs):
        whX = np.where(  ((im[:-1,:]==nPair[0]) & (im[1:,:]==nPair[1])) |
                         ((im[:-1,:]==nPair[1]) & (im[1:,:]==nPair[0]))  )
        whY = np.where(  ((im[:,:-1]==nPair[0]) & (im[:,1:]==nPair[1])) |
                         ((im[:,:-1]==nPair[1]) & (im[:,1:]==nPair[0]))  )
        for j in range(len(whX[0])):
            x,y = whX[1][j]-0.5 , whX[0][j]+0.5
            _=plt.plot([x,x+1],[y,y],colors[i],linewidth=2)
        for j in range(len(whY[0])):
            x,y = whY[1][j]+0.5 , whY[0][j]-0.5
            _=plt.plot([x,x],[y,y+1],colors[i],linewidth=2)

def CountourPlotFromCVLS(cVLS):
    for cvls in cVLSByFrame[0]:
        _=dplot(cvls[2])

def GetPairLengthsAndBadSpots(pairs,pairsByFrame,cVLSByFrame,tAx):
    pairLengths = []
    badSpots = []
    for pair in pairs:
        pairLengths.append([])
        badSpots.append([])
        for i,cVLS in enumerate(cVLSByFrame):
            broken = False
            for cvls in cVLS:
                if cvls[0]==list(pair):
                    pairLengths[-1].append(cvls[1])
                    broken=True
                    break
            if not broken:
                pairLengths[-1].append(0)
            elif pair not in pairsByFrame[i]:
                badSpots[-1].append([tAx[i],pairLengths[-1][-1]])
    return pairLengths,badSpots
