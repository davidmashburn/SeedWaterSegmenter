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
import FilenameSort as FS
import ExcelHelper
import ImageContour
from ValueReceived import imshow_vr
import GifTiffLoader as GTL
import SeedWaterSegmenter as SWS

from np_utils import flatten,limitInteriorPoints,totuple,removeDuplicates,deletecases,floatIntStringOrNone

#from MyPath import cprint
#from pylab import *

def norm(a,framesExcludeAfter=None):
    if framesExcludeAfter==None:
        return a/np.mean(a)
    else:
        return a/np.mean(a[:framesExcludeAfter])

def openCSV(filename,colDelim=',',rowDelim='\n',allowNone=True):
    '''A very basic .csv reader'''
    fid=open(filename,'r')
    r=fid.read()
    fid.close()
    
    return [[floatIntStringOrNone(i) for i in j.split(colDelim)] for j in r.replace('\r','').replace(' ','').split(rowDelim)]

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
        #cprint(str(i))    
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

def polyArea(vertexList):
    '''Gotta love (almost) one-liners! (and Gauss' Law...)'''
    a=np.array(vertexList,dtype=np.float)
    return abs(np.sum(np.cross(a,np.roll(a,-1,axis=0)))/2)

def CheckForMalformedRegions(watershed):
    '''This goes through all frames and values and uses shapely to test if regions are disjoint in any way'''
    
    ### THIS FUNCTION IMPORTS SHAPELY DIRECTLY ... DO NOT LINK TO IT!
    ### NO NEED TO ADD ANOTHER DEPENDENCY FOR EVERYONE ELSE!
    #import shapely.geometry
    # Actually, shapely is not needed for this operation after all (now uses polyArea)
    
    for frame in range(len(watershed)):
        if watershed[frame]!=None:
            allVals = np.unique(watershed[frame])
            for v in allVals:
                ic = ImageContour.GetContourInPlace(watershed[frame],v)
                #p=shapely.geometry.asPolygon(ic.tolist())
                if (watershed[frame]==v).sum() != polyArea(ic): #p.area:
                     print 'frame:',frame,'value:',v,'-- Value has disjoint regions!'
                #if not p.is_valid: # (temporarily?) removed
                #    print 'frame:',frame,'value:',v,'--Value has regions connected by point only!'


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

def GetContourValuesLengthsAndSubContoursAndOrderOfSubContoursByFrame(watershed,allValsByFrame):
    '''Identify every contour in the watershed segmentation for each frame as defined by each pair of cell IDs (vals).
       Also, return the order of contours that is needed to reconstruct the border surrounding each cell.'''
    cVLSByFrame = []
    orderOfSCsByValueByFrame = [] # List of dicts
    for frame in range(len(watershed)):
        identifier=0 # unique id for each cVLS
        cVLS = [] # for each subcountour: [value, length, subcontour points]
        orderOfSCsByValue = {} # For each cellID, an ordered list of indexes to the subcontours (into cVLS) that reconstruct the contour (negative values mean the sc needs to be flipped)
        for v in allValsByFrame[frame]:
            boundingRect=ImageContour.GetBoundingRect(watershed[frame],v)
            perimeterVals,perimeterList,subContours = ImageContour.GetPerimeterByNeighborVal(watershed[frame],v,boundingRect=boundingRect,getSubContours=True)
            numSCs=len(perimeterVals)
            subContoursAdj = [(np.array(sc)+[boundingRect[0][0],boundingRect[1][0]]).tolist() for sc in subContours] # Will need to - 0.5 to line up on an overlay
            if len(perimeterList)>0:
                orderOfSCsByValue[v] = []
                for i in range(numSCs):
                    newSC = [sorted([v,perimeterVals[i]]),perimeterList[i],subContoursAdj[i],identifier]
                    matchingCVLS = [ sc for sc in cVLS if sc[0]==newSC[0] ]
                    matchingCVLS = [ sc for sc in matchingCVLS if sorted([newSC[2][0],newSC[2][-1]]) == sorted([sc[2][0],sc[2][-1]]) ]  # Should only possibly find 1 match...
                    if matchingCVLS==[]: # New contour!
                        cVLS.append(newSC)
                        orderOfSCsByValue[v].append(identifier)
                        identifier+=1
                    else:
                        matchingCVLS[0][1] = min(matchingCVLS[0][1],newSC[1]) # keep the minimum perimeter length...
                        orderOfSCsByValue[v].append(-matchingCVLS[0][3]) # Minus means it is backwards!        
        cVLS.sort()
        IDs = [i[3] for i in cVLS]
        cVLS = [i[:3] for i in cVLS] # Remove identifiers from cVLS
        # Reindex after sorting...
        for v in allValsByFrame[frame]:       #cmp(x,0) -> sign(x) in python
            orderOfSCsByValue[v] = [cmp(ID,0)*IDs.index(abs(ID)) for ID in orderOfSCsByValue[v]]
        cVLSByFrame.append(cVLS)
        orderOfSCsByValueByFrame.append(orderOfSCsByValue)
    return cVLSByFrame,orderOfSCsByValueByFrame

        ## NOT NEEDED! KEEPING FOR REFERENCE!
        #for i in range(len(cVLS)-1,0,-1):
        #    for j in range(i-1,-1,-1): # Loop backwards through the sorted list of cvls's... if the value pair matches, check the endpoints (they will always be reversed for adjacent regions (always go ccw...))
        #        if cVLS[i][0]!=cVLS[j][0]: # once we no longer match the value pair, we know there are no more matches in the list...
        #            break
        #        ######## VERIFY THIS ACTUALLY WORKS THE SAME WAY!!!
        #        elif (cVLS[i][2][-1],cVLS[i][2][0]]) == (cVLS[j][2][0],cVLS[j][2][-1]):  # if 2 subcoutours are the same,
        #            if cVLS[j][1]>cVLS[i][1]:
        #                cVLS[j],cVLS[i] = cVLS[i],cVLS[j] #swap!
        #            shortest = min(cVLS[j][1],cVLS[i][1])                              # keep only the one with the minimum length computation
        #            
        #            cVLS[j][1] = shortest
        #            del(cVLS[i])
        #            break

def GetXYListAndPolyListFromCVLS(cVLS,allValsByFrame,orderOfSCsByValueByFrame):
    numFrames = len(cVLS)
    xyList = [ sorted(list(set(  tuple(pt) for c in cvlsByVal for pt in c[2]  ))) for cvlsByVal in cVLS ]
    polyList = []
    
    for t in range(numFrames):
        polyList.append({})
        for v in allValsByFrame[t]:
            subContours = [ ( cVLS[t][-index][2][::-1] if index<0 else cVLS[t][index][2] )
                           for index in orderOfSCsByValueByFrame[t][v] ] # reconstruct the sc's, flipping if index is negative
            polyList[-1][v] = [ xyList[t].index(totuple(pt)) for sc in subContours for pt in sc ]
            polyList[-1][v] = removeDuplicates(polyList[-1][v])+[polyList[-1][v][0]] # Remove interior duplication...
    return xyList,polyList

def GetCVLSWithLimitedPointsBetweenNodes(cVLS,allValsByFrame,splitLength=1,fixedNumInteriorPoints=None):
    allValues = list(set(tuple(flatten(allValsByFrame))))
    allPairs = sorted(list(set([tuple(c[0]) for cVLSByFrame in cVLS for c in cVLSByFrame])))

    if fixedNumInteriorPoints:
        numInteriorPoints = {p:fixedNumInteriorPoints for p in allPairs}
    else:
        minLength = {}
        for p in allPairs:
            minLength[p] = min([ len(c[2]) for cVLSByFrame in cVLS for c in cVLSByFrame if tuple(c[0])==p ])
        numInteriorPoints = {}
        for p in allPairs:
            numInteriorPoints[p] = minLength[p]//splitLength
    cVLS2 = deepcopy(cVLS)
    for cvlsByValue in cVLS2:
        for c in cvlsByValue:
            c[2] = limitInteriorPoints(c[2],numInteriorPoints[tuple(c[0])])
    return cVLS2

def GetXYListAndPolyListWithLimitedPointsBetweenNodes(cVLS,allValsByFrame,orderOfSCsByValueByFrame,splitLength=1,fixedNumInteriorPoints=None):
    cVLS2 = GetCVLSWithLimitedPointsBetweenNodes(cVLS,allValsByFrame,splitLength,fixedNumInteriorPoints)
    return GetXYListAndPolyListFromCVLS(cVLS2,allValsByFrame,orderOfSCsByValueByFrame)

def ContourPlotFromImage(im,neighborPairs):
    _=imshow_vr(im,interpolation='nearest',cmap=plt.cm.gray)
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

def ContourPlotFromCVLS(cVLSByFrame,frame=0):
    for cvls in cVLSByFrame[frame]:
        cvls=np.array(cvls[2])
        _=plt.plot( cvls[:,0], cVLS[:,1] )

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

# If you need to use the cVLS's anyway, why compute them twice??
# BTW, this works fine for static but DOES NOT give good connection from frame to frame
#   1: No connection from point to point through time
#   2: In-between points can come and go willy-nilly!!
def MakePolygonNetworkFromWaterSeg(waterArr,minSplit=30,allValsByFrame=None,cVLS=None):
    if allValsByFrame==None:
        allValsByFrame = [ list(set(w.flat)) for w in waterArr ]
    if cVLS==None:
        cVLS = SWHelpers.GetContourValuesLengthsAndSubContoursByFrame(waterArr,allValsByFrame)
    allPtsList = []
    allSegsList = []
    for t in range(len(cVLS)):
        contours = [i[2] for i in cVLS[t]]
        nc = len(contours)
        lc = [len(i) for i in contours]
        div = [i//30 for i in lc]
        sep = [ (lc[i] if div[i]==0 else lc[i]//div[i])  for i in range(nc)]
        #subPtsOld = [ ([] if div[i]==0 else contours[i][ ::(len(contours[i])//div[i]) ]  )
        #          for i in range(nc) ]
        inds = [  ( [0,lc[i]-1] if (lc[i]-1 < 2*minSplit) else
                    range(0,lc[i]-2,sep[i])+[lc[i]-1] )
                for i in range(nc)  ]
        
        subPts = [[ tuple(contours[i][j]) for j in inds[i] ] 
                                          for i in range(nc) ]
        allPts = sorted(list(set(flatten( subPts ))))
        subPtsInds = [[ allPts.index(j) for j in i ] 
                                        for i in subPts ]
        allSegsInds = flatten([ zip(i[:-1],i[1:]) for i in subPtsInds ])
        allPtsList.append(allPts)
        allSegsList.append(allSegsInds)
        
    return allPtsList,allSegsList

def MultNumOrNone(x,y):
    if x==None or y==None:
        return None # None's -> None's
    else:
        return x*y

def SquareNumOrNone(x):
    if x==None:
        return None
    else:
        return x**2

def GetNormalizedArea(areaFromCSV):
    NormArea = []
    for i in range(areaFromCSV.shape[1]):
        if None in areaFromCSV[:,i].tolist():
            NormArea.append(areaFromCSV[:,i])
        else:
            NormArea.append( areaFromCSV[:,i] / np.mean(areaFromCSV[:,i]) )
    NormArea=np.array(NormArea).T
    return NormArea

def GetWoundValues(d):
    if 'ManualInputs.py' in os.listdir(d):
        exec( open(d+'/ManualInputs.py','r').read().replace('\r','').split('\n')[0] )
        return woundVals
    else:
        print 'ManualInputs.py does not exist!'
def GetNeighborValues(d,wv,ind):
    if 'Neighbors.py' in os.listdir(d):
        exec( open(d+'/Neighbors.py','r').read() )
        allVals = []
        for v in wv:
            n=neighbors[ind][v-2]
            if n!=None:  allVals+=n
        allVals = list(set(allVals).difference(wv))
        return allVals
    else:
        print 'Neighbors.py does not exist!'

#####################################################################
# Integrated intensity functions:
#####################################################################
def OutlineMasksPlot(arr,waterArr,cellID,t=0):
    '''Plot the image and 3 outline masks to demonstrate the effects of different numbers of dilate/erodes'''
    wv = (waterArr[t]==cellID)
    bright = np.zeros_like(arr[t])
    
    wh=np.where(wv)
    bright[wh] = arr[t][wh]
    
    plt.clf()
    plt.subplot(231)
    plt.imshow(arr[t],cmap=plt.cm.gray)
    
    plt.subplot(232)
    plt.imshow(bright,cmap=plt.cm.gray)
    
    for i in range(1,5):
        wo = ndimage.binary_dilation(wv,iterations=i) -  \
             ndimage.binary_erosion(wv,iterations=i)
        wh=np.where(wo)
        bright = np.zeros_like(arr[t])
        bright[wh]=arr[t][wh]
        
        plt.subplot(2,3,i+2)
        plt.imshow(bright,cmap=plt.cm.gray)

def DilateMinusErode(binaryArr,structure=None,iterations=1):
    if iterations<1:
        print 'You must specify at least one iteration!!!'
        return
    else:
        return ndimage.binary_dilation(binaryArr,structure=structure,iterations=iterations) -  \
                ndimage.binary_erosion(binaryArr,structure=structure,iterations=iterations)

def GetCellVolume(waterArr,cellID,structure=None,iterations=0):
    sums = []
    for t in range(len(waterArr)):
        w = (waterArr[t]==cellID)
        if iterations>0: #if iterations is >0, w turns into the outline region instead
            w = DilateMinusErode(w,structure=structure,iterations=iterations)
        sums.append( w.sum() )
    return np.array(sums)

def GetTotalIntegratedBrightness(arr,waterArr,cellID,structure=None,iterations=0):
    '''Sum the brightness of the original image under the final segmented region (2D+time or 3D+time)'''
    sums = []
    bright = np.zeros_like(arr[0])
    for t in range(len(arr)):
        bright[:]=0
        w = (waterArr[t]==cellID)
        if iterations>0: #if iterations is >0, w turns into the outline region instead
            w = DilateMinusErode(w,structure=structure,iterations=iterations)
        wh = np.where(w)
        bright[wh]=arr[t][wh]
        sums.append( bright.sum() )
    return np.array(sums)

def GetEdgeIntegratedBrightness(arr,waterArr,cellID,thickness=2,skip3DTopAndBottom=False):
    '''Sum the brightness of the original image under the edges of the final segmented region, with a specificed thickness (even integers only) (2D+time or 3D+time)'''
    if thickness%2==1 or thickness<2:
        print 'thickness parameter must be an even integer > 0!'
        return
    iterations = thickness//2
    if not skip3DTopAndBottom:
        return GetTotalIntegratedBrightness(arr,waterArr,cellID,iterations=iterations)
    elif arr.ndim==4:
        structure = np.array([ [[0,1,0],[1,1,1],[0,1,0]] ],np.bool)
        return GetTotalIntegratedBrightness(arr,waterArr,cellID,structure=structure,iterations=iterations)
    else:
        'You can only use skip3DTopAndBottom=True on a 4D array!'
        return

def GetVolumesAndIntegratedBrightnesses(testArrays,waterArr,cellID,thickness,skip3DBottomAndTop=True):
    volume = GetCellVolume(waterArr,cellID)
    volumeEdge = GetCellVolume(waterArr,cellID,structure=[[[0,1,0],[1,1,1],[0,1,0]]],iterations=thickness//2)
    
    intensities = [ GetTotalIntegratedBrightness(i,waterArr,cellID=cellID)
                   for i in testArrays ]
    edgeIntensities = [ GetEdgeIntegratedBrightness(i,waterArr,cellID,thickness=thickness,skip3DTopAndBottom=skip3DBottomAndTop)
                       for i in testArrays ]
    
    intensitiesPerPixel = [ i/volume for i in intensities ]
    edgeIntensitiesPerPixel = [ i/volumeEdge for i in edgeIntensities ]
    
    return [volume,volumeEdge],intensities,edgeIntensities,intensitiesPerPixel,edgeIntensitiesPerPixel

#####################################################################


#####################################################################
# An old function from CellFunctions.py for plotting some weird gradient angles...
# Just put it here to keep from losing the code
#####################################################################
def GradientAnglePlot(f='/media/sda1/Documents and Settings/Owner/My Documents/VIIBRE--ScarHealing/ActiveData/Stack_Zproject_GBR_DC.gif'):
    r=VolumeGif.LoadMonolithicGif(f)[0]
    
    clf()
    numSub=5
    g=4
    subplot(1,numSub,1)
    test=ndimage.gaussian_filter(1.*r[120:180,120:180],g)
    pdx=ndimage.convolve(test,n.array([[-1],[0],[1]]))
    pdy=ndimage.convolve(test,n.array([[-1,0,1]]))
    imshow(test)
    subplot(1,numSub,2)
    imshow(pdx)
    subplot(1,numSub,3)
    imshow(pdy)
    subplot(1,numSub,4)
    imshow(n.arctan(pdy/pdx))
    subplot(1,numSub,5)
    imshow(ndimage.convolve(test,n.array([[-1,0,1],[0,0,0],[1,0,-1]])))

