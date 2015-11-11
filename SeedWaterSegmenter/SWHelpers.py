#!/usr/bin/env python
''''A useful collection of helper functions for SWS'''
import os
import imp

import numpy as np
from scipy import ndimage
import matplotlib.pyplot as plt

import FilenameSort as FS
import ImageContour
import GifTiffLoader as GTL
from GifTiffLoader import Image # either pillow or PIL
import mahotas

#from list_utils:
from np_utils import deletecases,floatIntStringOrNone,ziptranspose,shape_multiply

#from MyPath import cprint
#from pylab import *

def norm(a,framesExcludeAfter=None):
    '''Standard norm (divide by the mean). Can ignore frames at the end when computing the mean'''
    if framesExcludeAfter==None:
        return a/np.mean(a)
    else:
        return a/np.mean(a[:framesExcludeAfter])

def openCSV(filename,colDelim=',',rowDelim='\n',allowNone=True):
    '''A very basic .csv reader'''
    fid=open(filename,'r')
    r=fid.read()
    fid.close()
    
    return [[floatIntStringOrNone(i)
              for i in j.split(colDelim)]
             for j in r.replace('\r','').replace(' ','').split(rowDelim)]

def GetWoutline(watershed,dilate):
    woutline = np.zeros(watershed.shape,np.int)
    woutline[:-1,:-1] = (np.diff(watershed,axis=0)[:,:-1]!=0) + (np.diff(watershed,axis=1)[:-1]!=0)
    
    if dilate>0:
        woutline = ndimage.binary_dilation(woutline,iterations=dilate)
    
    return woutline

def GetSeedArray(watershed,SeedsPyFile,frameNumber):  # THIS FUNCTION IS CURRENTLY BROKEN
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
        pass
        ################## This function no longer exists:
        ##################SWS.PointsToArray(seedList[frameNumber], seedVals[frameNumber],seedArray)
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
def GetTimeAxisAdj(SegmentsDirectory):
    sList = []
    timeAdj = 0
    count=0
    for i in FS.getSortedListOfFiles(SegmentsDirectory,'*.tif'):
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
        sList.append(s)
    return sList

def SaveOverlayImages(tifStackData,SeedsPyFile,SegmentsDirectory,red=True):
    # This function is more-or-less just a template right now, has MAJOR issues...
    r = GTL.LoadFileSequence(SegmentsDirectory, wildcard='*.tif')
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
    '''Gotta love (almost) one-liners! (and Green's Theorem...)'''
    a=np.array(vertexList,dtype=np.float)
    return abs(np.sum(np.cross(a,np.roll(a,-1,axis=0)))/2)

def ConvertOutlinesToWatershed(origArr,outlines,useDilate=False,structure=[[0,1,0],[1,1,1],[0,1,0]]):
    if useDilate:
        outlines=ndimage.morphology.binary_dilation(origArr,structure=structure)
    labels = ndimage.label(1-outlines)[0]
    wh=np.where(labels==0)
    labels[wh]=-1
    labels+=1
    labels = labels.astype(np.uint16)
    labels[0]=1
    labels[-1]=1
    labels[:,0]=1
    labels[:,-1]=1
    
    water = mahotas.cwatershed(origArr,labels)
    return water


def CheckForMalformedRegions(watershed,usePrint=True):
    '''This goes through all frames and values and uses a poylogon area
       algorithm to test if regions are disjoint in any way.'''
    
    outputString = ''
    
    for frame in range(len(watershed)):
        print 'Checking Frame',frame
        if watershed[frame]!=None:
            allVals = np.unique(watershed[frame])
            if len(allVals)<2:
                continue # no point in checking if it's empty or just a big solid rectangle...
            for v in allVals:
                if v==1:
                    # For the background, take the inverse instead
                    filledRegion = (watershed[frame]!=1).astype(np.int)
                else:
                    filledRegion = (watershed[frame]==v).astype(np.int)
                ic = ImageContour.GetContourInPlace(filledRegion,1)
                if filledRegion.sum() != polyArea(ic):
                    s = ' '.join(['frame:',str(frame),', value:',str(v)])
                    outputString = outputString + ' ' + s + '\n'
    
    if outputString=='':
        outputString = 'No malformed regions found!'
    
    if usePrint:
        print outputString
    else:
        return outputString
        

def MakeCellNetworkJsonFiles(waterArr,d):
    print 'Making static (simple networks) json file:'
    ImageContour.SubContourTools.GetCellNetworkListStatic(waterArr,d,forceRemake=True)
    print 'Made static json file!'
    
    print 'Making matched networks json file:'
    _,_,allMatched = ImageContour.SubContourTools.GetMatchedCellNetworkListsPrevNext(waterArr,d,forceRemake=True)
    print 'Made matched networks json file!'
    
    return allMatched

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
        woundNeighborOuterPairs = [] # Everything making up the perimeter around 1st ring of wound neighbors
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
            if n!=None:
                allVals+=n
        allVals = list(set(allVals).difference(wv))
        return allVals
    else:
        print 'Neighbors.py does not exist!'

#####################################################################
# Function for converting an outline array to a watershed array:
#####################################################################
def WatershedFillDoubleSizedOutlineArray(arr):
    '''Take an image array similar to the "Outlines"
       (but with lines=0, cells=1 instead)
       and turn it into a (double-sized) fully segmented image'''
    arrEdit = ndimage.measurements.label(arr)[0]
    arrDouble = shape_multiply(arrEdit,[2,2]).astype(np.uint16)
    arrInvDouble = (1-shape_multiply(arr,[2,2])).astype(np.uint16)
    waterArr = mahotas.cwatershed(arrInvDouble,arrDouble)
    return waterArr

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
        wo = ( ndimage.binary_dilation(wv,iterations=i) -
               ndimage.binary_erosion(wv,iterations=i)  )
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
        return ( ndimage.binary_dilation(binaryArr,structure=structure,iterations=iterations) -
                 ndimage.binary_erosion(binaryArr,structure=structure,iterations=iterations)  )

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
    '''Sum the brightness of the original image under the edges of the final segmented region,
       with a specificed thickness (even integers only) (2D+time or 3D+time)'''
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
        print 'You can only use skip3DTopAndBottom=True on a 4D array!'
        return

def GetVolumesAndIntegratedBrightnesses(testArrays,waterArr,cellID,thickness,skip3DBottomAndTop=True):
    '''Takes testArrays (one or more sets of raw data), waterArr (segmentation data), and cellID to make measurements of a cell
       Returns the volume under the cell and the volume of the xy-rim around the cell with a certain thicknes
           and for each of these for each testArray: the integrated intensity and average intensity per pixel of the raw array under these same regions
       again, testArrays is a list of arrays and the return values: "intensities","edgeIntensities","intensitiesPerPixel","edgeIntensitiesPerPixel"
           all return lists of measurements for each testArray
       Uses GetCellVolume in both direct mode (iterations=0) and edge mode (runs dilate-erode <iterations> times) to get a band around the edge'''
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
'''
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
'''
