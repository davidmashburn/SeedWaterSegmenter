#!/usr/bin/env python
"""SeedWater Segmenter is a watershed segmentation / registration
program for image stacks developed at Vanderbilt University using
Mahotas, WxPython, Matplotlib, Numpy, Scipy, and PIL.
SeedWater enables users to generate seeds automatically and then edit
them by hand in Matplotlib using custom interactions.  Seeds for
subsequent frames are generated as centroids of regions in the previous
frame.

Also, GifTiffLoader is a wrapper to automatically load Tiff and Gif
files as numpy arrays using PIL.
GifTiffLoader also relies on FilenameSort and cmpGen.

Sponsored by the NSF and HFSP through the Shane Hutson Laboratory, part of
Vanderbilt Institute of Integrative Biosystems Research (VIIBRE)."""

__author__ = "David N. Mashburn <david.n.mashburn@gmail.com>"

import os,sys
from copy import copy,deepcopy
from time import time, sleep
import imp

import numpy as np
from numpy import nonzero
from numpy.random import rand

import scipy.ndimage
import scipy.sparse
import Image
import GifTiffLoader as GTL
from cmpGen import cmpGen

import wx
import math
from scipy.ndimage import center_of_mass,gaussian_filter,median_filter

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.nxutils import points_inside_poly
from mpl_polygon_lasso import PolyLasso

import mahotas

import ImageContour
import EllipseFitter
import ExcelHelper

class Timer:
    def __init__(self):
        self.timeStart=time()
    def test(self):
        print time()-self.timeStart
        self.timeStart=time()
    def reset(self):
        self.timeStart=time()
timer = Timer()

try:
    import cv
    HAS_CV=True
    def cv2array(im):
        depth2dtype = {
                       cv.IPL_DEPTH_8U: 'uint8',
                       cv.IPL_DEPTH_8S: 'int8',
                       cv.IPL_DEPTH_16U: 'uint16',
                       cv.IPL_DEPTH_16S: 'int16',
                       cv.IPL_DEPTH_32S: 'int32',
                       cv.IPL_DEPTH_32F: 'float32',
                       cv.IPL_DEPTH_64F: 'float64',
                      }
        arrdtype=im.depth
        a = np.fromstring(
            im.tostring(),
            dtype=depth2dtype[im.depth],
            count=im.width*im.height*im.nChannels)
        a.shape = (im.height,im.width,im.nChannels)
        return a

    def array2cv(a):
        dtype2depth = {
                       'uint8':   cv.IPL_DEPTH_8U,
                       'int8':    cv.IPL_DEPTH_8S,
                       'uint16':  cv.IPL_DEPTH_16U,
                       'int16':   cv.IPL_DEPTH_16S,
                       'int32':   cv.IPL_DEPTH_32S,
                       'float32': cv.IPL_DEPTH_32F,
                       'float64': cv.IPL_DEPTH_64F,
                      }
        try:
            nChannels = a.shape[2]
        except:
            nChannels = 1
        cv_im = cv.CreateImageHeader((a.shape[1],a.shape[0]),
                                     dtype2depth[str(a.dtype)],
                                     nChannels)
        cv.SetData(cv_im, a.tostring(),
                   a.dtype.itemsize*nChannels*a.shape[1])
        return cv_im

    def cvWater(a,m):
        im = cv.fromarray(np.array(a,dtype=np.uint8))
        imC = cv.CreateImage(cv.GetSize(im), 8, 3) # Color
        cv.CvtColor(im, imC, cv.CV_GRAY2BGR)
        markers = array2cv(np.array(m,dtype=np.int32))
        #cv.CreateImage(cv.GetSize(im), cv.IPL_DEPTH_32S, 1)
        
        cv.Watershed(imC, markers)
        w=cv2array(markers)[:,:,0]
        
        return w
except ImportError:
    HAS_CV=False

letterKeys = '''a: Change to Add/Delete Mode
b: Background Subtraction
c: Change to Center Mode (select wound center for post-processing)
d: Change to Draw Mode
e: Change to Extra Seed Mode
f: Force Remake of Seeds From Centroids of Previous Frame
   (Shift-F directly copies the seeds from the previous frame)
g: Gaussian Filter
h: Sharpen Filter
i: Invert Image
j: Join selected seeds (after lasso)
k: Median Filter
l: Change to lasso mode
m: Change to Move Mode (Allows pan/zoom from toolbar)
n: Next Frame
o: Open Images For Segmentation
p: Previous Frame (Shift-P goes to the last visited frame)
q: Quit
r: Reset all image processing (notably Gaussian)
s: Save
t: Toggle overlays (useful when outline is blocking image)
u: Undo Seed Action
v: Change the value of a region (highlighted in e or x mode)
w: Run watershed (automatic now) - Refreshes Map Plot also...
x: Change to Switch Regions Mode
y: Special "Do Nothing" watershed...  For use on bad frames
   so frame after after will retain seeds from frame before.'''

# Options for specific users... let me know if you want a
# specific set of parameters, too...
username = os.path.split(os.path.expanduser('~'))[-1]

if username in ['mashbudn']:
    DONT_PANIC=False
    DEFAULT_SEED_SIZE=2
    USE_DEBUG_PRINT=True
    STEAL_B_KEY = False
elif username in ['Holley']:
    DONT_PANIC=True
    DEFAULT_SEED_SIZE=2
    USE_DEBUG_PRINT=False
    STEAL_B_KEY = False
elif username in ['Aroshan']:
    DONT_PANIC=False
    DEFAULT_SEED_SIZE=3
    USE_DEBUG_PRINT=False
    STEAL_B_KEY = True
    
else:
    DONT_PANIC=False
    DEFAULT_SEED_SIZE=3
    USE_DEBUG_PRINT=False
    STEAL_B_KEY = False

def dprint(a):
    """Debug print"""
    if USE_DEBUG_PRINT:
        print a

txtHeader=["# SeedWater notes file version 1",
           "# This file is the generated notes from SeedWater.",
           "# If you edit it by hand, do not change the frame headers!"]
txtsepA="# ---------------- Frame "
txtsepB=" ----------------"

def deletecases(l,valList):
    if not hasattr(valList,'__iter__'):
        valList=[valList]
    l=copy(l) # work on a copy
    for v in valList:
        for i in range(l.count(v)):
            l.remove(v)
    return l

def GetPointsAtRadius(arr_shape,x,y,r):
    """Create a binary image of a circle radius r"""
    w,h=arr_shape
    r2 = int(r**2)
    l=[]
    for i in range(int(x-r),int(x+r+1)):
        for j in range(int(y-r),int(y+r+1)):
            if (i-x)**2+(j-y)**2 < r2:
                if 0<=i<w and 0<=j<h:
                    l+=[[i,j]]
    return l

def ImageCircle(r):
    """Create a binary image of a circle radius r"""
    im = np.zeros([2*r-1,2*r-1],dtype=np.int)
    r2 = r**2
    for i in range(2*r-1):
        for j in range(2*r-1):
            if (i-r+1)**2+(j-r+1)**2 < r2:
                im[i,j] = 1
            else:
                im[i,j]=0
    return im

def blitCircleToArray(arr,x,y,r,val):
    xL,xH = np.clip(x-r+1,0,arr.shape[0]), np.clip(x+r,0,arr.shape[0])
    yL,yH = np.clip(y-r+1,0,arr.shape[1]), np.clip(y+r,0,arr.shape[1])
    xcL,xcH = xL-(x-r+1), 2*r-1 + xH-(x+r)
    ycL,ycH = yL-(y-r+1), 2*r-1 + yH-(y+r)
    #print (xL,xH,yL,yH),(xcL,xcH,ycL,ycH)
    #print xH-xL,yH-yL,xcH-xcL,ycH-ycL
    c=ImageCircle(r)[xcL:xcH,ycL:ycH]
    #print arr[xL:xH,yL:yH].shape,c.shape
    arr[xL:xH,yL:yH] *= (1-c)
    arr[xL:xH,yL:yH] += (val*c)

def MaxMinFinder(x,useMin=True):
    sl  = [slice(None,-1),slice(None),slice(1,None)]
    sl2 = [slice(1,None),slice(None),slice(None,-1)]
    # Create 9 arrays offset +-1 in each directions (neighbors)
    xNeigh = np.zeros([9]+list(x.shape))
    for i in range(9):
        xNeigh[i,sl[i//3],sl[i%3]] = x[sl2[i//3],sl2[i%3]]
    xNeigh[4]=0 # set self to 0 -- 4 because 3*1+1=4 (1,1)
    
    maxMap = (x*0+1).astype(np.bool)
    for i in range(9):
        if i!=4:
            if useMin:
                maxMap = maxMap * (x<xNeigh[i])
            else:
                maxMap = maxMap * (x>xNeigh[i])
    return maxMap

# Thanks Wikipedia!!
# A line-drawing algorithm by the pixel...
def BresenhamFunction(p0,p1):
     [x0,y0], [x1,y1] = p0, p1
     steep = abs(y1-y0) > abs(x1-x0)
     if steep:
         x0,y0 = y0,x0 # swap
         x1,y1 = y1,x1 # swap
     if x0 > x1:
         x0,x1 = x1,x0 # swap
         y0,y1 = y1,y0 # swap
     deltax = x1-x0
     deltay = abs(y1-y0)
     error = deltax/2
     y = y0
     if y0 < y1:
         ystep=1
     else:
         ystep=-1
     l=[]
     for x in range(x0,x1+1):
         if steep:
             l.append([y,x])
         else:
             l.append([x,y])
         error = error-deltay
         if error < 0:
             y = y+ystep
             error = error+deltax
     return l

def CalcPolyPerimeter(pts):
    d = 0
    for i in range(len(pts)):
        d += np.sqrt( (pts[i][0]-pts[i-1][0])**2 + (pts[i][1]-pts[i-1][1])**2 )
    return d

def GetMinPoly(pts):
    minDist = 1e10
    minPoly=[]
    # Brute Force; hey, at least it works!
    for p in itertools.permutations(pts):
        newDist = CalcPolyPerimeter(p)
        if newDist<minDist:
            minDist=newDist
            minPoly = p
    return minPoly

def dilateByCircle(l,arr_shape,r):
    l_old=copy(l)
    for i in l_old:
        print i
        els = GetPointsAtRadius(arr_shape,i[0],i[1],r)
        for el in els:
            if not el in l:
                l.append(el)

def CreateOutlines(arr,walgorithm='PyMorph'):
    if walgorithm in ['OpenCV']: # This is now deprecated and will not produce correct results on saving!
        return (arr==-1)
    elif walgorithm in ['PyMorph','Dummy']: # Later add CV2PM
        newArr=np.zeros(arr.shape)
        newArr[:-1,:-1] = (np.diff(arr,axis=0)[:,:-1]!=0) + \
                          (np.diff(arr,axis=1)[:-1]!=0)
        return newArr
    else:
        raise ValueError, walgorithm+' is not a valid watershed algorithm!'

IorN = lambda i: None if i=='None' else int(i)
ForN = lambda f: None if f=='None' else float(f)

def GetReportPixel(sfORwd): # create and return the report_pixel function...
    def report_pixel(x,y):
        if sfORwd.__class__ == SegmenterFrame:
            wd=sfORwd.wd
        elif sfORwd.__class__ == WatershedData:
            wd=sfORwd
        else:
            raise TypeError("sfORwd must be a WatershedData or SegmenterFrame Object!")
        
        s=wd.watershed[wd.index].shape
        if 0<x<s[1] and 0<y<s[0]:
            return "value=" + str(wd.watershed[wd.index][y,x]) + \
                   "  x=" + str(x) + " y=" + str(y)
        else:
            return "x=" + str(x) + " y=" + str(y)
    return report_pixel

def sumN(l):
    """Sum that counts None as 0's"""
    return sum([(0 if i==None else i) for i in l])
def meanN(l):
    """Mean that counts None as 0's"""
    return sumN(l)*1./len(l)
def meanN_skip(l):
    return sumN(l)*1./sum([(0 if i==None else 1) for i in l ])

def CreateTimeAxis(timeIntervals,gapIntervals,numberFramesPerSeries,numberOfFrames):
    if not ( len(timeIntervals)==(len(gapIntervals)+1)==len(numberFramesPerSeries) ):
        print 'intervalList and numberWintervalList must have the same length to be valid!!!'
        print 'gapIntervals must be one element shorter than the others!'
        return
    #Example:
    #timeIntervals = [10,20]
    #gapIntervals=[21] # Should be one less
    #numberFramesPerSeries=[30,60]
    
    if sum(numberFramesPerSeries)<numberOfFrames:
       numberFramesPerSeries[-1] = numberOfFrames - sum(numberFramesPerSeries[:-1])
    
    timeAxis = []
    for i in range(len(timeIntervals)):
        ran = np.arange(0,timeIntervals[i]*numberFramesPerSeries[i],timeIntervals[i])
        if i==0:
            add=0
        else:
            add=timeAxis[-1]+gapIntervals[i-1]
        timeAxis += (ran+add).tolist()
    return timeAxis

def GetPlacementName(num):
    if num==1:   return '1st'
    elif num==2: return '2nd'
    elif num==3: return '3rd'
    else:        return str(num)+'th'

# Trying to use this for both area and major/minor axis causes issues...
# wound components are added and ring components are averaged...
# Just something to keep in mind
def GetBinnedValueForExcel(value,name,bins,timeAxis,woundVals):
    # THE GOOD WAY :)
    value=np.array(value)
    m=[]
    #time,wound,ring1,ring#
    valueBinned = [[["Wound "+name,None,"Time"]]+[[p,None,timeAxis[i]] for i,p in enumerate(map(sumN,value[:, woundVals-np.r_[2] ]))]]
    for b in range(len(bins)):
        if bins[b]!=[]:
            m.append(map(meanN, value[:, bins[b] ] ))
            #[["Test1"],["Test2"],["Test3"]]+
            valueBinned.append(value[:, bins[b] ].tolist())
            for j in range(len(valueBinned[-1])): # loop through time...
                valueBinned[-1][j]+=[None,timeAxis[j],m[-1][j]]
    firstSheet = [ [timeAxis[i],p]+[m[b][i] for b in range(len(bins))] for i,p in enumerate(map(sumN,value[:, woundVals-np.r_[2] ]))]
    valueBinned.insert(0,[["time","wound"]+["ring"+str(i+1) for i in range(len(bins))]]+firstSheet)
    return valueBinned

    ## Crappy way... :(
    #oldIndex = wd.index
    #areaBinned = []
    #for i in range(wd.length):
    #    if wd.seedList[i] in [None,[]]:
    #        break
    #    if i>=len(wd.woundCenters):
    #        print 'No Wound Centers After Frame',i,'!'
    #        break
    #    wd.index=i
    #    
    #    wd.UpdateValuesList()
    #    area=wd.CalculateArea(doUpdate=False)
    #    areaBinned.append([])
    #    for b in bins:
    #        inds=[]
    #        for v in b:
    #            if v in wd.valList:
    #                inds.append(wd.valList.index(v))
    #        areaBinned[-1].append( [area[i] for i in inds] )
    #wd.index=oldIndex

# operates on the neighbor list in place
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

def GetAllValsByFrame(neighborPairsByFrame):
    allValsByFrame = []
    for neighborPairs in neighborPairsByFrame:
        allValsByFrame.append(sorted(list(set([j for i in neighborPairs for j in i]))))
    return allValsByFrame


# Implement these later... a partially-changed version is in Python folder...
class WoundContours:
    WNv=[]
    WNl=[]
    NNv=[]
    NNl=[]
    NOv=[]
    NOl=[]
    def Extend(self):
        for i in [self.WNv,self.WNl,self.NNv,self.NNl,self.NOv,self.NOl]:
            i.append([])
    def Fill(self,wc):
        for i in [self.WNv,self.WNl,self.NNv,self.NNl,self.NOv,self.NOl]:
            if len(i)==0:
                print 'Lists have no elements!'
                return
        self.WNv = wc.WNv
        self.WNl = wc.WNl
        self.NNv = wc.NNv
        self.NNl = wc.NNl
        self.NOv = wc.NOv
        self.NOl = wc.NOl

def SeedListToSparse(seeds,vals,shape):
    if seeds==None:
        return None
    else:
        row,col = np.array(seeds).T
        sparse = scipy.sparse.coo_matrix((vals,(row,col)), shape=shape).tolil()
        return sparse

def GetMapPlotRandomArray():
    np.random.seed(0)
    mapPlotRandomArray=np.array([np.random.random(10000)*236+20,
                                 np.random.random(10000)*236+20,
                                 np.random.random(10000)*236+20], dtype=np.uint8)
    mapPlotRandomArray[:,0]=255
    mapPlotRandomArray[:,1]=255
    return mapPlotRandomArray

#def MakeCmapFromArray(a):
#    r,b,g = [],[],[]
#    for i in range(len(a[0])):
#        p=i*1./(len(a[0])-1)
#        r.append((p,a[0][i],a[0][i]))
#        g.append((p,a[1][i],a[1][i]))
#        b.append((p,a[2][i],a[2][i]))
#    segmentdata = {'red':tuple(r),'green':tuple(g),'blue':tuple(b)}
#    return matplotlib.colors.LinearSegmentedColormap('rand4SW',segmentdata)

# This class now deals with image stacks as well as images...
class WatershedData:
    def __init__(self,arrayIn,previousSeeds=None):
        # Change this to allow orig, min-max filtered, and sharpened data
        # Need a way to mark the "Un-seeds" that exist in the background only...
        if len(arrayIn.shape)==2:
            arrayIn=arrayIn.reshape([1,arrayIn.shape[0],arrayIn.shape[1]])
        self.length=arrayIn.shape[0] # number of frames
        self.shape = arrayIn.shape
        # I AM CHANGING THE DEFAULT BEHAVIOR!!!!!!!!!
        # I don't really think this will change the watersheds, but it
        #  really might change the brightness of the image...
        if arrayIn.dtype==np.uint16:
            self.origData = arrayIn
            self.origData *= int( (2**16-1)/arrayIn.max() )
        elif arrayIn.dtype==np.uint8:
            self.origData = np.array(arrayIn,dtype=np.uint16)
            self.origData *= int( (2**16-1)/arrayIn.max() )
        #    self.origData *= 2**8
        else:
            self.origData = GTL.DivideConvertType(arrayIn,bits=16,maxVal=2**15-1,
                                                  zeroMode='stretch',maxMode='stretch')
        # 1 full copy...
        
        #bg=gaussian_filter(arrayIn,64)
        #self.origData = arrayIn+bg.max()-bg
        #self.origData = gaussian_filter(self.origData,2) + gaussian_filter(self.origData,6).max()-gaussian_filter(self.origData,6)
        #self.origData = self.origData-self.origData.min()
        
        self.gaussSigma=0
        self.filterData = np.array(self.origData)  # 2 full copy
        
        self.sparseList = [None]*self.length # will be a list of sparse matrices (lil format)
        self.seedSelections = [None]*self.length # will also be a sparse matrix list
        self.notes = [None]*self.length
        self.woundCenters = [None]*self.length
        
        self.background = True # Set to True for better results...
        self.watershed=np.zeros(self.shape,dtype=np.int16)  # 3 full copy
        self.seedArray = np.zeros(self.shape[1:],dtype=np.int)
        self.woutline=np.zeros(self.shape[1:],dtype=np.int)
        self.index=0
        self.lastFrameVisited=0
        # walgorithm should be 'PyMorph', Dummy, ['OpenCV', 'CV2PM']
        self.walgorithm = ['PyMorph']*self.length
        self.framesVisited=[False]*self.length
        self.framesVisited[self.index]=True
        self.rgba=None # Will hold the data for the color plot...
        self.rgbaH=None # Will hold the data for the highlight plot...
        self.rgbM=None # Will hold the data for the map plot...
        self.bgPlot=None
        self.overlayVisible=True
        self.colorplot=None
        self.mapPlot=None
        self.mapCentroids=[]
        self.hiplot=None
        self.showWoundCenters=False
        self.redDot=None
        self.bwDot=None
        self.selectionVals=[]
        self.pointSize=DEFAULT_SEED_SIZE
        self.point_mode=False
        self.mapPlotRandomArray=GetMapPlotRandomArray()
        self.mapPlotRandomFloatArray = np.array(self.mapPlotRandomArray,dtype=np.float)/255.
        self.mapPlotCmap = matplotlib.colors.ListedColormap(self.mapPlotRandomFloatArray.T)
        #self.mapPlotCmap = MakeCmapFromArray(self.mapPlotRandomFloatArray)
        
        self.previousDrawPoint=None
        
        # For undo...
        self.oldSparseList = None
        self.oldSeedSelections = None
        self.old_point_mode=self.point_mode
        #self.UpdateSeeds()
    
    def SetUndoPoint(self):
        # Set old state to be current state before any action
        # This allows undo to work...
        self.oldSparseList = deepcopy(self.sparseList[self.index])
        self.oldSeedSelections = deepcopy(self.seedSelections[self.index])
        self.old_point_mode=self.point_mode
        self.previousDrawPoint=None
    def RemoveUndo(self):
        # Remove old state
        self.oldSparseList = None
        self.oldSeedSelections = None
        self.old_point_mode=self.point_mode
    def Undo(self):
        # Swap old and current state
        if None not in [self.oldSparseList,self.oldSeedSelections]:
            print 'Undoing'
            print 'SL'
            print self.sparseList[self.index].nnz
            print 'oSL'
            print self.oldSparseList.nnz
            self.sparseList[self.index], self.oldSparseList = self.oldSparseList, self.sparseList[self.index] # SWAP
            self.seedSelections[self.index], self.oldSeedSelections = self.oldSeedSelections, self.seedSelections[self.index] # SWAP
            self.point_mode, self.old_point_mode = self.old_point_mode,self.point_mode # SWAP
            self.previousDrawPoint=None
    def SaveTxt(self,filename):
        fid = open(filename,'w')
        for i in txtHeader:
            fid.write(i+'\r\n')
        
        for i,t in enumerate(self.notes):
            if t not in ['',u'',None]:
                fid.write(txtsepA+str(i)+txtsepB+'\r\n')
                t2=t.replace(os.linesep,'\n').replace('\r','').replace('\n','\r\n')
                fid.write(t2+'\r\n')
    def LoadTxt(self,filename):
        print 'Load Notes'
        fid = open(filename,'r')
        for i in range(3):
            line=fid.readline().replace(os.linesep,'').replace('\n','').replace('\r','')
            if line==txtHeader[i]:
                pass
            else:
                print "Invalid Header"
                print line,"is not",txtHeader[i]
                return
        fnLast=-1
        newFrame=False
        self.notes=[]
        lA=len(txtsepA)
        lB=len(txtsepB)
        for line in fid:
            line=line.replace(os.linesep,'').replace('\n','').replace('\r','')
            # Check for Header line
            if len(line)>lA+lB:
                if line[:lA]==txtsepA and line[-lB:]==txtsepB and line[lA:-lB].isdigit():
                    fn=int(line[lA:-lB])
                    if fn>fnLast:
                        self.notes+=['']*(fn-fnLast) # Add entries for each frame to notes
                        fnLast=fn
                        newFrame=True
                    else:
                        print 'Invalid Frame Separators'
                        print 'Frame Separators must be in chronological order!'
                        return
            if fnLast==-1:
                print 'Should be a Frame Separator right after the Header!'
                return
            if newFrame:
                newFrame=False
            else:
                self.notes[-1]+=line+os.linesep
        self.notes+=[None]*(self.length-len(self.notes))
        print self.notes
    def Save(self,d,saveOutlines=True,saveMapImages=True): # Ignore Undo
        print 'Saving To ',d
        # TODO: FOR SOME WEIRD REASON, SAVE/LOAD only recreates the seeds
        #       for the first and last saved frame...
        # TODO: I need a simple indicator to show whether or not a frame
        #       has been initiated...
        if not os.path.exists(d):
            os.mkdir(d)
        
        segmentsD = os.path.join(d,'Segments')
        outlinesD = os.path.join(d,'Outlines')
        mapsD = os.path.join(d,'Maps')
        
        if not os.path.exists(segmentsD):
            os.mkdir(segmentsD)
        if not os.path.exists(outlinesD):
            os.mkdir(outlinesD)
        if not os.path.exists(mapsD):
            os.mkdir(mapsD)
        
        segmentsBase = os.path.join(segmentsD,'Segment')
        # Py makes it easy to load ;)
        seedPointsFile = os.path.join(d,'Seeds.py')
        notesFile = os.path.join(d,'Notes.txt')
        outlinesBase = os.path.join(outlinesD,'Outline')
        mapBase = os.path.join(mapsD,'Map')
        
        GTL.SaveFileSequence(self.watershed,basename=segmentsBase,
                             format='tif',sparseSave=self.framesVisited)
        
        if saveOutlines:
            getOutlines = lambda arr: CreateOutlines(arr,walgorithm='PyMorph') # This will break the deprecated OpenCV watershed
            GTL.SaveFileSequence(self.watershed,basename=outlinesBase,
                                 format='tif',sparseSave=self.framesVisited,functionToRunOnFrames=getOutlines)
        
        # Still just want an easy format to save and load...
        cooList = [  ( None if i==None else i.tocoo() )  for i in self.sparseList  ]
        seedList = [  ( None if i==None else np.array([i.row,i.col]).T.tolist())  for i in cooList  ]
        seedVals = [  ( None if i==None else i.data.astype(np.int).tolist() ) for i in cooList  ]
        
        seedListStr = repr(seedList)#.replace('[[[','['+os.linesep+'[[') \
                                         #.replace('[[','[ [') \
                                         #.replace(']]]','] ]'+os.linesep+']') \
                                         #.replace(']], ','] ],'+os.linesep) \
                                         #.replace('], ','],'+os.linesep+'  ')
        seedValsStr = repr(seedVals)#.replace('[[','['+os.linesep+'[') \
                                         #.replace(']',']'+os.linesep+']') \
                                         #.replace(', ',','+os.linesep+'  ')
        walgorithmStr = repr(self.walgorithm)
        woundCentersStr = repr(self.woundCenters)
        
        fid=open(seedPointsFile,'w')
        fid.write('seedList = '+seedListStr.replace('\r','').replace('\n',''))
        fid.write('\r\n')
        fid.write('seedVals = '+seedValsStr.replace('\r','').replace('\n',''))
        fid.write('\r\n')
        fid.write('walgorithm = '+walgorithmStr.replace('\r','').replace('\n',''))
        fid.write('\r\n')
        fid.write('woundCenters = '+woundCentersStr.replace('\r','').replace('\n',''))
        fid.write('\r\n')
        fid.close()
        print 'Saving seeds for '+str(self.length)+' frames:'
        for i,s in enumerate(self.sparseList):
            if s==None:
                print i,'--'
            elif s.nnz==0:
                print 'Empty!'
            else:
                print i,'initialized'
        
        # Save Notes
        self.SaveTxt(notesFile)
        
        # Save Stats
        self.SaveAllStats(d)
        
        # Save Map Images
        if saveMapImages:
            oldIndex=self.index
            plt.ioff()
            for i in range(self.length):
                self.index=i
                if self.framesVisited[i]:
                    self.MapPlot(saveFile=mapBase+str(i)+'.png')
            self.index = oldIndex
        
        self.framesVisited=[False]*self.length
        self.framesVisited[self.index]=True
        
    def Open(self,d):
        segmentsD = os.path.join(d,'Segments')
        seedPointsFile = os.path.join(d,'Seeds.py')
        notesFile = os.path.join(d,'Notes.txt')
        
        self.watershed[:]=0
        self.woutline[:]=0
        
        if not os.path.exists(segmentsD):
            os.mkdir(segmentsD)
        watershedTemp = GTL.LoadFileSequence(segmentsD,'Segment*')
        if watershedTemp!=None:
            for i in range(min(len(self.watershed),len(watershedTemp))):
                self.watershed[i] = watershedTemp[i]
                #self.CreateOutlines(i)
        
        print self.shape
        print self.watershed[0].max()
        
        #fid = open(seedPointsFile)
        #try:
        #    exec(fid.read().replace('\r','')) # Loads seedList and seedVals from the file...
        #except SyntaxError:
        #    print 'Invalid Syntax in Seeds.py'
        #fid.close()
        #Seeds.seedList = seedList
        #Seeds.seedVals = seedVals
        
        
        # use python's load_module function to get the data instead...
        fid=open(seedPointsFile,'U')
        Seeds = imp.load_module('Seeds',fid,'Seeds.py',('.py','U',1))
        fid.close()
        
        # A clunkier way to do it...
        #oldPath=sys.path
        #sys.path = [os.path.split(seedPointsFile)[0]]
        #imp.reload(Seeds)
        #sys.path=oldPath
        
        if len(Seeds.seedList)==self.length:
            for i in range(self.length):
                if None in [Seeds.seedList[i],Seeds.seedVals[i]]:
                    self.sparseList[i] = None
                else:
                    row,col = np.array(Seeds.seedList[i]).T
                    vals = Seeds.seedVals[i]
                    self.sparseList[i] = scipy.sparse.coo_matrix((vals,(row,col)), shape=self.shape[1:], dtype=np.uint16).tolil() # I guess I could change the dtype later if I need to...
            try: # Since walgorithm is not part of early versions, allow it to be optional
                Seeds.walgorithm
            except:
                Seeds.walgorithm=['PyMorph']*self.length
            
            try: # Since woundCenters was not part of version 0.2 and earlier, allow it to be optional
                Seeds.woundCenters
            except:
                Seeds.woundCenters=[None]*self.length
            
            self.walgorithm = Seeds.walgorithm
            self.woundCenters = Seeds.woundCenters
            
            print 'Loading seeds for '+str(self.length)+' frames:'
            for i,s in enumerate(self.sparseList):
                if s==None:
                    print i,'--'
                elif s.nnz==0:
                    print 'Empty!'
                else:
                    print i,'initialized'
            
            self.seedSelections = [None]*self.length
            for i in range(self.length):
                if self.sparseList[i]!=None:
                    self.seedSelections[i] = scipy.sparse.lil_matrix(self.shape[1:],dtype=np.bool)
            
            # This is technically unnecessary because it only gets called
            # on new WatershedData's anyway and those are always initialized to frame 0
            lastActiveFrame = self.GetLastActiveFrame()
            if self.index>lastActiveFrame:
                self.index=lastActiveFrame
            
            self.UpdateSeeds()
            
            # Load Notes
            if os.path.exists(notesFile):
                self.LoadTxt(notesFile)
            
            # Load Stats
            self.LoadStats(d)
            
            self.ColorPlot()
            self.MapPlot()
        else:
            print 'These seeds do not match the images open!'
            print 'Images length:',self.length,'  Seeds length:',len(Seeds.seedList)
        
        self.framesVisited=[False]*self.length
        self.framesVisited[self.index]=True
        
        self.RemoveUndo()
        
    def UpdateSeeds(self,force=False):
        '''Update seedArray from either local minima or seedList'''
        if force==True:
            self.SetUndoPoint()
            maxMap = MaxMinFinder(gaussian_filter(self.origData[self.index],self.gaussSigma)) # Max pts of gauss filtered image (single pixels with values of 1)
            row,col = np.where(maxMap)
            self.seedArray[:]=0
            for i in range(len(row)):
                self.seedArray[row[i],col[i]] = i+2 # 0 and 1 are reserved for unfilled region and background value respectively
            self.seedArray = scipy.ndimage.morphology.grey_dilation(self.seedArray,footprint=ImageCircle(DEFAULT_SEED_SIZE))
            self.sparseList[self.index] = scipy.sparse.lil_matrix(self.seedArray,dtype=np.uint16)
            self.seedSelections[self.index] = scipy.sparse.lil_matrix(self.shape[1:],dtype=np.bool)
        else:
            self.seedArray = self.sparseList[self.index].toarray()
        
        if self.background:
            val = 2 #(2 if self.walgorithm=='cv' else 1) # I could do this, but what's the need?
            self.seedArray[:val,:]  = 1
            self.seedArray[-val:,:] = 1
            self.seedArray[:,:val]  = 1
            self.seedArray[:,-val:] = 1
        self.Watershed()
    def Watershed(self,algorithm='PyMorph'): # Ignore Undo
        if not HAS_CV and algorithm in 'OpenCV':
            print 'OpenCV is not installed!'
            print 'Use PyMorph algorithm instead!'
            return
        if algorithm in ['OpenCV','PyMorph','Dummy']: # Later add CV2PM
            self.walgorithm[self.index] = algorithm
        else:
            print algorithm,'is not a valid watershed algorithm!'
            return
        # Run the real watershed algorithm
        if self.walgorithm[self.index]=='OpenCV':
            self.watershed[self.index] = cvWater(GTL.DivideConvertType(self.filterData[self.index],8),
                                             self.seedArray)
        elif self.walgorithm[self.index]=='PyMorph':
            #self.watershed[self.index] = pmWatershed.cwatershed(
            # Change to use mahotas ... should be faster...
            self.watershed[self.index] = mahotas.cwatershed(
                                             self.filterData[self.index],
                                             self.seedArray)
        elif self.walgorithm[self.index]=='Dummy':
            # Just set watershed directly to seeds
            self.watershed[self.index]=self.seedArray
            # Change 0 (unfilled) to 1 (background)
            w=np.where(self.watershed[self.index]==0)
            self.watershed[self.index][w]=1
        else:
            print self.walgorithm[self.index],'is not a valid watershed algorithm!'
        self.woutline = CreateOutlines(self.watershed[self.index],
                                       walgorithm=self.walgorithm[self.index])
    def UpdateValuesList(self,index=None):
        if index==None:
            index=self.index
        self.valList=[]
        for v in range(2,self.watershed[index].max()+1):
            if np.any(self.watershed[index]==v):
                self.valList.append(v)
    
    def GetCentroids(self,index,doUpdate=True): # Returns a floating point value...
        if doUpdate:
            self.UpdateValuesList(index)
        cmList = []
        for v in self.valList:
            cmList.append( center_of_mass((self.watershed[index]==v).astype(np.float)) )
        return cmList
    
    def GetTripleJunctionImage(self,index):
        'GetTripleJunctionImage generates a boolean array True at all points\n'
        'where 4 pixels meet with at least 3 distinct values.'
        sh = [i-1 for i in self.watershed[index].shape]
        tripleJ = np.zeros(sh,dtype=np.bool)
        for i in range(sh[0]):
            for j in range(sh[1]):
                if len(set([  self.watershed[index][i][j],
                              self.watershed[index][i+1][j],
                              self.watershed[index][i][j+1],
                              self.watershed[index][i+1][j+1]  ])) > 2:
                    tripleJ[i][j] = 1
        return tripleJ
    def GetTripleJunctionPointsAndVals(self,index):
        tj = np.array(np.where(self.GetTripleJunctionImage(index))).T
        idsByTJ=[]
        idsByTJArr = np.zeros([len(tj),4],np.int)
        for t in range(len(tj)):
            idsByTJ.append( list(set([self.watershed[index][tj[t][0]+x][tj[t][1]+y]
                                       for x,y in [[0, 0],[1, 0],[0, 1],[1, 1]]])))
            idsByTJArr[t][:len(idsByTJ[-1])] = idsByTJ[-1]
        self.UpdateValuesList(index=index)
        idList = self.valList
        tjsByID=[]
        for v in idList:
            tjsByID.append(np.where(idsByTJArr==v)[0].tolist())
        
        #      List of XY Coords of Triple Junctions
        #          List of all existing cellID's in this frame
        #                  For each triple junction, the 3 or 4 cellID's touch it
        #                           For each cell (by ID) which triple junctions touch it
        return tj, idList, idsByTJ, tjsByID # These last 2 index into each other...
    
    # If I ever decided to speed this up, the best way would be to pre-calculate all the pair distances
    # and then find the shortest combined path from summing these up
    # And I'm also sure someone has an algorithm that does it even better...
    def GetPolygonsAndPerimeters(self):
        lastActiveFrame = self.GetLastActiveFrame()
        tjPolys = []
        polyPerimeters = []
        edgeCells = []
        for index in range(lastActiveFrame+1):
            tj, idList, idsByTJ, tjsByID = self.GetTripleJunctionPointsAndVals(index) # Runs an Update
            tjPolys.append([None]*len(self.valList))
            polyPerimeters.append([None]*len(self.valList))
            for i,v in enumerate(self.valList): # same as idList
                tjs = [ tj[j] for j in tjsByID[i] ]
                tjPolys[-1][i] = GetMinPoly(tjs)
                polyPerimeters[-1][i] = CalcPolyPerimeter(tjPolys[-1][i])
         # so these, as with area, perimeter, etc, come out indexed by value *index*, not by value
         # use tjPolys[i], NOT tjPolys[v]
        return tjPolys,polyPerimeters
    
    def MakeSeedsFromPrevious(self):
        if self.index>0:
            self.SetUndoPoint()
            self.seedArray[:] = 0
            cmList = self.GetCentroids(self.index-1)
            for cm in cmList:
                cm = [int(round(cm[0])),int(round(cm[1]))] # Convert to int
                self.seedArray[cm[0],cm[1]] = self.watershed[self.index-1][cm[0],cm[1]]
            self.seedArray = scipy.ndimage.morphology.grey_dilation(self.seedArray,footprint=ImageCircle(DEFAULT_SEED_SIZE))
            self.sparseList[self.index] = scipy.sparse.lil_matrix(self.seedArray)
            self.seedSelections[self.index] = scipy.sparse.lil_matrix(self.shape[1:],dtype=np.bool)
            self.UpdateSeeds()
            self.Watershed()
    def CopySeedsFromPrevious(self):
        if self.index>0:
            self.SetUndoPoint()
            self.sparseList[self.index]=deepcopy(self.sparseList[self.index-1])
            self.seedSelections[self.index]=deepcopy(self.seedSelections[self.index-1])
            self.UpdateSeeds()
            self.Watershed()
    def GetLastActiveFrame(self):
        lastActiveFrame=self.length-1
        for i in range(self.length):
            if self.sparseList[i]==None:
                lastActiveFrame = i-1
                break
        return lastActiveFrame
    
    def NextFrame(self,doPlots=True):
        if self.index+1<self.length:
            self.lastFrameVisited = self.index
            self.index+=1
            
            print 'Move to Frame',self.index
            if self.sparseList[self.index]==None:
                self.MakeSeedsFromPrevious()
            else:
                self.UpdateSeeds()
            self.RemoveUndo()
            if doPlots: # If we skip over it, no need to save...
                self.framesVisited[self.index]=True
                self.MapPlot()
                self.ColorPlot()
        else:
            print 'No more frames!!!!'
    
    def PreviousFrame(self,doPlots=True):
        if self.index>0:
            self.lastFrameVisited = self.index
            self.index-=1
            print 'Move to Frame',self.index
            self.UpdateSeeds()
        else:
            print 'At the beginning.'
        self.RemoveUndo()
        if doPlots: # If we skip over it, no need to save...
            self.framesVisited[self.index]=True
            self.MapPlot()
            self.ColorPlot()
    def MoveToFrame(self,newIndex):
        if newIndex<0:
            newIndex=0
        if newIndex>=self.length:
           newIndex=self.length-1

        lastActiveFrame = self.GetLastActiveFrame()
        if newIndex<lastActiveFrame+1:
            self.lastFrameVisited = self.index
            self.index=newIndex # hop ahead, then find new seeds...
            print 'Move to Frame',self.index
            self.UpdateSeeds()
            self.RemoveUndo()
        elif newIndex==lastActiveFrame:
            self.lastFrameVisited = self.index
            self.index=newIndex # hop ahead, then find new seeds...
            print 'Move to Frame',self.index
            self.MakeSeedsFromPrevious()
            self.RemoveUndo()
        else:# newIndex>lastActiveFrame+1
            print "You can't move that far ahead yet!"
        
        self.framesVisited[self.index]=True
        self.MapPlot()
        self.ColorPlot()
    def GoToLastVisitedFrame(self):
        self.MoveToFrame(self.lastFrameVisited)
    def Invert(self):
        self.filterData[self.index] = self.filterData[self.index].max()-self.filterData[self.index]
        self.RemoveUndo()
    def Gauss(self,rad):
        '''Apply Gaussian Blur to data'''
        self.gaussSigma=rad
        self.RemoveUndo()
    def ResetData(self):
        '''Return data to original data'''
        self.gaussSigma=0
        self.filterData[self.index] = np.array(self.origData[self.index])
        self.RemoveUndo()
    def BGSubtract(self,rad):
        '''Apply Gaussian Blur to data'''
        bg = gaussian_filter(self.filterData[self.index],rad)
        self.filterData[self.index] = self.filterData[self.index] + bg.max() - bg
        self.filterData[self.index] -= self.filterData[self.index].min()
        self.RemoveUndo()
    def Median(self,filterSize):
        self.filterData[self.index] = median_filter(self.filterData[self.index],filterSize)
        self.RemoveUndo()
    def Sharpen(self,r1,r2):
        '''Apply Gaussian Blur to data'''
        g1 = gaussian_filter(self.filterData[self.index],r1)
        g2 = gaussian_filter(self.filterData[self.index],r2)
        self.filterData[self.index] =  g1 + g2.max() - g2
        self.filterData[self.index] -= self.filterData[self.index].min()
        self.RemoveUndo()
    def OutOfBounds(self,point):
        s=self.origData[self.index].shape
        return ( point[0]<0 or point[1]<0 or point[0]>=s[0] or point[1]>=s[1] )
        
    def UpdatePointsWithVal(self,wh,val): # Raw update given points and values
        self.seedArray[wh] = val
        self.sparseList[self.index] = scipy.sparse.lil_matrix(self.seedArray,dtype=np.uint16)
        self.seedSelections[self.index] = scipy.sparse.lil_matrix(self.shape[1:],dtype=np.bool)
        self.UpdateSeeds()
    def NewSeed(self,pointIn,val=None):
        '''Add a new seed point'''
        if self.OutOfBounds(pointIn):
            return False
        self.SetUndoPoint()
        
        newPoints = GetPointsAtRadius(self.seedArray.shape,pointIn[0],pointIn[1],self.pointSize)
        
        if val==None:
            m = 2 # 0 is unassigned and 1 is for background...
            for l in self.sparseList:
                if l!=None:
                    if l.nnz!=0:
                        m = max( m, l.tocoo().data.max() )
            val = m+1
        
        points = np.array(newPoints).T
        wh = (points[0],points[1])
        self.UpdatePointsWithVal(wh,val)
        
        return True
    def ExtraSeed(self,pointIn,val):
        return self.NewSeed(pointIn,val=val)
    
    def ExtraSeedLine(self,point0,point1,val):
        newPoints = BresenhamFunction(point0, point1)
        if self.pointSize==1: # Make 1-pixel lines
            pass
        else:                 # Make 2-pixel lines
            dilateByCircle(newPoints,self.watershed[0].shape,1.9)
        
        points = np.array(newPoints).T
        wh = (points[0],points[1])
        self.UpdatePointsWithVal(wh,val)
    
    def DeleteSeedByRegion(self,point):
        '''Remove the seed point from the list'''
        self.SetUndoPoint()
        if self.watershed[self.index]!=None:
            val=self.watershed[self.index,point[0],point[1]]
            wh = np.where(self.seedArray==val)
            self.UpdatePointsWithVal( wh, 0)
            return True
        else:
            return False
    def DeleteSeed(self,point):
        '''Remove the seed point from the list'''
        self.SetUndoPoint()
        if self.seedArray[point[0],point[1]]!=0:
            self.UpdatePointsWithVal( ( np.array([point[0]]) , np.array([point[1]]) ), 0)
            return True
        else:
            return False
    def DeleteSelectedSeeds(self,invertSelections=False):
        self.SetUndoPoint()
        wh = np.where( self.seedSelections[self.index].toarray() ^ invertSelections )
        self.UpdatePointsWithVal(wh,0)
        
    def MergeSelectedSeeds(self):
        self.SetUndoPoint()
        
        # get the val to change to:
        wh = np.where(self.seedArray * self.seedSelections[self.index].toarray())
        newVal = min(self.seedArray[wh])
        self.UpdatePointsWithVal(wh,newVal)
    
    def ChangeRegionValue(self,vOld,vNew):
        # I think it is probaly the right thing to add undo here, but it sure
        # Does make things quirky b/c only changes in seeds are recorded in Undo...
        #self.SetUndoPoint()
        self.RemoveUndo()
        wh = np.where(self.seedArray==vOld)
        self.UpdatePointsWithVal(wh,vNew)
    def SwitchRegionValues(self,v1,v2):
        # I think it is probaly the right thing to add undo here, but it sure
        # Does make things quirky b/c only changes in seeds are recorded in Undo...
        #self.SetUndoPoint()
        
        self.RemoveUndo()
        wh1 = np.where(self.seedArray==v1)
        wh2 = np.where(self.seedArray==v2)
        self.UpdatePointsWithVal(wh1,v2)
        self.UpdatePointsWithVal(wh2,v1)
    #def DrawOrigImage(self): # deprecated
    def MapPlot(self,saveFile=None,useText=False):
        plt.figure(2)#;cla()
        
        wi = np.array(self.watershed[self.index],dtype=np.int)
        
        if saveFile!=None and useText: # Normally skip all the text labels... they are painfully slow!
            for i in self.mapCentroids:
                if i!=None:
                    i.set_visible(False)
            
            for i in self.GetActiveSeedValues():
                if i>0:
                    cm = center_of_mass((wi==i).astype(np.float))
                    if math.isnan(cm[0]) or math.isnan(cm[1]):
                        print 'CM is NAN'
                        print 'val',i
                        print cm[0]
                        print cm[1]
                    else:
                        if i>=len(self.mapCentroids):
                            self.mapCentroids+=[None]*(1+i-len(self.mapCentroids))
                        x,y = int(round(cm[0])),int(round(cm[1]))
                        if self.mapCentroids[i]!=None:
                            self.mapCentroids[i].set_position([y-5,x+5])
                            self.mapCentroids[i].set_visible(True)
                        else:
                            self.mapCentroids[i] = plt.text(y-5,x+5,str(i),fontsize=10)
                            #self.mapCentroids[i].set_visible(True)
                        # Never actually remove the objects, just make them invisible...
                        # Sneaky...
                        
            plt.draw()
            plt.savefig(saveFile)
            plt.cla()
        elif saveFile!=None: # I may have broken this... sorry...
            f1=np.vectorize(lambda x: self.mapPlotRandomArray[0][x])
            f2=np.vectorize(lambda x: self.mapPlotRandomArray[1][x])
            f3=np.vectorize(lambda x: self.mapPlotRandomArray[2][x])
            
            if self.rgbM==None:
                self.rgbM = np.array([f1(wi),f2(wi),f3(wi)]).transpose(1,2,0)
            else:
                # No more dependence on cython function; save will be a little slower, but whatever...
                #convToRandColors(self.mapPlotRandomArray,self.rgbM,
                #                 wi,wi.shape[0],wi.shape[1])
                # This step is still really slow...
                self.rgbM[:,:,0] = f1(self.watershed[self.index])
                self.rgbM[:,:,1] = f2(self.watershed[self.index])
                self.rgbM[:,:,2] = f3(self.watershed[self.index])
            
            im = Image.fromarray(self.rgbM)
            im.save(saveFile)
            #imsave(saveFile, map)
        else:
            if self.mapPlot!=None:
                ###self.mapPlot.set_data(self.rgbM)
                self.mapPlot.set_data(self.watershed[self.index])
                
            else:
                ###self.mapPlot = plt.imshow(self.rgbM,interpolation='nearest',animated=True)
                self.mapPlot = plt.imshow(self.watershed[self.index],animated=True,
                                          interpolation='nearest',cmap=self.mapPlotCmap,
                                          norm=matplotlib.colors.NoNorm() )
            self.DrawBWDot()
            #plt.draw() # Now DrawBWDot calls this instead...
    def MapPlotWTracks(self):
        # Disable the old mapPlot and redraw every time (will be slow...)
        # Later, could possibly store the tracks part for animations, too...
        self.mapPlot=None
        self.MapPlot()
        if hasattr(self,'centroidX'):
            x=np.array(self.centroidX)
            y=np.array(self.centroidY)
        else:
            self.UpdateValuesList(index)
            x=[]
            y=[]
            for index in range(self.length):
                if self.sparseList[index]==None:
                    break
                elif self.sparseList[index].nnz==0:
                    break
                centroid = self.GetCentroids(index,doUpdate=False)
                x.append([i[0] for i in centroid])
                y.append([i[1] for i in centroid])
            x=np.array(x)
            y=n.array(y)
        plt.plot(y[0],x[0],'kx')
        plt.plot(y,x,'k')
        # initiate the mouse-over printing...
        plt.gca().format_coord = GetReportPixel(self)
    def ToggleOverlaysVisible(self):
        self.overlayVisible = not self.overlayVisible
        self.ColorPlot()
    def ColorPlot(self):
        plt.figure(1)
        im = self.filterData[self.index] # NOT 8- bit!!!
        seed = (self.seedArray!=0).astype(np.uint8)
        outl = self.woutline.astype(np.uint8)
        seedORoutl = np.array(np.array(seed+outl)>0,np.uint8)
        if self.rgba==None:
            self.rgba=np.array([0*seed,
                          255*seed,
                          255*outl,
                          self.overlayVisible*255*seedORoutl]).transpose(1,2,0)
        else:
            self.rgba[:,:,0]=0*seed
            self.rgba[:,:,1]=255*seed
            self.rgba[:,:,2]=255*outl
            self.rgba[:,:,3]=self.overlayVisible*255*seedORoutl
        # 0.18s
        if self.colorplot!=None and self.bgPlot!=None:
            self.bgPlot.set_data(im)
            self.colorplot.set_data(self.rgba)
        else:
            self.bgPlot = plt.imshow(im,cmap=plt.cm.gray,interpolation='nearest',animated=True)
            self.colorplot = plt.imshow(self.rgba,interpolation='nearest',animated=True)
            plt.xlim(0,self.watershed[0].shape[1])
            plt.ylim(self.watershed[0].shape[0],0)
        
        dprint(['wc',self.woundCenters[self.index],self.showWoundCenters])
        # Need to force this to plot so that axes will flip, just set to invisible...
        if self.redDot==None:
            self.redDot=plt.plot([0],[0],'ro')
            plt.xlim(0,self.watershed[0].shape[1])
            plt.ylim(self.watershed[0].shape[0],0)
        
        if self.showWoundCenters and self.woundCenters[self.index]!=None:
            [x,y] = self.woundCenters[self.index]
            self.redDot[0].set_data([x],[y])
            self.redDot[0].set_visible(True) # Turn On Visible
        else:
            self.redDot[0].set_visible(False) # Turn Off Visible
        #plt.xlim(0,self.watershed[0].shape[1])
        #plt.ylim(self.watershed[0].shape[0],0)
        # AAARRRGGG!!!
        # NEED TO FORCE COORDINATES TO STAY IMAGE-STYLE!
        
        # This is by FAR the slowest step, about 0.25s
        self.HighlightRegion()
        # plt.draw() # Highlight Region also calls draw, so this doesn't need to!
    def UpdateSelection(self,point=None,append=False):
        if append==False:
            self.selectionVals=[]
        if point!=None:
            self.selectionVals+=[self.watershed[self.index,point[0],point[1]]]
        
        self.point_mode=False
        
    def DrawBWDot(self):
        fn=plt.gcf().number
        plt.figure(2)
        
        if self.point_mode==False:
            cm=np.zeros([len(self.selectionVals),2],dtype=np.float)
            for i,v in enumerate(self.selectionVals):
                if v>1:
                    # centroid of selected region v
                    cm[i] = center_of_mass((self.watershed[self.index]==v).astype(np.float))
            if self.bwDot==None and self.point_mode==False:
                self.bwDot=plt.plot(cm[:,1],cm[:,0],'ko')
                self.bwDot[0].set_markeredgecolor('w')
                self.bwDot[0].set_markeredgewidth(1.2)
            else:
                self.bwDot[0].set_data(cm[:,1],cm[:,0])
        
        plt.draw()
        # Set the figure back to what it was
        plt.figure(fn)
    
    def HighlightRegion(self):
        a=np.zeros(self.shape[1:],dtype=np.uint8)
        if self.point_mode:
            a[:] = self.seedSelections[self.index].toarray().astype(np.uint8)*255
            
            if self.rgbaH==None:
                self.rgbaH=np.array([a,a*0,a,a]).transpose(1,2,0)
            else:
                self.rgbaH[:,:,0]=a
                self.rgbaH[:,:,1]=a*0
                self.rgbaH[:,:,2]=a
                self.rgbaH[:,:,3]=a
        else:
            for i in self.selectionVals:
                if i>0:
                    a += 255*(self.watershed[self.index]==i)
            if self.rgbaH==None:
                self.rgbaH=np.array([a,a,a*0,a//2]).transpose(1,2,0)
            else:
                self.rgbaH[:,:,0]=a
                self.rgbaH[:,:,1]=a
                self.rgbaH[:,:,2]=a*0
                self.rgbaH[:,:,3]=a//2
            self.DrawBWDot()
        
        if self.hiplot!=None:
            self.hiplot.set_data(self.rgbaH)
        else:
            self.hiplot = plt.imshow(self.rgbaH,interpolation='nearest',animated=True)
        
        plt.draw()
        
    def LassoCallback(self,verts):
        verts = [ v[::-1] for v in verts ]
        
        self.SetUndoPoint()
        
        # Clear selections...
        self.seedSelections[self.index] = scipy.sparse.lil_matrix(self.shape[1:],dtype=np.bool)
        
        cooMat = self.sparseList[self.index].tocoo()
        seedList = np.array([cooMat.row,cooMat.col]).T
        wh = np.where(points_inside_poly(seedList, verts))
        row = np.array([cooMat.row[i] for i in wh])
        col = np.array([cooMat.col[i] for i in wh])
        
        # This should be faster for big selections
        sel = self.seedSelections[self.index].toarray()
        sel[(row,col)] = 1
        self.seedSelections[self.index] = scipy.sparse.lil_matrix(sel)
        
        self.point_mode=True
        
        self.ColorPlot()
        
        plt.figure(1).canvas.draw_idle()
        plt.figure(1).canvas.widgetlock.release(self.lasso)
        del self.lasso
    def CompressSeedValues(self):
        oldIndex = self.index
        
        numFrames=self.length
        for i in range(self.length):
            if self.sparseList[i]==None:
                numFrames = i
                break
        
        usedVals=[]
        for v in range(2,self.watershed.max()+1):
            if np.any(self.watershed==v):
                usedVals.append(v)
        usedVals.sort() # sort the list
        
        print usedVals
        
        newVal = 2
        for val in usedVals:
            if val!=newVal:
                for i in range(numFrames):
                    print 'frame',i
                    self.index=i
                    
                    if self.oldSparseList!=None:
                        oldSeedArray = self.oldSparseList.toarray()
                        oldSeedArray[ np.where(oldSeedArray==val) ] = newVal
                        self.oldSparseList = scipy.sparse.lil_matrix(oldSeedArray,dtype=np.uint16)
                    
                    self.seedArray[ np.where(self.seedArray==val) ] = newVal
                    self.sparseList[i] = scipy.sparse.lil_matrix(self.seedArray,dtype=np.uint16)
                    
                    if self.selectionVals not in [[],None]:
                        while val in l:
                            l[l.index(val)]=newVal
                    
                    self.UpdateSeeds() # takes care of watershed...
                    self.UpdateValuesList()
                    self.UpdateSelection()
                    
                    self.framesVisited[i]=True
            newVal+=1
        self.index = oldIndex
    def AutoCenterWound(self,woundVals):
        oldIndex=self.index
        for self.index in range(self.length): # For each frame
            if self.sparseList[self.index]!=None:
                newWVal = 1e6 # I don't think anyone will ever create a million cells...
                wi = np.array(self.watershed[self.index],dtype=np.int)
                for v in woundVals:
                    wi[np.where(wi==v)] = newWVal
                test = (wi==newWVal)
                if np.any(test):
                    self.woundCenters[self.index] = center_of_mass(test.astype(np.float))[::-1]
        self.index=oldIndex
    def GetBoundingRectangles(self,doUpdate=True):
        if doUpdate:
            self.UpdateValuesList()
        boundsList = []
        for v in self.valList:
            x,y=np.where(self.watershed[self.index]==v)
            boundsList.append( [[min(x),max(x)],[min(y),max(y)]] )
        return boundsList
    def CalculateArea(self,doUpdate=True):
        if doUpdate:
            self.UpdateValuesList()
            #boundsList = self.GetBoundingRectangles() # not really used...
        areaList = []
        for v in self.valList:
            areaList.append( np.sum(self.watershed[self.index]==v) )
        return areaList
    def CalculatePerimeter(self,boundsList=None):
        if boundsList==None:
            self.UpdateValuesList()
            boundsList = self.GetBoundingRectangles()
        perimeterList = []
        for i,v in enumerate(self.valList):
            perimeterList.append(ImageContour.GetIJPerimeter(self.watershed[self.index],v,boundsList[i]))
        return perimeterList
    def CalculateBestFitEllipse(self,boundsList=None):
        if boundsList==None:
            self.UpdateValuesList()
            boundsList = self.GetBoundingRectangles()
        ellipseList = []
        for i,v in enumerate(self.valList):
            [[xm,xM],[ym,yM]] = boundsList[i]
            ellipseList.append(EllipseFitter.EllipseFitter(self.watershed[self.index][xm:xM+1,ym:yM+1]==v))
        return ellipseList
    def CalculateCentroids(self,doUpdate=True):
        cmList = self.GetCentroids(self.index,doUpdate=doUpdate)
        return cmList
    def CalculateWoundDistance(self,cmList=None):
        # possibly open a user window to select wound location???
        # should run from the SegmenterFrame and not WatershedData,though...
        if cmList==None:
            cmList = self.GetCentroids(self.index)
        
        # Need to flip these
        woundY,woundX = self.woundCenters[self.index]
        
        rList = []
        thetaList = []
        for [x,y] in cmList:
            r = np.sqrt((x-woundX)**2+(y-woundY)**2)
            theta = np.arctan2(1.*(x-woundX),1.*(y-woundY))
            rList.append(r)
            thetaList.append(theta)
        return rList, thetaList
    def CalculateNeighbors(self,boundsList=None):
        if boundsList==None:
            self.UpdateValuesList()
            boundsList = self.GetBoundingRectangles()
        neighborList=[]
        for i,v in enumerate(self.valList):
            neighborList.append([])
            [[xm,xM],[ym,yM]] = boundsList[i]
            if xm==0 or ym==0 or xM==self.watershed.shape[1]-1 or yM==self.watershed.shape[2]-1:
                # cell is on the boundary
                # this should never happen, because background is forced onto all of the edges
                neighborList[i].append(0)
            waterSmaller=self.watershed[self.index][xm-1:xM+2,ym-1:yM+2]
            binary = (waterSmaller==v)
            neighborPixels = (scipy.ndimage.binary_dilation(binary)-binary)*waterSmaller
            #Loop over neighborPixels just like in UpdateValuesList, but also look at background
            for neighborVal in range(1,neighborPixels.max()+1):
                if np.any(neighborPixels==neighborVal):
                    neighborList[i].append(neighborVal)
        return neighborList
    def GetMaxSeedVal(self):
        maxSeedVal=2
        for l in self.sparseList:
            if l!=None:
                if l.nnz!=0:
                    maxSeedVal = max( maxSeedVal, l.tocoo().data.max() )
        return maxSeedVal
    def CollectAllStats(self):
        maxSeedVal = self.GetMaxSeedVal()
        
        self.valuesByFrame=[]
        self.centroidX=[]
        self.centroidY=[]
        self.xmin=[]
        self.xmax=[]
        self.ymin=[]
        self.ymax=[]
        self.area=[]
        self.perimeter=[]
        self.major=[]
        self.minor=[]
        self.angle=[]
        self.woundDistance=[]
        self.woundAngle=[]
        self.neighbors=[]
        oldIndex=self.index
        adjLength = self.GetLastActiveFrame()+1
        for index in range(adjLength):
            self.index=index
            
            print 'Calculating for Frame',index
            print 'Get Values'
            self.UpdateValuesList()
            self.valuesByFrame.append(self.valList)
            print 'Get Centroids'
            centroid = self.CalculateCentroids(doUpdate=False)
            self.centroidX.append([i[0] for i in centroid])
            self.centroidY.append([i[1] for i in centroid])
            print 'Get Bounds'
            boundsList=self.GetBoundingRectangles(doUpdate=False)
            br=np.array(boundsList)
            self.xmin.append(br[:,0,0].tolist())
            self.xmax.append(br[:,0,1].tolist())
            self.ymin.append(br[:,1,0].tolist())
            self.ymax.append(br[:,1,1].tolist())
            print 'Get Area'
            self.area.append(self.CalculateArea(doUpdate=False))
            print 'Get Perimeter'
            self.perimeter.append(self.CalculatePerimeter(boundsList=boundsList))
            print 'Get BF Ellipse'
            bfe = np.array(self.CalculateBestFitEllipse(boundsList=boundsList))
            self.major.append(bfe[:,1].tolist())
            self.minor.append(bfe[:,2].tolist())
            self.angle.append(bfe[:,0].tolist())
            rList,thetaList=self.CalculateWoundDistance(cmList=centroid)
            self.woundDistance.append(rList)
            self.woundAngle.append(thetaList)
            print 'Get Neighbors'
            self.neighbors.append(self.CalculateNeighbors(boundsList=boundsList))
            
            print 'Fix Offsets'
            # Make it so that there is an entry for every value, regardless of whether that cell was
            # in that frame or not
            # This makes all the quantities line up...
            for v in range(2,maxSeedVal+1): # This can get ugly if there is too much deleting and what-not...
                                            # Use the CompressSeedVals Button to remedy the situation...
                for i in range(len(self.valuesByFrame)):
                    if not v in self.valuesByFrame[i]:
                        self.valuesByFrame[i].insert(v-2,v)
                        self.centroidX[i].insert(v-2,None)
                        self.centroidY[i].insert(v-2,None)
                        self.xmin[i].insert(v-2,None)
                        self.xmax[i].insert(v-2,None)
                        self.ymin[i].insert(v-2,None)
                        self.ymax[i].insert(v-2,None)
                        self.area[i].insert(v-2,None)
                        self.perimeter[i].insert(v-2,None)
                        self.major[i].insert(v-2,None)
                        self.minor[i].insert(v-2,None)
                        self.angle[i].insert(v-2,None)
                        self.woundDistance[i].insert(v-2,None)
                        self.woundAngle[i].insert(v-2,None)
                        self.neighbors[i].insert(v-2,None)
        
        print 'Done Collecting!'
        self.index=oldIndex
    def SaveCSV(self,name,directory,data):
        fid=open(os.path.join(directory,name+'.csv'),'w')
        fid.write(repr(data).replace('\r','').replace('\n','')
                            .replace('[[','') \
                            .replace(']]','') \
                            .replace('[','') \
                            .replace('], ','\r\n'))
        fid.close()
    def SaveAllStats(self,directory=None):
        if directory==None:
            directory = wx.DirSelector()
        
        if hasattr(self,'centroidX'): # If it has centroid, it should have everything else, too...
            if self.centroidX not in [None,[]]:  # This could change usefulness, so beware...
                data = [self.valuesByFrame, self.centroidX, self.centroidY,
                        self.xmin, self.xmax, self.ymin, self.ymax,
                        self.area, self.perimeter,
                        self.major, self.minor, self.angle,
                        self.woundDistance, self.woundAngle]
                        
                names = ['CellValuesByFrame','CentroidX','CentroidY',
                         'Xmin','Xmax','Ymin','Ymax',
                         'Area','Perimeter',
                         'MajorAxis','MinorAxis', 'Angle',
                         'DistanceToWound','AngleToWound']
                
                adjLength = self.GetLastActiveFrame()+1
                if adjLength > 256:
                    print 'Over 256 frames!'
                    print 'Ask user how to save xls file...'
                    teDialog=wx.TextEntryDialog(None,
        'There are too many frames ('+str(adjLength)+') for one excel sheet (max 256)!\n'
        'Enter the number of frames to skip between saving.\n'
        '(A value of 1 means save all the data, 2 is save every other frame, etc...)\n'
        'Otherwise, the data will be broken into 256 frame files.',
                                              caption='Skip Frames?',defaultValue='1')
                    dia=teDialog.ShowModal()
                    if dia==wx.ID_CANCEL:
                        print 'Cancel, save all frames'
                        val = 1
                    elif dia==wx.ID_OK:
                        val = teDialog.GetValue()
                        if val.isdigit():
                            val=int(val)
                            if val>1:
                                print 'Only save every',val,'frames.'
                            else:
                                val=1
                                print 'Only save every',val,'frames.'
                        else:
                            print 'Bad entry!  Just save all frames!'
                            val=1
                    if val>1:
                        data = deepcopy(data) # let's not overwrite our existing attributes, shall we?
                        for i in range(len(data)):
                            data[i] = data[i][::val]
                    numBreakups=1
                    if len(data[0])>256: # This will break the excel write, so make break into pieces...
                        numBreakups = (len(data[0])-1)/256 + 1
                        # this breaks the data into somewhat even chunks...
                        #   nf=(len(data[0])-1)//numBreakups+1 # number of frames of a typical breakup
                        # indexList = [nf*i for i in range(numBreakups)]+[len(data[0])]
                        
                        # but I think it might be just as well to do this instead:
                        indexList = [256*i for i in range(numBreakups)]+[len(data[0])]
                        
                        dataSets=[[]]*numBreakups # make a slot for each excel file...
                        for f in range(numBreakups): # each file
                            for d in data: # each parameter
                                dataSets[f].append(d[indexList[f]:indexList[f+1]])
                    else:
                        dataSets = [data]
                    for i,d in enumerate(dataSets): # write a file for each set, max 256 frames...
                        cutByStr = ('' if val==1 else 'CutBy'+str(val))
                        partStr = ('' if numBreakups==1 else 'Part'+str(i+1))
                        ExcelHelper.excelWrite(d,names,os.path.join(directory,'Calculations'+cutByStr+partStr+'.xls'),flip=True)
                else:
                    ExcelHelper.excelWrite(data,names,os.path.join(directory,'Calculations.xls'),flip=True)
                for i,name in enumerate(names):
                    self.SaveCSV(name,directory,data[i])
                
                # This won't excel-ify very well since it's a 3D data set with variable numbers of links...
                # So just save it as a .py file instead...
                neighborsStr = repr(self.neighbors)
                
                fid=open(os.path.join(directory,'Neighbors.py'),'w')
                fid.write('neighbors = '+neighborsStr)
                fid.close()
    def CheckForMalformedRegions(self):
        '''This goes through all frames and values and uses shapely to test if regions are disjoint in any way'''
        ### THIS FUNCTION IMPORTS SHAPELY DIRECTLY ... DO NOT LINK TO IT!
        ### NO NEED TO ADD ANOTHER DEPENDENCY FOR EVERYONE ELSE!
        import shapely.geometry
        for frame in range(wd.GetLastActiveFrame()+1):
            wd.UpdateValuesList(frame)
            for v in wd.valList:
                ic = ImageContour.GetContourInPlace(wd.watershed[frame],v)
                p=shapely.geometry.asPolygon(ic.tolist())
                if (wd.watershed[frame]==v).sum()!=p.area:
                     print 'frame:',frame,'value:',v,'-- Value has disjoint regions!'
                if not p.is_valid:
                    print 'frame:',frame,'value:',v,'--Value has regions connected by point only!'
    def SaveTripleJunctions(self,directory=None):
        '''This is horribly slow, so keep it out of RunCalc and RunCalc2!'''
        if directory==None:
            directory = wx.DirSelector()
        adjLength = self.GetLastActiveFrame()+1
        
        fid = open(os.path.join(directory,'TripleJunctions.py'),'w')
        fid.write('[')
        for i in range(adjLength):
            tj = np.array(np.where(self.GetTripleJunctionImage(i))).T.tolist()
            fid.write(repr(tj))
            if i!=adjLength-1:
                fid.write(', ')
        fid.write(']')
        fid.close()
    def GetActiveSeedValues(self):
        return sorted(list(set( self.sparseList[self.index].tocoo().data.tolist() )))
    def SaveTripleJunctionsWithCellIDs(self,directory=None):
        '''This is horribly slow, so keep it out of RunCalc and RunCalc2!'''
        if directory==None:
            directory = wx.DirSelector()
        adjLength = self.GetLastActiveFrame()+1
        bigData=[]
        ColNames = ['X','Y','Cell ID 1','Cell ID 2','Cell ID 3','(Cell ID 4)']
        for i in range(adjLength):
            tj = np.where(self.GetTripleJunctionImage(i))
            l=[]
            for t in range(len(tj[0])):
                cellIDs = list(set([self.watershed[i][tj[0][t]+x][tj[1][t]+y] for x,y in [[0, 0],[1, 0],[0, 1],[1, 1]]]))
                cellIDs = map(str,cellIDs)
                if len(cellIDs)==3:
                    cellIDs = cellIDs+['']
                l.append([str(tj[0][t]),str(tj[1][t])]+cellIDs)
            bigData.append([ColNames]+l)
        ExcelHelper.excelWrite(bigData,[str(i+1) for i in range(adjLength)],
                               os.path.join(directory,'TripleJunctionsWithCellIDs.xls'))
    #NOT TESTED
    def GetNeighborPairsByFrame(self,woundVals): # Merges wound vals... make sure you run LoadStats or CollectAllStats first!
        neighborPairsByFrame = []
        for frame in range(self.length):
            MergeWoundVals(self.neighbors[frame],woundVals)
            neighborPairs = set()
            for i in range(len(self.neighbors[frame])):
                if self.neighbors[frame][i]!=None:
                    if len(self.neighbors[frame][i])>0:
                        neighborPairs.update( [tuple(sorted([i+2,j])) for j in self.neighbors[frame][i] ] )
            neighborPairsByFrame.append(sorted(list(neighborPairs)))
        return neighborPairsByFrame
    #NOT TESTED
    def GetContourValuesLengthsAndSubContoursByFrame(self,allValsByFrame,woundValsSorted):
        ######### NOT DONE!!!!!!! ####################
        # THIS IS NOT THIS SIMPLE... STILL NEED TO SOMEHOW MERGE WOUND
        cVLSByFrame = []
        for frame in range(self.length):
            cVLS = [] # for each subcountour: [value, length, subcontour points]
            for v in allValsByFrame[frame]:
                wi = self.watershed[frame]
                for v in woundValsSorted[1:]:
                    watershed[np.where(watershed==v)]=woundValsSorted[0]
                boundingRect=ImageContour.GetBoundingRect(wi,v)
                # No longer needed: #contour,turns,vals = ImageContour.GetContour(wi,v,boundingRect=boundingRect,byNeighbor=True)
                perimeterVals,perimeterList,subContours = ImageContour.GetPerimeterByNeighborVal(wi,v,boundingRect=boundingRect,getSubContours=True)
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
    #NOT TESTED
    def SaveSubContours(self,directory=None):
        '''This is horribly slow, so keep it out of RunCalc and RunCalc2!'''
        if directory==None:
            directory = wx.DirSelector()
        
        MI = self.GetManualInputs(d)
        if MI==None:
            return    
        woundValsSorted = sorted(MI.woundVals)
        self.CollectAllStats() # Safer but slower way to do it...
        neighborPairsByFrame = self.GetNeighborPairsByFrame(woundValsSorted)
        allValsByFrame = GetAllValsByFrame(neighborPairsByFrame)
        
        cVLSByFrame = self.GetContourValuesLengthsAndSubContoursByFrame(allValsByFrame,woundValsSorted)
        fid=open(os.path.join(directory,'CellBoundarySubContours.py'),'w')
        fid.write('cVLSByFrame = '+repr(cVLSByFrame))
        fid.close()
    def LoadStats(self,directory=None):
        if directory==None:
            directory = wx.DirSelector()
        
        #calculationsXls = os.path.join(directory,'Calculations.xls')
        neighborsPy = os.path.join(directory,'Neighbors.py')
        
        names = ['CellValuesByFrame','CentroidX','CentroidY',
                 'Xmin','Xmax','Ymin','Ymax',
                 'Area','Perimeter',
                 'MajorAxis','MinorAxis', 'Angle',
                 'DistanceToWound','AngleToWound'] # Always copied from SaveAllStats
        # Initialise the variables...
        self.valuesByFrame, self.centroidX, self.centroidY = [],[],[]
        self.xmin, self.xmax, self.ymin, self.ymax = [],[],[],[]
        self.area, self.perimeter = [],[]
        self.major, self.minor, self.angle = [],[],[]
        self.woundDistance, self.woundAngle = [],[]
        
        vars = [self.valuesByFrame, self.centroidX, self.centroidY,
                self.xmin, self.xmax, self.ymin, self.ymax,
                self.area, self.perimeter,
                self.major, self.minor, self.angle,
                self.woundDistance, self.woundAngle]
        
        for i,name in enumerate(names):
            nameCsv = os.path.join(directory,name+'.csv')
            if os.path.exists(nameCsv):
                fid=open(os.path.join(directory,name+'.csv'),'r')
                data = fid.read()
                fid.close()
                
                # Remove all the '\r' first
                data = data.replace('\r','')
                
                if '\n' in data:
                    data = data.split('\n')
                else:
                    print 'File must have either \r\n (Windows) or \n (Unix) file endings!'
                    continue
                if name in ['CellValuesByFrame','Xmin','Xmax','Ymin','Ymax']:
                    data = [map(IorN,d.replace(' ','').split(',')) for d in data]
                else:
                    dprint(name)
                    data = [map(ForN,d.replace(' ','').split(',')) for d in data]
                vars[i][:] = data
            else:
                print 'No file named',nameCsv
                print 'Will not load',name
            
        # Changed this to use the csv files instead...
        #data,names = ExcelHelper.excelRead(calculationsXls,flip=True)
        # Crawl the data and replace all the emtpy strings '' with None
        #for i in range(len(data)):
        #    for j in range(len(data[i])):
        #        for k in range(len(data[i][j])):
        #            if data[i][j][k]=='':
        #                data[i][j][k]=None
        #[self.valuesByFrame, self.centroidX, self.centroidY,
        # self.xmin, self.xmax, self.ymin, self.ymax,
        # self.area, self.perimeter,
        # self.major, self.minor, self.angle,
        # self.woundDistance, self.woundAngle] = data
        
        if os.path.exists(neighborsPy):
            # Load Neighbors from Neighbors.py
            #fid=open(neighborsPy,'r')
            #exec(fid.read().replace('\r',''))
            #self.neighbors = neighbors
            #fid.close()
            fid=open(neighborsPy,'U')
            Neighbors = imp.load_module('Neighbors',fid,'Neighbors.py',('.py','U',1))
            self.neighbors = Neighbors.neighbors
            fid.close()
    
    def PlotAreasAndPerimeters(self): # Must run CollectAllStats first!
        area=np.array(self.area)
        per=np.array(self.perimeter)
        
        fn=plt.gcf().number
        for v in range(area.shape[1]):
            plt.figure(4)
            plt.plot(area[:,v],label=str(v))
            plt.figure(5)
            plt.plot(per[:,v],label=str(v))
        plt.figure(4)
        plt.legend()
        plt.figure(5)
        plt.legend()
        plt.figure(fn)
                
    def TestCalculations(self):
        self.UpdateValuesList()
        print self.valList
        print self.GetBoundingRectangles()
        print self.CalculateCentroids()
        fn=plt.gcf().number
        plt.figure(2)
        for [[xm,xM],[ym,yM]] in self.GetBoundingRectangles():
            plt.plot([ym,ym,yM,yM,ym],[xm,xM,xM,xm,xm],'k-')
        plt.figure(fn)
        
        print self.CalculateArea()
        print self.CalculatePerimeter()
        print self.CalculateBestFitEllipse()
        #print self.CalculateWoundDistance()
        print self.CalculateNeighbors()
    def CreateBinsWDistance(self,frameToCompare,woundVals,binSize=30,initR=0):
        if not hasattr(self,'woundDistance'):
            print 'Error! Need have collected stats for the data before running CreateBinsWDistance!'
            return
        d=np.array(self.woundDistance)
        # Exchange all the None's for NaN's in the wound distances...
        dB = np.array(d)
        for i in range(dB.shape[0]):
            for j in range(dB.shape[1]):
                if dB[i,j]==None:
                   dB[i,j]=np.nan
        # Sort the wound distances based on the distances at a particular time
        dBm = np.mean(dB,0)
        d_comp = d[frameToCompare] # pick a frame to use in sorting
        l = range(len(dBm))
        l.sort(cmpGen(lambda x: d_comp[x]))
        # create bins of equal size...
        # bins will each contain an index list
        nbins=int( np.ceil( (d_comp.max()-initR) / binSize ) )
        bins=[[] for i in range(nbins)]
        for i in range(len(d_comp)):
            if d_comp[i]!=None and i+2 not in woundVals: # The wound...
                ind = int( np.floor( (d_comp[i]-initR) / binSize ) )
                bins[ind].append(i)
        while bins[-1]==[]:
            del(bins[-1])
        return bins
    def CreateBinsWNeighbors(self,frameToCompare,woundVals):
        if not hasattr(self,'neighbors'):
            print 'Error! Need have collected stats for the data before running CreateBinsWNeighbors!'
            return
        nei = self.neighbors[frameToCompare]
        bins = []
        wm2 = [w-2 for w in woundVals]
        dprint(wm2)
        while True:
            bins.append([])
            iter = (wm2 if len(bins)<2 else bins[-2])
            dprint(iter)
            for w in iter:
                dprint(nei[w])
                if nei[w]!=None:
                    for v in nei[w]:
                        if v>1: # skip background
                            bins[-1].append(v-2)
            bins[-1]=list(set(bins[-1]).difference(*( [wm2]+map(set,bins[:-1])) )) # remove duplicates...
            if bins[-1]==[]:
                break
        while bins[-1]==[]:
            del(bins[-1])
        return bins
    def GetWoundNeighborContourLengths(self,index,woundVals,printPerimeterCheck=False):
        """Find the boundary lengths of the wound and it's neighboring cells"""
        # Warning! No index checking...
        newWVal = 1e6 # I don't think anyone will ever create a million cells...
        wi = np.array(self.watershed[index],dtype=np.int)
        for v in woundVals:
            wi[np.where(wi==v)] = newWVal
        # get contour lengths and contour length values for the wound-neighbor interfaces...
        WNv,WNl = ImageContour.GetPerimeterByNeighborVal(wi,newWVal)
        if printPerimeterCheck:
            totalPer = ImageContour.GetIJPerimeter(wi,newWVal)
            print 'sum of contours',sumN(WNl)
            print 'total perimeter',totalPer
        
        NWl = []
        NWv = []
        NNl = []
        NNv = []
        NOl = []
        NOv = []
        for val in WNv:
            # get contour lengths and contour length values for the wound's neighbor-neighbor interfaces...
            clv,cl = ImageContour.GetPerimeterByNeighborVal(wi,val)
            for i in range(len(clv)):
                if clv[i]==newWVal:
                    NWl.append(cl[i])
                    NWv.append(val)
                elif clv[i] in WNv:
                    NNl.append(cl[i])
                    NNv.append([val,clv[i]])
                else:
                    NOl.append(cl[i])
                    NOv.append([val,clv[i]])
        
        if printPerimeterCheck:
            print 'sum of contours (B)',sumN(NWl)
            print 'sum of contours Out',sumN(NOl)

        for i in range(len(WNv)):
            # just for good measure ;)
            if WNv != NWv: # WHAT DOES THIS REALLY DO???
                print 'Error!!! Wound-Neighbor contours do not match Neighbor-Wound contours!'
                break
            WNl[i] = min(WNl[i],NWl[i])
        NWl=WNl # Just for completeness...
        
        if printPerimeterCheck:
            print 'sum of countours (C)',sumN(WNl)

        # So, WNv and NWv are equal
        # and I should average together WNl and NWl to get the best answer for the contour length...

        for i in range(len(NNv)):
            for j in range(len(NNv)):
                if i>=len(NNv) or j>=len(NNv):
                    break
                if NNv[i]==NNv[j][::-1]:
                    NNl[i] = min(NNl[i],NNl[j])
                    del(NNv[j])
                    del(NNl[j])
        
        return WNv,WNl,NNv,NNl,NOv,NOl
    
    def GetWNOContourLengthsSorted(self,woundVals):
        WNv,WNl,NNv,NNl,NOv,NOl = [],[],[],[],[],[]
        WNsortedVals,NNsortedVals,NOsortedVals = [],[],[]
        for index in range(self.length):
            dprint('index: '+str(index))
            if self.sparseList[index]!=None:
                if self.sparseList[index].nnz!=0:
                    for i in [WNv,WNl,NNv,NNl,NOv,NOl]:
                        i.append([])
                    WNv[-1],WNl[-1],NNv[-1],NNl[-1],NOv[-1],NOl[-1] = \
                                     self.GetWoundNeighborContourLengths(index,woundVals)
                    for i in WNv[-1]:
                        if i not in WNsortedVals:
                            WNsortedVals.append(i)
                    for i in NNv[-1]:
                        if (i not in NNsortedVals) and (i[::-1] not in NNsortedVals):
                            NNsortedVals.append(i)
                    for i in NOv[-1]:
                        if i not in NOsortedVals:
                            NOsortedVals.append(i)
        WNsortedVals.sort()
        NNsortedVals.sort()
        NOsortedVals.sort()
        WNl_byVal=[]
        for v in WNsortedVals:
            WNl_byVal.append([])
            for frame in range(len(WNv)):
                if v in WNv[frame]:
                    ind=WNv[frame].index(v)
                    WNl_byVal[-1].append(WNl[frame][ind])
                else:
                    WNl_byVal[-1].append(None)
        NNl_byVal=[]
        for v in NNsortedVals:
            NNl_byVal.append([])
            for frame in range(len(NNv)):
                forw = (v in NNv[frame])
                back = (v[::-1] in NNv[frame])
                if forw and back:
                    print "Both shouldn't be in the list..."
                    break
                if forw:
                    ind=NNv[frame].index(v)
                    NNl_byVal[-1].append(NNl[frame][ind])
                if back:
                    ind=NNv[frame].index(v[::-1])
                    NNl_byVal[-1].append(NNl[frame][ind])
                if not forw and not back:
                    NNl_byVal[-1].append(None)
        NOl_byVal=[]
        for v in NOsortedVals:
            NOl_byVal.append([])
            for frame in range(len(NOv)):
                if v in NOv[frame]:
                    ind=NOv[frame].index(v)
                    NOl_byVal[-1].append(NOl[frame][ind])
                else:
                    NOl_byVal[-1].append(None)
        return WNv,WNl,NNv,NNl,NOv,NOl,WNsortedVals,WNl_byVal,NNsortedVals,NNl_byVal,NOsortedVals,NOl_byVal

    def GetBinnedContours(self,woundVals,WNv,WNl,NNl,NOl,usePrint=True):
        """Return contours for Wound Perimeter (3 versions)""" + \
        """                    Length of Tangential Boundaries from Wound""" + \
        """                    Total Perimeter Around First Ring Of Cells Around Wound (2 versions)"""
        WP,WPc2,spokes,WNP,WNPc = [],[],[],[],[]
        
        # The easy way to calculate these:
        WPc = [sumN(i) for i in WNl]
        spokes = [sumN(i) for i in NNl]
        WNPc = [sumN(i) for i in NOl]
        
        # Verification ;)
        for frame in range(self.length):
            if self.sparseList[frame]==None:
                break
            elif self.sparseList[frame].nnz==0:
                break
            # get the wound perimeter the good way ;)
            newWVal = 1e6 # I don't think anyone will ever create a million cells...
            wi = np.array(self.watershed[frame],dtype=np.int)
            for v in woundVals:
                wi[np.where(wi==v)] = newWVal
            WP.append(ImageContour.GetIJPerimeter(wi,newWVal))
            WPc2.append(sumN(ImageContour.GetPerimeterByNeighborVal(wi,newWVal)[1]))
            #spokes.append(sumN(NNl[frame]))
            
            # Total P of all NNs
            totalWN_Per = sumN([self.perimeter[frame][v-2] for v in WNv[frame]]) # I think for this we really need the local vals for each frame...

            # Perimeter of second ring around wound
            WNP.append(totalWN_Per - sumN(WNl[frame]) - 2*sumN(NNl[frame]))
            #WNPc.append(sumN(NOl[frame]))
        
        if usePrint:
            print 'Wound Perimeters'
            print WP
            print 'Wound Perimeter by Contours'
            print WPc
            print 'Wound Perimeter by Contours (B version)'
            print WPc2
            print 'Spokes'
            print spokes
            print 'Perimeter of Second Ring Around Wound'
            print WNP
            print 'Perimeter of Second Ring Around Wound (B)'
            print WNPc
        
        return WP,WPc,WPc2,spokes,WNP,WNPc
    def GetManualInputs(self,d):
        if not os.path.exists(d):
            print 'ManualInputs.py not found!'
            return
        ManualInputs = os.path.join(d,"ManualInputs.py")
        if not os.path.exists(ManualInputs):
            te = wx.TextEntryDialog(None,'Enter values for the following;\nThey will then be saved in ManualInputs.py:','ManualInputs.py',
                                    'woundVals=[]\n'
                                    'timeIntervals=[]\n'
                                    'numberFramesPerSeries=[]\n'
                                    'gapIntervals=[]\n',style=wx.TE_MULTILINE)
            te.ShowModal()
            fid = open(ManualInputs,'w')
            fid.write(te.GetValue())
            fid.close()
        #fid = open(ManualInputs,'r')
        #exec(fid.read().replace('\r',''))
        #fid.close()
        fid=open(ManualInputs,'U')
        MI = imp.load_module('ManualInputs',fid,'ManualInputs.py',('.py','U',1))
        fid.close()
        try:
            MI.woundVals
            MI.timeIntervals
            MI.gapIntervals
            MI.numberFramesPerSeries
        except:
            print "All the variables (woundVals,timeIntervals,gapIntervals,numberFramesPerSeries) must be defined!"
            return
        return MI
    def RunCalculations2(self,d):
        if not hasattr(self,'centroidX'):
            print 'You must do CollectAllStats before RunCalculations2!'
            return
        MI = self.GetManualInputs(d)
        if MI==None:
            return    
        
        woundVals=MI.woundVals
        if woundVals==[]:
            print 'Error! You need to define the woundVals in ManualInputs.py!'
            return
        if MI.timeIntervals==[] or MI.numberFramesPerSeries==[]:
            print 'Error! You need to define the time information in ManualInputs.py!'
            return
        for i in range(self.length):
            if self.sparseList[i]==None:
                break
            else:
                hasWound=False
                for wv in woundVals:
                    if wv in self.GetActiveSeedValues():
                        hasWound=True
                        break
                if not hasWound:
                    print 'Error! Frame '+str(i)+' does not have any wound!'
                    print 'Error! You must have a wound region in all frames to use RunCalculations2!!!'
                    print 'Option 1: If you messed up the wound values, correct them.'
                    print 'Option 2: Go through each initialized frame and add at least a small bit of wound where it might be.'
                    print '          Then re-run RunCalculations and then try RunCalculations2.'
                    return
        # Perform Processing and output steps here:
        ###############################
        #BTW, just hard-code:
        binSize=None
        initR=None
        useNei=True
        timeAxis = CreateTimeAxis(MI.timeIntervals,MI.gapIntervals,MI.numberFramesPerSeries,len(self.centroidX))
        
        BinningFrame = 3
        if useNei:
            bins = self.CreateBinsWNeighbors(BinningFrame,woundVals)
        else:
            bins = self.CreateBinsWDistance(BinningFrame,woundVals,binSize=binSize,initR=initR)
        # etc...
        
        # Added information about which cells are in each ring
        # Save a file with the length of each ring and the values of the cells in each ring
        fid = open(os.path.join(d,'cellsPerRing.csv'),'w')
        fid.write("The number of cells in each ring by ring:\n")
        for i in range(len(bins)):
            fid.write(str(len(bins[i]))+',')
        fid.write('\n')
        fid.write("The actual cell values for each ring (each column):\n")
        # now write out each of the cell IDs...
        for i in range(len(bins)):
            fid.write(','.join([str(b+2) for b in bins[i]]))
            fid.write('\n')
        fid.close()
        
        # Hitting Error Here:
        WNv,WNl,NNv,NNl,NOv,NOl,WNsortedVals,WNl_byVal,NNsortedVals,NNl_byVal,NOsortedVals,NOl_byVal = self.GetWNOContourLengthsSorted(woundVals)
        WP,WPc,WPc2,spokes,WNP,WNPc = self.GetBinnedContours(woundVals,WNv,WNl,NNl,NOl,usePrint=False)
        
        ExcelHelper.excelWrite( [[["Time Axis"]+timeAxis[:len(WP)],
                                  ["Wound Perimeter"]+WP,
                                  ["Wound Perimeter (Countour Method)"]+WPc,
                                  ["Wound Perimeter (Countour Method 2)"]+WPc2,
                                  ["Total Radial Contour Length"]+spokes,
                                  ["Perimeter Around 1st Ring"]+WNP,
                                  ["Perimeter Around 1st Ring (Method 2)"]+WNPc]],
                                  ['Contour Lengths'],
                                  os.path.join(d,'ContourLengths.xls'),flip=True)
        
        areaBinned = GetBinnedValueForExcel(self.area,'Area',bins,timeAxis,woundVals)
        #areaBinned -> [binNumber][frameNumber][ValueWithinBin]
        # So, for i in areaBinned:
        #   Create an excel sheet with rows denote ValueWithinBin and columns denote time
        ExcelHelper.excelWrite(areaBinned,
                               sheetnames=['Summary','Wound']+[GetPlacementName(i-1)+' Ring of Cells Around Wound' for i in range(2,len(areaBinned))],
                               filename=os.path.join(d,'BinnedAreas.xls'))
        # Output a binned file for major and minor, too...
        majorBinned = GetBinnedValueForExcel(self.major,'Minor Axis',bins,timeAxis,woundVals)
        minorBinned = GetBinnedValueForExcel(self.minor,'Major Axis',bins,timeAxis,woundVals)
        ExcelHelper.excelWrite(majorBinned,
                               sheetnames=['Summary','Wound']+[GetPlacementName(i-1)+' Ring of Cells Around Wound' for i in range(2,len(areaBinned))],
                               filename=os.path.join(d,'BinnedMajor.xls'))
        ExcelHelper.excelWrite(minorBinned,
                               sheetnames=['Summary','Wound']+[GetPlacementName(i-1)+' Ring of Cells Around Wound' for i in range(2,len(areaBinned))],
                               filename=os.path.join(d,'BinnedMinor.xls'))
        
        # I'm going to refer to binned averages here (like m up above...)
        # put each quantity in each sheet
        majorMeans = [["time"]+["ring"+str(i+1) for i in range(len(bins))]]
        minorMeans = [["time"]+["ring"+str(i+1) for i in range(len(bins))]]
        aspectRatioMeans = [["time"]+["ring"+str(i+1) for i in range(len(bins))]]
        angleMeans = [["time"]+["ring"+str(i+1) for i in range(len(bins))]]
        ax_ththMeans = [["time"]+["ring"+str(i+1) for i in range(len(bins))]]
        ax_rrMeans = [["time"]+["ring"+str(i+1) for i in range(len(bins))]]
        rrdththMeans = [["time"]+["ring"+str(i+1) for i in range(len(bins))]]
        
        for frame in range(self.length):
            if self.sparseList[frame]==None:
                break
            elif self.sparseList[index].nnz==0:
                break
            
            if i>=len(self.woundCenters):
                print 'No Wound Centers After Frame',i,'!'
                break
            
            majorMeans.append([timeAxis[frame]])
            minorMeans.append([timeAxis[frame]])
            aspectRatioMeans.append([timeAxis[frame]])
            angleMeans.append([timeAxis[frame]])
            ax_ththMeans.append([timeAxis[frame]])
            ax_rrMeans.append([timeAxis[frame]])
            rrdththMeans.append([timeAxis[frame]])
            
            for b in range(len(bins)):
                woundAngles = [self.woundAngle[frame][i] for i in bins[b]]
                cx,cy = self.woundCenters[frame]
                majorBin = []
                minorBin = []
                aspectRatioBin = []
                angleBin = []
                ax_ththBin = []
                ax_rrBin = []
                rrdththBin = []
                for i,wa in enumerate(woundAngles):
                    #print 'Value',bins[b][i]+2
                    majorBin.append( self.major[frame][bins[b][i]] )
                    minorBin.append( self.minor[frame][bins[b][i]] )
                    if None in (majorBin[-1],minorBin[-1]):
                        aspectRatioBin.append(None)
                        angleBin.append(None)
                        ax_ththBin.append(None)
                        ax_rrBin.append(None)
                        rrdththBin.append(None)
                    else:
                        if minorBin[-1]==0: # Should NEVER happen
                            aspectRatioBin.append(1e6)
                        else:
                            aspectRatioBin.append(majorBin[-1]/minorBin[-1])
                        
                        angleBin.append( self.angle[frame][bins[b][i]] )
                        
                        ax_ththBin.append(None)
                        ax_rrBin.append(None)
                        ax_ththBin[-1],ax_rrBin[-1] = EllipseFitter.GetRotatedMajorMinor(majorBin[-1],minorBin[-1],angleBin[-1],wa*180/np.pi)
                        
                        if ax_ththBin[-1]==0: # Should NEVER happen
                            rrdththBin.append(1e6)
                        else:
                            rrdththBin.append(ax_rrBin[-1]/ax_ththBin[-1])
                        
                        
                    #print 'pa_angle and woundAngle:',angleBin[-1],wa*180/np.pi
                    #print 'principle axes:',majorBin[-1],minorBin[-1]
                    #print '"radial" axes',ax_rrBin[-1],ax_ththBin[-1]
                    #print 'square comparison (for good measure...)', \
                    #      majorBin[-1]**2+minorBin[-1]**2,ax_rrBin[-1]**2+ax_ththBin[-1]**2
                majorMeans[-1].append( meanN(majorBin) )
                minorMeans[-1].append( meanN(minorBin) )
                aspectRatioMeans[-1].append( meanN(aspectRatioBin) )
                angleMeans[-1].append( meanN(angleBin) )
                ax_ththMeans[-1].append( meanN(ax_ththBin) )
                ax_rrMeans[-1].append( meanN(ax_rrBin) )
                rrdththMeans[-1].append( meanN(rrdththBin) )
        
        excelOut = [ majorMeans, minorMeans, aspectRatioMeans, angleMeans,
                     ax_ththMeans, ax_rrMeans, rrdththMeans ]
        
        sheetNames = ["Major","Minor","Aspect Ratio","Angle","Irr","Iqq","Irr D Iqq"]
        ExcelHelper.excelWrite( excelOut,sheetNames,os.path.join(d,'PrincipalAxes.xls'),flip=False)
        
        oldIndex=self.index
        self.index = 0
        self.MapPlot(saveFile=os.path.join(d,'MapWithCellIDs.png'),useText=True)
        self.mapCentroids=[]
        self.mapPlot=None
        self.index = oldIndex
        self.MapPlot()
        
        ####################################
        print 'Finished'
        return
                
    def WCPlot(self,ind,woundVals,WNv,NNv):
        wi = np.array(self.watershed[ind],dtype=np.int)
        for v in woundVals:
            wi[np.where(wi==v)]=1e6
        plt.figure(10)
        plt.clf()
        plt.ylim(plt.gca().get_ylim()[::-1])
        for v in WNv[ind]:
            for c in ImageContour.GetBoundaryLine(wi,1e6,v):
                plt.plot(*(np.array( c )[:,::-1].T),marker='.')
                sleep(0.4)
                plt.draw()
        for v in NNv[ind]:
            for c in ImageContour.GetBoundaryLine(wi,v[0],v[1]):
                plt.plot(*(np.array( c )[:,::-1].T))
                sleep(0.4)
                plt.draw()
    # Later, change this to use ExendedWoundContours instead of wcList,swcList
    def PerimeterTests(self,frame,testInd,woundVals,wcList,swcList,useSorted=False):
        """perimeter check for NO"""
        [WNv,WNl,NNv,NNl,NOv,NOl] = wcList
        [WNsortedVals,WNl_byVal,NNsortedVals,NNl_byVal,NOsortedVals,NOl_byVal] = swcList
                       
        
        if not hasattr(self,'perimeter'):
            print 'Error! Must collect stats for the data before running PerimeterTests!'
            return
        
        # get the wound perimeter the good way ;)
        newWVal = 1e6 # I don't think anyone will ever create a million cells...
        wi = np.array(self.watershed[frame],dtype=np.int)
        for v in woundVals:
            wi[np.where(wi==v)] = newWVal
        WP = ImageContour.GetIJPerimeter(wi,newWVal)
        WPc = sumN(ImageContour.GetPerimeterByNeighborVal(wi,newWVal)[1])
        print 'Wound Perimeter'
        print WP
        print 'Wound Perimeter by Contours'
        print WPc
        print 'total P of all NNs'
        print sumN([self.perimeter[frame][v-2] for v in WNv[frame]])
        
        print 'sum of WNlengths, NN lengths and NO lengths'
        if useSorted: print sumN(np.array(WNl_byVal)[:,frame]) + \
                            2*sumN(np.array(NNl_byVal)[:,frame]) + \
                            sumN(np.array(NOl_byVal)[:,frame])
        else:         print sumN(WNl[frame]) + \
                            2*sumN(NNl[frame]) + \
                            sumN(NOl[frame])
        print 'Num of WN vals:'
        if useSorted: print len(WNsortedVals)
        else:         print len(WNv[frame])
        
        if useSorted: testVal = WNsortedVals[testInd]
        else:         testVal = WNv[frame][testInd]
        
        print "Total P for",testVal,"w/GP"
        print self.perimeter[frame][testVal-2]
        vals, pers = ImageContour.GetPerimeterByNeighborVal(self.watershed[frame],testVal)
        print "Total P for",testVal,"w/GPBN"
        print sumN(pers)
        print "WN part"
        
        if useSorted: s1=WNl_byVal[testInd][frame]
        else:         s1=WNl[frame][testInd]
        
        print s1
        
        if s1==None:
            print 'This value (',testVal,') is not a wound neighbor on this frame!'
        else:
            s2 = 0
            if useSorted:
                for i,vL in enumerate(NNsortedVals):
                    if testVal in vL:
                        if NNl_byVal[i][frame]!=None:
                            s2+=NNl_byVal[i][frame]
            else:
                for i,vL in enumerate(NNv[frame]):
                    if testVal in vL:
                        if NNl[frame][i]!=None:
                            s2+=NNl[frame][i]
            print 'NN part'
            print s2
            s3 = 0
            if useSorted:
                for i,vL in enumerate(NOsortedVals):
                    if testVal in vL:
                        if NOl_byVal[i][frame]!=None:
                            s3+=NOl_byVal[i][frame]
            else:
                for i,vL in enumerate(NOv[frame]):
                    if testVal in vL:
                        if NOl[frame][i]!=None:
                            s3+=NOl[frame][i]
        print "NO part"
        print s3
        print "Sum of last 3:"
        print s1+s2+s3
    def RingsPlot(self,bins,frame,woundVals,useNei=False,binSize=None,initR=0):
        # b,g,r,c,m,y,k
        colors=[[0,0,255],[0,255,0],[255,0,0],[0,255,255],[255,0,255],[255,255,0]]*10 # that ought to do it ;-P
        rgbMap4Rings = 255*np.ones(list(self.watershed[0].shape)+[3],dtype=np.uint8)
        for i,w in enumerate(woundVals):
            rgbMap4Rings[np.where(self.watershed[frame]==w)]=colors[i]
        for i in range(len(bins)):
            for j in range(len(bins[i])):
                ind=i+2
                val=bins[i][j]+2
                rgbMap4Rings[np.where(self.watershed[frame]==val)]=colors[ind-2+len(woundVals)]
        plt.imshow(rgbMap4Rings)
        a=plt.gca()
        if not useNei:
            cir = plt.Circle( (self.woundCenters[0][0],self.woundCenters[0][1]), radius=0.1,Fill=False)
            a.add_patch(cir)
            for i in range(len(bins)):
                cir = plt.Circle( (self.woundCenters[0][0],self.woundCenters[0][1]), radius=initR+i*binSize,Fill=False)
                a.add_patch(cir)
        plt.draw()
        return rgbMap4Rings

mouseModeHelpTxt = ["Move Mode Help:\n"
                    "  Left and Right Click defer to toolbar\n"
                    "  Try out the pan and zoom in the toolbar",
                    "Add/Delete Mode Help:\n"
                    "  Left Click =  Add\n"
                    "  Right Click = Delete",
                    "Lasso Mode Help:\n"
                    "  Left Click =  Draw (can also drag)\n"
                    "  Right Click = Close loop, select points",
                    "Extra Mode Help:\n"
                    "  Left Click =  Select Region\n"
                    "  Right Click = Create Extra Seed For Selected Region ",
                    "Switch Mode Help:\n"
                    "  Left Click =  Select Region\n"
                    "  Right Click = Swap With Another Region",
                    "Draw Mode Help:\n"
                    "  Left Click = Select Region\n"
                    "  Right Click = Draw!",
                    "Center Mode Help:\n"
                    "  Left Click = Select Wound Center\n"
                    "  Right Click = Select Wound Center for ALL Frames After"]

class SegmenterFrame(wx.Frame):
    def __init__(self, *args, **kwds):
        kwds["style"] = wx.DEFAULT_FRAME_STYLE
        wx.Frame.__init__(self, *args, **kwds)
        
        self.MainPanel = wx.Panel(self, -1)
        self.FilenameText = wx.TextCtrl(self.MainPanel, -1, "")
        # Needs Tooltips instead!!
        self.MouseModeRadioBox = wx.RadioBox(self.MainPanel, -1, "Mouse Mode",
                                             choices=["Move Mode (Default)",
                                                      "Add/Delete Mode",
                                                      "Lasso Mode",
                                                      "Extra Seed (Highlight) Mode",
                                                      "Switch Regions Mode",
                                                      "Draw Mode",
                                                      "Center Mode"],
                                             majorDimension=0, style=wx.RA_SPECIFY_ROWS)
        self.MouseModeHelpText = wx.StaticText(self.MainPanel, -1, mouseModeHelpTxt[3])
        self.FrameNumber = wx.SpinCtrl(self.MainPanel, -1, "", min=0, max=1)
        self.FrameNumberText = wx.StaticText(self.MainPanel, -1, "Frame Number")
        self.ToggleViewCheckbox = wx.CheckBox(self.MainPanel, -1, "Show Segmentation Overlays?")
        self.NotesTextBox = wx.TextCtrl(self.MainPanel, -1, "Type Notes Here (per frame)",
                                        style=wx.TE_PROCESS_ENTER|wx.TE_MULTILINE)
        self.ReSaveAllFramesButton = wx.Button(self.MainPanel, -1, "Re-Save All Frames")
        self.ReRunAllWatershedsButton = wx.Button(self.MainPanel, -1, "Re-Run All Watersheds")
        self.CompressSeedValuesButton = wx.Button(self.MainPanel, -1, "Compress Seed Values")
        self.AutoCenterWoundButton = wx.Button(self.MainPanel, -1, "Auto Center Wound")
        self.RunCalculationsButton = wx.Button(self.MainPanel, -1, "Run Calculations for All Frames")
        self.RunCalculations2Button = wx.Button(self.MainPanel, -1, "Run Calculations (2) for All Frames")
        self.StatusText = wx.StaticText(self.MainPanel, -1, "Status:\n    Initializing...")
        self.SummaryText = wx.StaticText(self.MainPanel, -1, '''
Summary of Key Commands:
'''+letterKeys+'''
1-9: Change size of drawing circle in Extra/Draw Mode
Delete/Backspace: Delete selected points after lasso
Shift + Delete/Backspace: Delete unselected points after lasso
Arrow Keys: Move selected seeds (after lasso)
  +     (up)                These keys are an alternative
[   ]   (left/right)        since arrow keys are having 
  "     (down)              problems on Windows.

(BTW, you have to use Ctrl in this window, but not in others)
''')

        self.SetTitle("SeedWater Segmenter")
        self.MouseModeRadioBox.SetSelection(0)
        
        self.menuBar = wx.MenuBar()
        self.FileMenu = wx.Menu()
        item = self.FileMenu.Append(-1, "O&pen Stack or Sequence (Gif/Tiff)\tCtrl-O", "Open Gif/Tiff")
        self.Bind(wx.EVT_MENU, (lambda event: self.HandleMPLandWxKeyEvents('o')), item)
        item = self.FileMenu.Append(-1, "E&xit\tCtrl-Q", "Exit demo")
        self.Bind(wx.EVT_MENU, (lambda event: self.HandleMPLandWxKeyEvents('q')), item)
        self.menuBar.Append(self.FileMenu, "&File")
        
        self.CommandMenu = wx.Menu()
        
        for keyCommand in letterKeys.split('\n'):
            key = keyCommand[0].upper()
            descr = keyCommand[3:]
            if key in 'QO':
                continue
            item = self.CommandMenu.Append(-1, '&'+key+' - '+descr+"\tCtrl-"+key, "")
            def Temp(event):
                ekey = self.CommandMenu.FindItemById(event.Id).Label[0].lower()
                self.HandleMPLandWxKeyEvents(ekey)
                print ekey
            self.Bind(wx.EVT_MENU, Temp,item)
        self.menuBar.Append(self.CommandMenu, "&Command")
            
        self.SetMenuBar(self.menuBar)
        
        mainSizer = wx.BoxSizer(wx.HORIZONTAL)
        panelVSizer = wx.BoxSizer(wx.VERTICAL)
        splitHSizer = wx.BoxSizer(wx.HORIZONTAL)
        controlsVSizer = wx.BoxSizer(wx.VERTICAL)
        frameNumberHSizer = wx.BoxSizer(wx.HORIZONTAL)
        reHSizer = wx.BoxSizer(wx.HORIZONTAL)
        compressCenterHSizer = wx.BoxSizer(wx.HORIZONTAL)
        panelVSizer.Add(self.FilenameText, 0, wx.EXPAND, 0)
        panelVSizer.Add((10,10), 0, 0, 0)
        frameNumberHSizer.Add(self.FrameNumber, 0, 0, 0)
        frameNumberHSizer.Add(self.FrameNumberText, 0, 0, 0)
        reHSizer.Add(self.ReSaveAllFramesButton, 0, 0, 0)
        reHSizer.Add((10,10), 0, 0, 0)
        reHSizer.Add(self.ReRunAllWatershedsButton, 0, 0, 0)
        compressCenterHSizer.Add(self.CompressSeedValuesButton, 0, 0, 0)
        compressCenterHSizer.Add((10,10), 0, 0, 0)
        compressCenterHSizer.Add(self.AutoCenterWoundButton, 0, 0, 0)
        
        controlsVSizer.Add(self.MouseModeRadioBox, 0, 0, 0)
        controlsVSizer.Add((10,10), 0, 0, 0)
        controlsVSizer.Add(self.MouseModeHelpText, 0, 0, 0)
        controlsVSizer.Add((10,10), 0, 0, 0)
        controlsVSizer.Add(frameNumberHSizer, 0, 0, 0)
        controlsVSizer.Add((10,10), 0, 0, 0)
        controlsVSizer.Add(self.ToggleViewCheckbox, 0, 0, 0)
        controlsVSizer.Add((10,10), 0, 0, 0)
        controlsVSizer.Add(self.NotesTextBox, 0, wx.EXPAND, 0)
        controlsVSizer.Add((10,10), 0, 0, 0)
        controlsVSizer.Add(reHSizer, 0, 0, 0)
        controlsVSizer.Add((10,10), 0, 0, 0)
        controlsVSizer.Add(compressCenterHSizer, 0, 0, 0)
        controlsVSizer.Add((10,10), 0, 0, 0)
        controlsVSizer.Add(self.RunCalculationsButton, 0, 0, 0)
        controlsVSizer.Add((10,10), 0, 0, 0)
        controlsVSizer.Add(self.RunCalculations2Button, 0, 0, 0)
        controlsVSizer.Add((10,10), 0, 0, 0)
        
        controlsVSizer.Add(self.StatusText, 0, wx.EXPAND, 0)
        
        
        splitHSizer.Add(controlsVSizer, 1, wx.EXPAND, 0)
        splitHSizer.Add((20,20), 0, 0, 0)
        splitHSizer.Add(self.SummaryText, 0, 0, 0)
        panelVSizer.Add(splitHSizer, 1, wx.EXPAND, 0)
        self.MainPanel.SetSizer(panelVSizer)
        mainSizer.Add(self.MainPanel, 1, wx.EXPAND, 0)
        
        self.SetSizer(mainSizer)
        mainSizer.Fit(self)
        self.Layout()
        
        self.shiftDown=False
        for i in [self,self.MainPanel,self.FilenameText,self.MouseModeRadioBox,
                  self.MouseModeHelpText,self.FrameNumber,self.FrameNumberText,
                  self.ToggleViewCheckbox,self.NotesTextBox,
                  self.ReSaveAllFramesButton,self.ReRunAllWatershedsButton,
                  self.CompressSeedValuesButton,self.AutoCenterWoundButton,self.RunCalculationsButton,self.RunCalculations2Button,
                  self.StatusText,self.SummaryText]:
            i.Bind(wx.EVT_KEY_DOWN, self.OnKeyDownWx)
        
        self.Bind(wx.EVT_TEXT_ENTER, self.FileNameTextCallback, self.FilenameText)
        self.Bind(wx.EVT_RADIOBOX, self.MouseModeRadioBoxCallback, self.MouseModeRadioBox)
        self.Bind(wx.EVT_SPINCTRL, self.FrameNumberCallback, self.FrameNumber)
        self.Bind(wx.EVT_CHECKBOX, self.ToggleViewCheckBoxCallback, self.ToggleViewCheckbox)
        #self.Bind(wx.EVT_TEXT, self.NotesTextBoxCallback, self.NotesTextBox)
        self.Bind(wx.EVT_BUTTON, self.ReSaveAllFramesCallback, self.ReSaveAllFramesButton)
        self.Bind(wx.EVT_BUTTON, self.ReRunAllWatershedsCallback, self.ReRunAllWatershedsButton)
        self.Bind(wx.EVT_BUTTON, self.CompressSeedValuesCallback,self.CompressSeedValuesButton)
        self.Bind(wx.EVT_BUTTON, self.AutoCenterWoundCallback,self.AutoCenterWoundButton)
        self.Bind(wx.EVT_BUTTON, self.RunCalculationsCallback, self.RunCalculationsButton)
        self.Bind(wx.EVT_BUTTON, self.RunCalculations2Callback, self.RunCalculations2Button)
        
    def OnInit(self,filename=None,setConnections=True):
        if setConnections:
            self.SetMPLKeyConnections()
            self.mouseMode='m'
            self.MouseModeHelpText.SetLabel(mouseModeHelpTxt[0])
            self.SetConnection(self.HandleMPLMouseEvents)
            
        self.SetStatus('Loading Image Data')
        self.filename=filename
        try:
            self.saveDir
        except:
            self.saveDir=os.path.split(self.filename)[0]
        
        self.FilenameText.SetValue(self.filename)
        
        g=GTL.LoadMonolithicOrSequenceSpecial(self.filename)
        
        #if GTL.GetShapeMonolithicOrSequence(self.filename)[3]: #isSequence
        #    g=GTL.LoadFileSequence(os.path.split(self.filename)[0])
        #else:
        #    g=GTL.LoadMonolithic(self.filename)
        
        self.wd=None # Release the old data before we load the new...
        self.wd = WatershedData(g)
        self.FrameNumber.SetValue(0)
        self.FrameNumber.SetRange(0,self.wd.length)
        self.FrameNumberText.SetLabel("Frame Number (last - "+str(self.wd.length-1)+")\n(last active - 0)")
        self.ToggleViewCheckbox.SetValue(self.wd.overlayVisible)
        
        #self.wd.DrawOrigImage()
        plt.figure(1); plt.imshow(self.wd.filterData[self.wd.index],cmap=plt.cm.gray)
        
        # Automatically ask to open previous seeds (just hit cancel if none...)
        self.SetStatus('Checking on Optional Watershed and Seed Data')
        d = wx.DirSelector('Optionally Pick Directory Where Seed Info is Stored (Cancel if 1st use of program)',defaultPath=self.saveDir)
        if d !=u'':
            self.SetStatus('Loading Watershed and Seed Data')
            self.saveDir=d
            self.wd.Open(d)
            
            self.UpdateFrameLabelText()
            if self.wd.notes[self.wd.index]==None:
                self.wd.notes[self.wd.index]=''
            self.NotesTextBox.SetValue(self.wd.notes[self.wd.index])
        else:
            self.HandleMPLandWxKeyEvents('g')
        self.SetStatus('Ready')
    def SetMPLKeyConnections(self):
        for i in range(1,3):
            plt.figure(i)
            plt.figure(i).canvas.mpl_connect('key_press_event',self.OnKeyDownMPL)
            plt.figure(i).canvas.mpl_connect('key_release_event',self.OnKeyReleaseMPL)
            plt.gca().format_coord = GetReportPixel(self) # Made it so we can also just pass self here...
    def SetConnection(self,f):
        plt.figure(1)
        self.connection = plt.connect('button_press_event',f)
    
    def FileNameTextCallback(self,event):
        f=event.GetValue()
        if os.path.exists(f):
            self.OnInit(f,setConnections=False)
            d = wx.DirSelector('Optionally Pick Directory Where Seed Info is Stored',defaultPath=self.saveDir)
            if d !=u'':
                self.saveDir=d
                self.wd.Open(d)
    def MouseModeRadioBoxCallback(self,event):
        if event.GetInt()==0:
            ckey='m'
        elif event.GetInt()==1:
            ckey='a'
        elif event.GetInt()==2:
            ckey='l'
        elif event.GetInt()==3:
            ckey='e'
        elif event.GetInt()==4:
            ckey='x'
        elif event.GetInt()==5:
            ckey='d'
        elif event.GetInt()==6:
            ckey='c'
        else:
            print 'Invalid Mouse Mode!!!'
        self.HandleMPLandWxKeyEvents(ckey)
    def UpdateFrameLabelText(self):
        self.FrameNumberText.SetLabel("Frame Number (last - "+str(self.wd.length-1)+")\n(last active - "+str(self.wd.GetLastActiveFrame())+")")
    def FrameNumberCallback(self,event):
        ind=event.GetInt()
        print 'FrameNumber Callback with value', ind
        self.wd.notes[self.wd.index]=self.NotesTextBox.GetValue()
        print self.NotesTextBox.GetValue()
        i2 = max(0,ind-1)
        if self.wd.sparseList[i2]!=None and ind<self.wd.length: # If frame before the one to move to is defined
            if self.wd.sparseList[ind]==None:
                self.SetStatus('Moving to New Frame '+str(ind))
            else:
                 self.SetStatus('Moving to Frame '+str(ind))
            if ind!=self.wd.index:
                self.wd.MoveToFrame(ind)
                print "I'm in here!!"
                self.UpdateFrameLabelText()
            #if ind < self.wd.index:
                #for i in range(self.wd.index-ind):
                #    self.wd.PreviousFrame(doPlots=False)
            #elif ind > self.wd.index:
                #for i in range(ind-self.wd.index):
                #    self.wd.NextFrame(doPlots=False)
            #self.wd.ColorPlot()
            #self.wd.MapPlot()
            #self.wd.framesVisited[self.wd.index]=True
        else:
            self.FrameNumber.SetValue(self.wd.index)
        if self.wd.notes[self.wd.index]==None:
            self.wd.notes[self.wd.index]=''
        self.NotesTextBox.SetValue(self.wd.notes[self.wd.index])
        self.SetStatus('Ready')
    def ToggleViewCheckBoxCallback(self,event):
        self.wd.ToggleOverlaysVisible()
    #def NotesTextBoxCallback(self,event):
    #    self.wd.notes[self.wd.index]=self.NotesTextBox.GetValue()
    def ReSaveAllFramesCallback(self,event):
        for i in range(self.wd.length):
            if self.wd.sparseList[i]!=None:
                self.wd.framesVisited[i]=True
        self.HandleMPLandWxKeyEvents('s')
    def ReRunAllWatershedsCallback(self,event):
        t = time()
        self.SetStatus('Test run of watershed (PyMorph algorithm)')
        oldAlgo = self.wd.walgorithm[self.wd.index] # make sure to put things back like we found them...
        self.wd.UpdateSeeds()
        self.wd.Watershed('PyMorph') # Ignore undo...
        t = time()-t # see how long that took...
        self.wd.ColorPlot()
        self.wd.MapPlot()
        print 'Watershed Done -- took',t,'seconds.'
        # Get the number of frames...
        numFrames=self.wd.length
        for i in range(self.wd.length):
            if self.wd.sparseList[i]==None:
                numFrames = i
                break
        
        self.SetStatus('Checking whether to run all frames')
        dlg = wx.MessageDialog(self, 'It will take approximately '+str(t*numFrames)+' seconds to run all watersheds.\n'
                                     'Are you sure you want to continue?\n',
                                     'Run Watershed on All Frames?',
                                     wx.OK | wx.ICON_INFORMATION | wx.CANCEL)
        
        self.wd.walgorithm[self.wd.index] = oldAlgo # set algorithm back (in case it was different from PyMorph...)
        
        if dlg.ShowModal() == wx.ID_OK:
            t = time()
            oldIndex = self.wd.index
            for i in range(numFrames):
                self.SetStatus('Running Watershed on Frame '+str(i)+' out of '+str(numFrames))
                self.wd.index=i
                self.wd.UpdateSeeds()
                self.wd.Watershed( algorithm = self.wd.walgorithm[self.wd.index] )
                self.wd.framesVisited[i]=True
            self.wd.index = oldIndex
            self.wd.UpdateSeeds()
            self.wd.ColorPlot()
            self.wd.MapPlot() 
            print 'All Watersheds Done -- took',time()-t,'seconds.'
        else:
            if oldAlgo!='PyMorph':
                self.SetStatus('Running Watershed on Current Frame')
                self.wd.Watershed( algorithm = oldAlgo )
                self.wd.ColorPlot()
                self.wd.MapPlot() 
                self.wd.framesVisited[self.wd.index]=True
                print 'Watershed Done'
        self.SetStatus('Ready')
    
    def CompressSeedValuesCallback(self,event):
        t = time()
        # 1. Get all seed values as a single list
        # 2. Check to see which seed values are actually used, build a list of the values that are
        # 3. Loop through all variables that involve seedVals (watershed, )
        #   * and starting with 2, move all the old values down to their new value
        #     (this should eliminate the same-seed convergence problem)
        # 4. Check that everything did right...
        # (4b. Remove the old seedValues???)
        # 5. Replot the current frame as the underlying data should have changed...
        self.wd.ColorPlot()
        self.wd.MapPlot()
        self.SetStatus('Checking whether to compress Seed Values')
        dlg = wx.MessageDialog(self, 'Compressing Seed Values removes unused seed values,\n'
                                     'but it will change the colors of some or all cells.\n'
                                     'Are you sure you want to continue?',
                                     'Compress (Remove Unused) Seed Values?',
                                     wx.OK | wx.ICON_INFORMATION | wx.CANCEL)
        if dlg.ShowModal() == wx.ID_OK:
            print 'Change all seeds to eliminate unused seed values'
            self.SetStatus('Changing all seeds to eliminate unused seed values')
            
            self.wd.CompressSeedValues()
            
            self.wd.ColorPlot()
            self.wd.MapPlot()
            
            print 'Finished Compressing the values...'
            self.SetStatus('Ready')
    
    def AutoCenterWoundCallback(self,event):
        self.SetStatus('Finding Centers for Wound Automatically...')
        print 'Asking for Wound Values'
        dlg = wx.TextEntryDialog(None, 'Which cells are wound (enter a comma-separated list of cell values)',
            'Wound Values', '#,#')
        if dlg.ShowModal() == wx.ID_OK:
            self.SetStatus('Running AutoCenterWound')
            print 'Run AutoCenterWound'
            val = dlg.GetValue()
            try:
                woundVals=map(int,val.split(','))
            except:
                print 'Invalid Entry!!'
                self.SetStatus('Ready')
                return
            
            if None not in woundVals:
                self.wd.AutoCenterWound(woundVals)
                self.wd.ColorPlot()
        self.SetStatus('Ready')
    
    def RunCalculationsCallback(self,event):
        print 'Run Calculations!'
        #self.wd.TestCalculations()
        self.SetStatus('Running Calculations')
        for index in range(self.wd.length):
            if self.wd.sparseList[index]==None:
                break
            elif self.wd.sparseList[index].nnz==0:
                break
            
            if self.wd.woundCenters[index]==None:
                s = 'You must use Center Mode to define a' + os.linesep + \
                    'center for all initialized frames!'
                wx.MessageBox(s)
                print s
                return
        self.wd.CollectAllStats()
        self.wd.PlotAreasAndPerimeters()
        self.HandleMPLandWxKeyEvents('s')
        self.SetStatus('Ready')
    def RunCalculations2Callback(self,event):
        print 'Run Calculations (2)!'
        #self.wd.TestCalculations()
        self.SetStatus('Running Calculations (2)')
        self.wd.RunCalculations2(self.saveDir)
        self.SetStatus('Ready')
    
    def SetStatus(self,text):
        # The 200 spaces ensure that any string will print completely after "Ready" (known problem on Linux)
        self.StatusText.SetLabel('Status:\n    '+text.replace(os.linesep,'\n').replace('\n','\n    ')+' '*200)
        self.Update()
    def HandleMPLMouseEvents(self,event):
        if None not in [event.xdata,event.ydata]:
            if self.mouseMode=='m':
                pass
            elif self.mouseMode=='a':
                x, y = int(event.xdata+0.5), int(event.ydata+0.5)
                if event.inaxes:
                    if event.button == 1:
                        print 'Add',x,y
                        self.wd.NewSeed([y,x])
                        #regionNum = out[int(event.ydata) , int(event.xdata)]
                        #place = np.where(inM==regionNum)
                        #wx,wy = place[0][0],place[1][0]
                    elif event.button == 3:
                        self.wd.DeleteSeedByRegion([y,x])
                    
                    plt.figure(1)#; cla()
                    self.wd.ColorPlot()#([x,y])
                    #figure(1); plot([x],[y],'ro')
            elif self.mouseMode=='e':
                x, y = int(event.xdata+0.5), int(event.ydata+0.5)
                if self.wd.watershed!=None and self.wd.woutline!=None:
                    if event.inaxes:
                        plt.figure(1)
                        if event.button == 1:
                            self.wd.UpdateSelection([y,x])
                            self.wd.ColorPlot()
                        elif event.button == 3 and len(self.wd.selectionVals)==1:
                            print 'Add Point with value',self.wd.selectionVals[0]
                            print 'At Location',[y,x]
                            self.wd.ExtraSeed([y,x],self.wd.selectionVals[0])
                            self.wd.ColorPlot()
            elif self.mouseMode=='l':
                x, y = int(event.xdata), int(event.ydata)
                if plt.figure(1).canvas.widgetlock.locked(): return
                if event.inaxes and event.button == 1:
                    self.wd.lasso = PolyLasso(event.inaxes, self.wd.LassoCallback)
                    # acquire a lock on the widget drawing
                    plt.figure(1).canvas.widgetlock(self.wd.lasso)
                    #self.wd.ColorPlot()
            elif self.mouseMode=='x':
                x, y = int(event.xdata+0.5), int(event.ydata+0.5)
                if self.wd.watershed!=None and self.wd.woutline!=None:
                    if event.inaxes:
                        plt.figure(1)
                        if event.button == 1:
                            self.wd.UpdateSelection([y,x])
                            self.wd.ColorPlot()
                        elif event.button == 3 and len(self.wd.selectionVals)==1:
                            print 'Switch all points with values:',self.wd.selectionVals[0],\
                                  'and', self.wd.watershed[self.wd.index][y,x]
                            self.wd.SwitchRegionValues(self.wd.selectionVals[0],self.wd.watershed[self.wd.index][y,x])
                            self.wd.ColorPlot()
                            self.wd.MapPlot()
            elif self.mouseMode=='d':
                x, y = int(event.xdata+0.5), int(event.ydata+0.5)
                if self.wd.watershed!=None and self.wd.woutline!=None:
                    if event.inaxes:
                        plt.figure(1)
                        if event.button == 1:
                            self.wd.UpdateSelection([y,x])
                            self.wd.previousDrawPoint=None
                            self.wd.ColorPlot()
                        elif event.button == 2 and len(self.wd.selectionVals)==1:
                            print 'middle click...'
                            print 'this should do 2 things:'
                            print ' 1: capture a second '
                            print ' 2: switch into a new sub-mode that lets me draw 2 lines'
                        #elif event.button == 3 and in new_special_mode:
                        #
                        elif event.button == 3 and len(self.wd.selectionVals)==1:
                            if self.wd.previousDrawPoint==None:
                                self.wd.SetUndoPoint()
                                print 'Start Drawing Seeds at Point',[y,x]
                                print 'With value',self.wd.selectionVals[0]
                                self.wd.ExtraSeedLine([y,x],[y,x],self.wd.selectionVals[0])
                            else: #draw a line on the image
                                self.wd.ExtraSeedLine(self.wd.previousDrawPoint,[y,x],self.wd.selectionVals[0])
                            self.wd.previousDrawPoint = [y,x]
                            
                            self.wd.ColorPlot()
                            #self.wd.MapPlot()
            elif self.mouseMode=='c':
                x, y = event.xdata, event.ydata # Do not round x and y!
                if self.wd.watershed!=None and self.wd.woutline!=None:
                    if event.inaxes:
                        plt.figure(1)
                        if event.button == 1:
                            self.wd.woundCenters[self.wd.index]=[x,y]
                        elif event.button == 3:
                            for i in range(self.wd.index,self.wd.length):
                                self.wd.woundCenters[i]=[x,y]
                            
                        self.wd.ColorPlot()
    def OnKeyDownWx(self,event):
        key = event.GetKeyCode()
        self.shiftDown = event.ShiftDown()
        if key == wx.WXK_DELETE:
            ckey='delete'
        elif key == wx.WXK_BACK:
            ckey='backspace'
        elif key == wx.WXK_UP:
            return 'up'
        elif key == wx.WXK_DOWN:
            return 'down'
        elif key == wx.WXK_LEFT:
            return 'left'
        elif key == wx.WXK_RIGHT:
            return 'right'
        elif key==wx.WXK_ADD:
            return '+'
        elif 0<=key<256:
            ckey=chr(key).lower()
        else:
            ckey=None
        
        if event.ControlDown():
            self.HandleMPLandWxKeyEvents(ckey)
        else:
            event.Skip()
    def OnKeyDownMPL(self,event):
        if event.key=='shift':
            self.shiftDown=True
        else:
            self.HandleMPLandWxKeyEvents(event.key.lower())
    def OnKeyReleaseMPL(self,event):
        if event.key=='shift':
            self.shiftDown=False
    def HandleMPLandWxKeyEvents(self,ckey): 
        if ckey=='o': # Open new file
            self.SetStatus('Checking for File to Open')
            f = wx.FileSelector(default_filename=self.filename)
            if os.path.exists(f):
                self.SetStatus('Opening File')
                self.OnInit(f,setConnections=False)
            print 'Finished Loading'
        elif ckey=='i': # Invert for watershed and MinFinder
            self.SetStatus('Inverting Image')
            self.wd.Invert() # NEEDS UNDO_RESET
            self.wd.ColorPlot()
        elif ckey=='g': # Gauss blur for MinFinder
            self.SetStatus('Checking for Gaussian')
            dlg = wx.TextEntryDialog(None, 'Sigma Value?',
                'Gauss Average Data', '4.0')
            if dlg.ShowModal() == wx.ID_OK:
                self.SetStatus('Running Gauss Filter')
                val = dlg.GetValue()
                try:
                    sigma = float(val)
                except:
                    sigma = None
                    print 'Invalid Number!!'
                
                if sigma != None:
                    self.wd.Gauss(sigma)
                    #self.wd.DrawOrigImage()
                    
                    doWater=False
                    if self.wd.sparseList[self.wd.index]==None:
                        doWater=True
                    elif self.wd.sparseList[index].nnz==0:
                        doWater=True
                    else:
                        self.SetStatus('Checking On Update Seeds')
                        dlg = wx.MessageDialog(self, 'Are you sure you want revert any manual seed changes?',
                                   'Revert Seeds Using Minima?',
                                   wx.OK | wx.ICON_INFORMATION | wx.CANCEL)
                        if dlg.ShowModal() == wx.ID_OK:
                            doWater=True
                    if doWater:
                        self.SetStatus('Running Watershed')
                        self.wd.UpdateSeeds(force=True)
                        print 'Run Watershed'
                        self.wd.Watershed()
                        self.wd.ColorPlot()
                        self.wd.MapPlot()
                        print 'Watershed Done'
            else:
                if self.wd.sparseList[self.wd.index] == None:
                    self.wd.sparseList[self.wd.index] = scipy.sparse.lil_matrix(self.wd.shape[1:],dtype=np.uint16)
                    # I guess I could change the dtype later if I need to...
                    self.wd.seedSelections[self.wd.index] = scipy.sparse.lil_matrix(self.shape[1:],dtype=np.bool)
        elif ckey=='r': # Reset all data
            self.SetStatus('Resetting All Filters')
            print 'Reset All Filters on Figures 1 and 2'
            self.wd.ResetData() # NEEDS UNDO_RESET
            #self.wd.DrawOrigImage()
            self.wd.ColorPlot()
        elif ckey=='b': # Background for watershed
            if not STEAL_B_KEY:
                # Can't Undo
                self.SetStatus('Checking for Background Subtraction')
                dlg = wx.TextEntryDialog(None, 'Sigma Value?',
                    'Background Subtraction', '64.0')
                if dlg.ShowModal() == wx.ID_OK:
                    self.SetStatus('Running Background Subtraction')
                    val = dlg.GetValue()
                    try:
                        sigma = float(val)
                    except:
                        sigma = None
                        print 'Invalid Number!!'
                    
                    if sigma != None:
                        self.wd.BGSubtract(sigma) # NEEDS UNDO_RESET
                        self.wd.ColorPlot()
            else:
                self.HandleMPLandWxKeyEvents('n')
                if self.wd.index>0:
                    print 'Copy seeds from previous frame'
                    self.SetStatus('Copying seeds from previous frame')
                    self.wd.CopySeedsFromPrevious()# NEEDS UNDO
                    self.wd.ColorPlot()
                    self.wd.MapPlot()
        elif ckey=='h': # Sharpen for watershed
            self.SetStatus('Checking for Sharpen')
            dlg = wx.TextEntryDialog(None, 'Two Comma-separated values for Sharpen Filter?',
                'Sharpen Data', '2.0,6.0')
            if dlg.ShowModal() == wx.ID_OK:
                self.SetStatus('Running Sharpen Filter')
                val = dlg.GetValue()
                try:
                    v1,v2=val.split(',')
                    s1 = float(v1)
                    s2 = float(v2)
                except:
                    s1 = None
                    s2 = None
                    print 'Invalid Number!!'
                
                if None not in [s1,s2]:
                    self.wd.Sharpen(s1,s2) # NEEDS UNDO_RESET
                    self.wd.ColorPlot()
        elif ckey=='k': # Median Filter for watershed
            self.SetStatus('Checking for Median')
            dlg = wx.TextEntryDialog(None, 'Median Filter Size?',
                'Median Filter Data', '3')
            if dlg.ShowModal() == wx.ID_OK:
                self.SetStatus('Running Median Filter')
                val = dlg.GetValue()
                try:
                    v = float(val)
                except:
                    v = None
                    print 'Invalid Number!!'
                
                if v!=None:
                    self.wd.Median(v) # NEEDS UNDO_RESET
                    self.wd.ColorPlot()
        elif ckey=='t': # Toggle color overlays for figure 2
            self.wd.ToggleOverlaysVisible()
            self.ToggleViewCheckbox.SetValue(self.wd.overlayVisible)
        elif ckey=='s': # Save
            self.SetStatus('Checking for Save')
            d = wx.DirSelector('Optionally Pick Directory Where Seed Info is Stored',
                           defaultPath=self.saveDir)
            if d!=u'':
                self.SetStatus('Saving')
                self.wd.notes[self.wd.index]=self.NotesTextBox.GetValue()
                self.saveDir=d
                self.wd.Save(d)
                print 'Finished Saving'
            else:
                print 'Save cancelled'
        elif ckey=='u': # Undo last manual seed change
            self.SetStatus('Undoing changes')
            print 'Undo'
            self.wd.Undo()
            self.wd.UpdateSeeds()
            self.wd.ColorPlot()
        elif ckey=='n': # Next frame
            ind = self.wd.index+1
            if ind<self.wd.length:
                if self.wd.sparseList[ind]==None:
                    self.SetStatus('Moving to Frame '+str(ind)+ ' (with Watershed)')
                else:
                    self.SetStatus('Moving to Frame '+str(ind))
                
                self.wd.notes[self.wd.index]=self.NotesTextBox.GetValue()
                self.wd.NextFrame()
                self.FrameNumber.SetValue(self.wd.index)
                self.UpdateFrameLabelText()
                if self.wd.notes[self.wd.index]==None:
                    self.wd.notes[self.wd.index]=''
                self.NotesTextBox.SetValue(self.wd.notes[self.wd.index])
                print 'Moved to Frame',self.wd.index
            else:
                print 'No more frames!'
        elif ckey=='f': # Force recreate from previous centroids
            # Can Undo
            if self.wd.index>0:
                if self.shiftDown:
                    print 'Copy seeds from previous frame'
                    self.SetStatus('Copying seeds from previous frame')
                    self.wd.CopySeedsFromPrevious()# NEEDS UNDO
                else:
                    print 'Remake seeds from previous frame'
                    self.SetStatus('Remaking seeds from last frame')
                    self.wd.MakeSeedsFromPrevious()# NEEDS UNDO
                self.wd.ColorPlot()
                self.wd.MapPlot()
        elif ckey=='p': # Previous frame
            if self.shiftDown:
                ind = self.wd.lastFrameVisited
            else:
                ind=self.wd.index-1
            
            if ind>=0:
                self.SetStatus('Moving to Frame '+str(ind))
                self.wd.notes[self.wd.index]=self.NotesTextBox.GetValue()
                self.wd.MoveToFrame(ind)
                self.FrameNumber.SetValue(self.wd.index)
                self.UpdateFrameLabelText()
                if self.wd.notes[self.wd.index]==None:
                    self.wd.notes[self.wd.index]=''
                self.NotesTextBox.SetValue(self.wd.notes[self.wd.index])
                print 'Moved to Frame',self.wd.index
            else:
                print 'Already at the beginning!'
        # OpenCV usage deprecated... never that great...
        # C key stolen for center/calculate mode...
        #elif ckey=='c': # Toggle cv on or off
        #    if HAS_CV:
        #        print 'Use the faster OpenCV watershed algorithm.'
        #        self.SetStatus('Run Watershed (OpenCV algorithm)')
        #        self.wd.Watershed('OpenCV') # Ignore undo...
        #        self.wd.ColorPlot()
        #        self.wd.MapPlot()
        #        print 'Watershed Done'
        #    else:
        #        print 'OpenCV is not installed'
        #        print 'Use PyMorph instead!'
        elif ckey=='w': # Run watershed
            print 'Run Watershed (PyMorph algorithm)'
            self.SetStatus('Running Watershed (PyMorph algorithm)')
            self.wd.Watershed('PyMorph') # Ignore undo...
            self.wd.ColorPlot()
            self.wd.MapPlot()
            print 'Watershed Done'
        elif ckey=='y': # Run Dummy watershed for bad frames...
            print 'Run Dummy Watershed for Bad Frames'
            self.SetStatus('Running Dummy Watershed')
            self.wd.Watershed('Dummy') # Ignore undo...
            self.wd.ColorPlot()
            self.wd.MapPlot()
        elif ckey in ['up','down','left','right','=','+',"'",'"','[',']','{','}']:
            self.SetStatus('Moving Seeds')
            u=[0,0]
            if ckey in ['up','=','+']:
                u=[-1,0]
            elif ckey in ['down',"'",'"']:
                u=[1,0]
            elif ckey in ['left','[','{']:
                u=[0,-1]
            elif ckey in ['right',']','}']:
                u=[0,1]
            ind=self.wd.index
            
            wh = np.where(self.wd.seedSelections[self.wd.index].toarray())
            oldSA = self.wd.seedArray[wh]
            sel = self.wd.seedSelections[self.wd.index].toarray()
            self.wd.seedArray[wh]=0
            sel[wh] = 0
            print 'Whr r u?',len(wh),len(u)
            wh = ( wh[0]+u[0] , wh[1]+u[1] )
            np.clip(wh[0],0,self.wd.shape[1]) # This is an ok solution, but there is probably a better one...
            np.clip(wh[1],0,self.wd.shape[2])
            self.wd.seedArray[wh] = oldSA
            sel[wh] = 1
            self.wd.seedSelections[self.wd.index] =scipy.sparse.lil_matrix(sel,dtype=np.bool)
            
            self.wd.sparseList[self.wd.index] = scipy.sparse.lil_matrix(self.wd.seedArray,dtype=np.uint16)
            
            # I broke this old behavior for nudge for now...
            # I think the right way to do it is basically a cut and paste for an overlay array...or something...
            #for i in range(len(self.wd.seedList[ind])):
            #    if self.wd.seedSelections[ind][i]:
            #        self.wd.seedList[ind][i][0] = self.wd.seedList[ind][i][0]+u[0]
            #        self.wd.seedList[ind][i][1] = self.wd.seedList[ind][i][1]+u[1]
            
            self.wd.UpdateSeeds()
            self.wd.ColorPlot()
        elif ckey in ['delete','backspace']: # Delete After Lasso
            # Need to implement the modes like point_mode
            # and region_select_mode
            # then d should delete the seeds, etc...
            # maybe m for merge, etc...
            # Should I also apply changes directly to the watershed map if possible?
            # I think that could get really confusing really quickly...
            
            if self.wd.point_mode:
                if self.shiftDown:
                    print 'Deleting Outside Lasso Now...'
                    self.SetStatus('Deleting Unselected Points')
                    print self.wd.sparseList[self.wd.index].nnz
                    print self.wd.oldSparseList.nnz
                    self.wd.DeleteSelectedSeeds(invertSelections=True) # NEEDS UNDO
                    print self.wd.sparseList[self.wd.index].nnz
                    print self.wd.oldSparseList.nnz
                    self.wd.ColorPlot()
                    print 'Finished Deleting'
                else:
                    print 'Deleting Now...'
                    self.SetStatus('Deleting Selected Points')
                    print self.wd.sparseList[self.wd.index].nnz
                    print self.wd.oldSparseList.nnz
                    self.wd.DeleteSelectedSeeds() # NEEDS UNDO
                    print self.wd.sparseList[self.wd.index].nnz
                    print self.wd.oldSparseList.nnz
                    self.wd.ColorPlot()
                    print 'Finished Deleting'
        elif ckey=='m':
            if self.wd.watershed!=None and self.wd.woutline!=None:
                print 'Change to Move Mode (Default)'
                self.mouseMode='m'
                self.MouseModeRadioBox.SetSelection(0)
                self.MouseModeHelpText.SetLabel(mouseModeHelpTxt[0])
                self.wd.showWoundCenters=False
        elif ckey=='a': # Activate Add/Delete (by region) mode
            if self.wd.watershed!=None and self.wd.woutline!=None:
                print 'Change to Add/Delete Mode'
                self.mouseMode='a'
                self.MouseModeRadioBox.SetSelection(1)
                self.MouseModeHelpText.SetLabel(mouseModeHelpTxt[1])
                self.wd.showWoundCenters=False
        elif ckey in ['1','2','3','4','5','6','7','8','9','0']:
            self.wd.pointSize = int(ckey)
        elif ckey=='l': # Activate Extra Seed mode
            print 'Change to Lasso Mode'
            self.mouseMode='l'
            self.MouseModeRadioBox.SetSelection(2)
            self.MouseModeHelpText.SetLabel(mouseModeHelpTxt[2])
            self.wd.point_mode=True
            self.wd.showWoundCenters=False
            self.wd.ColorPlot()
        elif ckey=='e': # Activate Extra Seed mode
            if self.wd.watershed!=None and self.wd.woutline!=None:
                print 'Change to Extras (Highlight) Mode'
                self.mouseMode='e'
                self.MouseModeRadioBox.SetSelection(3)
                self.MouseModeHelpText.SetLabel(mouseModeHelpTxt[3])
                self.wd.point_mode=False
                self.wd.showWoundCenters=False
                self.wd.ColorPlot()
        elif ckey=='x': # Activate Switch mode
            if self.wd.watershed!=None and self.wd.woutline!=None:
                print 'Change to Switch Mode'
                self.mouseMode='x'
                self.MouseModeRadioBox.SetSelection(4)
                self.MouseModeHelpText.SetLabel(mouseModeHelpTxt[4])
                self.wd.point_mode=False
                self.wd.showWoundCenters=False
                self.wd.ColorPlot()
        elif ckey=='d': # Draw Mode
            if self.wd.watershed!=None and self.wd.woutline!=None:
                print 'Change to Draw Mode'
                self.mouseMode='d'
                self.MouseModeRadioBox.SetSelection(5)
                self.MouseModeHelpText.SetLabel(mouseModeHelpTxt[5])
                self.wd.point_mode=False
                self.wd.previousDrawPoint=None
                self.wd.showWoundCenters=False
                self.wd.ColorPlot()
        elif ckey=='c':
            print 'Change to Center Select Mode--For Postprocessing'
            self.mouseMode='c'
            self.MouseModeRadioBox.SetSelection(6)
            self.MouseModeHelpText.SetLabel(mouseModeHelpTxt[6])
            self.wd.showWoundCenters=True
            self.wd.ColorPlot()
        elif ckey=='j':
            if self.shiftDown and USE_DEBUG_PRINT: # Aka, if it's a developer...
                te=wx.TextEntryDialog(None,'Inject Code... very dangerous, but fun!',style=wx.TE_MULTILINE|wx.OK|wx.CANCEL)
                if te.ShowModal()==wx.ID_OK:
                    exec(te.GetValue())
            else:
                if self.wd.point_mode:
                    if self.wd.watershed!=None and self.wd.woutline!=None:
                        print 'Join Seeds to be same type'
                        self.SetStatus('Joining Seeds')
                        self.wd.MergeSelectedSeeds() # NEEDS UNDO
                        self.wd.ColorPlot()
        elif ckey=='v':
            if not self.wd.point_mode:
                if self.wd.watershed!=None and self.wd.woutline!=None:
                    print 'Change the value of a watershed region'
                    dlg = wx.TextEntryDialog(None, 'New Value?',
                        'Pick a new value for the region.', '')
                    if dlg.ShowModal() == wx.ID_OK:
                        newVal = int(dlg.GetValue())
                        self.wd.ChangeRegionValue(self.wd.selectionVals[0],newVal) # NEEDS UNDO
                        self.wd.ColorPlot()
                        self.wd.MapPlot()
        elif ckey=='q': # Quit
            dlg = wx.MessageDialog(self, 'Are you sure you want to Quit?',
                               'Exit?',
                               wx.OK | wx.ICON_INFORMATION | wx.CANCEL
                               )
            if dlg.ShowModal() == wx.ID_OK:
                wx.Exit()
        
        if ckey not in ['d','t','i']:
            self.wd.previousDrawPoint=None
        
        self.SetStatus('Ready')
            
class SegmenterApp(wx.App):
    def OnInit(self):
        wx.InitAllImageHandlers()
        
        self.InitMPL()
        
        self.frame = SegmenterFrame(None, -1, "")
        
        self.SetTopWindow(self.frame)
        
        self.frame.Show()
        return 1
    def InitMPL(self):
        # Initialize the three figures in reverse order and disable keyboard navigation keys
        for j in range(2,0,-1):
            f=plt.figure(j)
            if j==2 and DONT_PANIC:
                plt.gca()
                plt.title("Don't Panic!!!")
            plt.draw()
            for i in f.canvas.callbacks.callbacks:
                if i=='key_press_event':
                    f.canvas.mpl_disconnect(f.canvas.callbacks.callbacks[i].keys()[0])

def InitializeMPL():
    plt.ion()

if __name__=='__main__':
    InitializeMPL()
    app = SegmenterApp(0)
    if len(sys.argv)>1:
        f = ' '.join(sys.argv[1:])
    else:
        f=wx.FileSelector()
    if os.path.exists(f):
        app.frame.OnInit(filename=f)
        app.MainLoop()
