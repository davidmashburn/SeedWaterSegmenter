#!/usr/bin/env python
''''A useful collection of helper functions for SWS'''
from copy import deepcopy

import numpy as np
import matplotlib.pyplot as plt

from cmpGen import cmpGen
import ImageContour


#from list_utils:
from np_utils import flatten,totuple,interpGen
#from np_utils:
from np_utils import limitInteriorPoints,limitInteriorPointsInterpolating

def GetValuesAroundSCPoint(watershed2d,point):
    '''Given any junction in a pixel array, get the set of unique values that surround it;
       there are always 4 pixels around a junction
       there is one more junction in each direction than there are pixels
       (but one less internal junction)
       "point" must be the index of an internal junction'''
    x,y = point
    if  0 < x <= watershed2d.shape[0] and  0 < y <= watershed2d.shape[1]:
        return tuple(np.unique( watershed2d[x-1:x+1,y-1:y+1] ).tolist())
    else:
        print "THIS POINT IS NOT INTERIOR TO THE ARRAY; THERE ARE NOT 4 POINTS AROUND IT!"
        return (None,None,None)

class SubContour:
    '''A class to hold the data for a single SubContour (basically a connected list of points)
       This was designed to replace the old and crufty cVLS (contour's values,length, and subcontour)
       These three fields map to values, adjusted_length, and points respectively (much clearer!)
       This is one sc, not a list of them (that's usually referred to as 'scList')
       '''
    points = [] # list of (x,y)'s
    numPoints = 0
    adjusted_length = 0 # length as computed by a perimeter-style algorithm
    values = (None,None) # always 2 values
    startPointValues = (None,None,None) # 3 or possibly 4 values around the start point ("triple junction")
    endPointValues = (None,None,None)   # 3 or possibly 4 values around the end point ("triple junction")
    identifier = None # only used sometimes for sorting purposes
    def __init__(self,**kwds):
        for k in kwds:
            self.__dict__[k] = kwds[k]
        if 'numPoints' not in kwds.keys():
            if 'points' in kwds.keys():
                self.numPoints = len(self.points)
    def cVLS(self): # for legacy purposes, returns a list
        return [self.values,self.adjusted_length,self.points]

class CellNetwork:
    '''Holds the critical information for a single frame in order to reconstruct any subcontour of full contour'''
    subContours = [] # list of SubContours
    orderOfSubContoursDict = {} # information about reconstructing full contours; negative values means reverse the contour when inserting
    allValues = []
    def __init__(self,**kwds):
        for k in kwds:
            self.__dict__[k] = kwds[k]
    
    def GetCvlsListAndOrdering(self): # for legacy purposes, returns a list
        cVLS_List = [ sc.cVLS for sc in subContours ]
        return cVLS_List,orderOfSubContoursDict
    
    def UpdateAllValues(self):
        self.allValues = sorted(list(set( [ v for sc in self.subContours for v in sc.values ] )))
    
    def GetAllPoints(self):
        '''Get a sorted set of all points in the subContours'''
        return sorted(list(set( [ tuple(pt) for sc in self.subContours for pt in sc.points ] )))
    
    def GetXYListAndPolyList(self,closeLoops=True):
        '''Get a list of points (xyList) and a dictionary of index lists (into xyList) with cellID keys (polyList)
           polyList contains the information that reconstructs each individual contour from points' indices
               (much like orderOfSubContoursDict does but using scs' indices instead)
           'closeLoops' determines if the first point is also appended to the end of each list to close the loop and makes plotting cleaner, but be cautious of this'''
        xyList = self.GetAllPoints()
        polyList = {}
        
        for v in self.allValues
            scPointsList = [ ( self.subContours[index].points if index>0 else # reconstruct the sc's, flipping if index is negative
                               self.subContours[-index].points[::-1] )
                            for index in orderOfSubContoursDict[v] ]
            polyList[v] = [ xyList.index(totuple(pt)) for scp in scPointsList for pt in scp[:-1] ] # skip each endpoint
            if closeLoops:
                polyList[v] = polyList[v]+[polyList[v][0]] # Tack on the first point back on at the end to close the loop
                                                           # VFMin doesn't like this format; make sure to remove this last
                                                           # point before saving to a file or passing to VFM...
            #polyList[-1][v] = removeDuplicates(polyList[-1][v])+[polyList[-1][v][0]] # Remove interior duplication... bad idea!
        
        return xyList,polyList
    
    def LimitPointsBetweenNodes(self,numInteriorPointsDict,interpolate=True):
        '''Operates IN-PLACE, so use cautiously...'''
        limIntPtsFunc = limitInteriorPointsInterpolating if interpolate else limitInteriorPoints
        for sc in self.subContours:
            sc.points = limIntPtsFunc(sc.points,numInteriorPointsDict[tuple(sc.values)])
            sc.numPoints = len(sc.points)
    
    def CleanUpEmptySubContours(self):
        '''If we deleted a bunch of contours, this reindexes everything.'''
        # First things first, make a mapping from old indices to new:
        scIndexMap = {}
        count = 0
        for i in range(len(self.subContours)):
            if self.subContours[i]!=None:
                scIndexMap[i] = count
                count+=1
                
        # Now go through and delete all the dead sc's
        self.subContours = [ sc for sc in self.subContours if sc!=None ]
        
        # Now, go in and reindex orderOfSubContoursDict
        for v in self.orderOfSubContoursDict.keys():
            orderOfSubContoursDict[v] = [ scIndexMap[i] for i in orderOfSubContoursDict[v] ]
    
    def RemoveValues(self,valuesToRemove):
        '''Remove all the values from all relevant attributes'''
        
        if 1 in valuesToRemove:
            raise ValueError("You can't eliminate the background (value=1)!")
            return
        
        # Collect a list of all SC's to be outright removed:
        scsToRemoveInternal = [ (i,sc) for i,sc in enumerate(self.subContours)
                               if len(set(sc.values).intersection(valuesToRemove))==2 ] # aka, this sc is is between two values we're removing
        scsToRemoveByBackground = [ (i,sc) for i,sc in enumerate(self.subContours)
                                   if sc.values in [(1,v) for v in valuesToRemove] ] # aka, this sc is between the background and a value to be removed
        # and remove them...
        for i,sc in (scsToRemoveInternal + scsToRemoveByBackground):
            scsToRemove[i]=None
        
        # And now, replace occurrences of valuesToRemove in the sc.values by 1 (background) instead
        for sc in self.subContours:
            if sc!=None:
                len_intersect = len(set(sc.values).intersection(valuesToRemove))
                if len_intersect==1:
                   sc.values = tuple(sorted([ (1 if v in valuesToRemove else v)
                                             for v in sc.values]))
                elif len_intersect==2:
                    print 'Now how did that happen? We just filtered those out!'
                    return
        
        # Remove the values from orderOfSubContoursDict and allValues
        orderOfSubContoursDict = { v:orderOfSubContoursDict[v] for v in orderOfSubContoursDict.keys() if v not in valuesToRemove }
        allValues = [ v for v in allValues if v not in valuesToRemove ]
        
        # Clean up and we're done!
        self.CleanUpEmptySubContours()
    
    def RemoveSubContour(self,index,useSimpleRemoval=True,leaveTinyFlipFlopContour=False):
        '''This removes a subcontour by one of 3 methods:
        
        useSimple=True:   Take any connecting contours and shift the connecting points all to the midpoint of the contour to be deleted
        useSimple=False:                    THESE TWO ARE NOT IMPLEMENTED
            leaveTinyFlipFlopContour=False: A lot like simple, except that the contour will contribute multiple points to the connecting contours
            leaveTinyFlipFlopContour=True:  A lot like above except that 2 3-junctions are created instead of a 4-junction;
                                            the two parallel contours are connected at the midpoint by a very tiny contour instead
        
        In any of the 3 cases, the last step is to delete the old sc
        
        THIS OPERATES ON DATA IN-PLACE, so be careful!'''
        
        scDel = self.subContours[index]
        
        if useSimpleRemoval:
            # Find the center point of the problematic subcontour:
            #npD2 = scDel.numPoints//2
            scDelMidpoint = interpGen(scDel.points,scDel.numPoints*0.5)
                         #( scDel.points[npD2] if scDel.numPoints%2==1 else
                         #  shallowMul(shallowAdd(scDel.points[npD2-1],scDel.points[npD2]),0.5) )
            
            # Find all the subcountours that share triple junctions with the start and/or end points of scDel:
            connectedSCsToStart = [ (i,sc)
                                   for i,sc in enumerate(self.subContours)
                                   if sc!=scDel and sc!=None and ( scDel.points[0] in [sc.points[0],sc.points[-1]] ) ]
            connectedSCsToEnd = [ (i,sc)
                                 for i,sc in enumerate(self.subContours)
                                 if sc!=scDel and sc!=None and ( scDel.points[-1] in [sc.points[0],sc.points[-1]] ) ]
            #print len(connectedSCsToStart),len(connectedSCsToEnd)
            
            for scDelPtInd,connectedSCs in ((0,connectedSCsToStart),(-1,connectedSCsToEnd)):
                for i,s in connectedSCs:
                    connPtInd = ( 0 if s.points[0]==scDel.points[scDelPtInd] else -1 ) # it has to be either the start or the end of s
                    self.subContours[i].points[connPtInd] = scDelMidpoint
        else:
            print 'NOT IMPLEMENTED'
            return
            if leaveTinyFlipFlopContour:
                pass
            else:
                pass
            ########################################################
            ## Can always try something like this if the simple skip-the-contour solution doesn't work...
            ## This is the more complex, but more flexible way to do this:
            #import ImageContour.ImageContour
            #reload(ImageContour.ImageContour)
            #ImageContour.AdjustPointsAwayFromLine = ImageContour.ImageContour.AdjustPointsAwayFromLine
            #
            #cL, cR = ImageContour.AdjustPointsAwayFromLine(np.array(scDel.points),0.2,pinch=True,usePlot=True)
            #print ind,scDelPtInd.points
            #print cL
            #print cR
            #scTmp=SWHelpers.SubContour(points=cL)
            #scTmp=SWHelpers.SubContour(points=cR)
            #del scTmp
            ########################################################
            
        for v in scDel.startPointValues + scDel.endPointValues: # Luckily, we only have to check values that were actually touching the deleted sc
            if index in self.orderOfSubContoursDict[v]:
                self.orderOfSubContoursDict[v] = [ i for i in self.orderOfSubContoursDict[v] if i!=index ]
        
        self.subContours[index] = None # This saves us from having to reindex orderOfSubContoursDict until later...
                                       # use CleanUpEmptySubContours to clean up
    
    def RemoveMultipleSubContours(indexList,useSimpleRemoval=True,leaveTinyFlipFlopContour=False):
        '''Remove a bunch of subcontours and then clean up after ourselves'''
        for i in sorted(list(set(indexList)))[::-1]:
            self.RemoveSubContour(i,useSimpleRemoval,leaveTinyFlipFlopContour)
        self.CleanUpEmptySubContours()
    
    def FindMatchesAndRemovals(other,searchInds=None): # "other" is a different CellNetwork
        '''Check all the subContours to see if they have a match in 'other' and return them
           Also check for flipped matches; in the event of a flip, both sc's are flagged for removal and these lists are also returned'''
        if not other.__class__!=self.__class__:
            raise TypeError('other must be a CellNetwork!')
            return
        
        if other==self:
            raise ValueError('other must be a different object!')
            return
        
        matchedInOther = []
        removeFromSelf = []
        removeFromOther = []
        
        if searchInds==None:
            searchInds = range(len(self.subContours))
        
        for ind in searchInds:
            sc = self.subContours[ind]
            # Look for each sc from A in B:
            matchTuples = [(i,scOther) for i,scOther in enumerate(other.subContours) if sc.values==scOther.values]
            matchInds,matches = zip(*matchTuples)
            
            if len(matches)==0:
                print 'sc in self but not in other:', sc.values, matches
                # get the values connected to the subcontour only at the corners
                opposingValues = tuple(sorted(list( set(sc.startPointValues+sc.endPointValues).difference(sc.values) )))
                # get the sc's from other that have these values as their main values (aka, sc switched to this/these)
                matchOppTuples = [ (i,scOther) for i,scOther in enumerate(other.subContours) if scOther.values==opposingValues ]
                matchOppInds,matchesOpp = zip(*matchOppTuples)
                
                if matchesOpp==[]:
                    print 'Not Recoverable!'
                else:
                    print 'Recoverable: sc in A at index',ind,sc.values,'matches to sc in B at index',matchOppInds[0],matchesOpp[0].values
                    print self.subContours[ind]
                    removeFromSelf.append(ind)            # actually DO the removals later so we don't muck up the indexing!
                    removeFromOther.append(matchOppInds[0])
                    matchedInOther.append(matchOppInds[0])
                    
                    if len(matchesOpp)>1:
                        print 'More than 1 match!'
                
            elif len(matches)>1:
                print "sc in A matches multiple sc's in B:",sca.values,matches
            else:
                matchedInOther.append(matchInds[0])
                
                # Also check all the start and end points:
                # This stuff isn't really being used right now...
                sp1, sp2 = sc.startPointValues, matches[0].startPointValues
                if sp1!=sp2:
                    print "start points don't match",sp1,sp2
                ep1, ep2 = sc.endPointValues, matches[0].endPointValues
                if ep1!=ep2:
                    print "end points don't match",ep1,ep2
        
        return matchedInOther,removeFromSelf,removeFromOther


def SubContourListfromCVLSList(cVLS_List,startPointValues_List=[],endPointValues_List=[]):
    '''Get a list of SubContour objects from an old list of cVLS's'''
    if startPointValues_List==[]:
        startPointValues_List = [[None,None,None] for c in cVLS_List]
    if endPointValues_List==[]:
        endPointValues_List = [[None,None,None] for c in cVLS_List]
    return [ SubContour(points = cvls[2],
                        # numPoints = len(cvls[2]), # happens automatically...
                        adjusted_length = cvls[1],
                        values = tuple(cvls[1]),
                        startPointValues = tuple(startPointValues_List[i]),
                        endPointValues = tuple(endPointValues_List[i]))
            for i,cvls in enumerate(cVLS_List)]

def GetCellNetwork(watershed2D,allValues=None):
    '''Basically a constructor for CellNetwork based on a watershed array'''
    if allValues==None:
        allValues = np.unique(watershed2D).tolist()
    identifier=1 # unique id for each subContour -- have to start with 1 b/c we use negatives and 0 can't be negative
    scList = []
    orderOfSCsByValue = {} # For each cellID, an ordered list of indexes to the scList that reconstruct the full contour
                           # (negative values mean the sc needs to be flipped)
    for v in allValsByFrame[frame]:
        boundingRect=ImageContour.GetBoundingRect(watershed[frame],v)
        # No longer needed: #contour,turns,vals = ImageContour.GetContour(watershed[0],v,boundingRect=boundingRect,byNeighbor=True)
        perimeterVals,perimeterList,scPointsList = ImageContour.GetPerimeterByNeighborVal(watershed[frame],v,boundingRect=boundingRect,getSubContours=True)
        numSCs=len(perimeterVals)
        scPointsListAdj = [ (np.array(scp)+[boundingRect[0][0],boundingRect[1][0]]).tolist()
                       for scp in scPointsList ] # Will need to - 0.5 to line up on an overlay
        if len(perimeterList)>0:
            orderOfSCsByValue[v] = []
            for i in range(numSCs):
                newSC = SubContour( points           = scPointsAdjList[i],
                                   # numPoints        = len(scPointsAdj[i]), # happens automatically
                                    adjusted_length  = perimeterList[i],
                                    values           = tuple(sorted([v,perimeterVals[i]])),
                                    startPointValues = GetValuesAroundSCPoint( watershed[frame], scPointsAdj[i][0] ),
                                    endPointValues   = GetValuesAroundSCPoint( watershed[frame], scPointsAdj[i][-1] ),
                                    identifier=identifier )
                matchingSCs = [ sc for sc in scList if sc.values==newSC.values ] # match any subcoutours in cVLS so far that are for the same pair of cells
                matchingSCs = [ sc for sc in matchingSCs if totuple(sc.points[::-1])==totuple(newSC.points) ] # Only keep subcoutours where the points match the reverse of the points in newSC
                                #sorted([newSC.points[0],newSC.points[-1]]) == sorted([sc.points[0],sc.points[-1]]) ] # Should only possibly find 1 match...
                if matchingSCs==[]: # This is a new subContour, not a duplicate!
                    scList.append(newSC)
                    orderOfSCsByValue[v].append(identifier)
                    identifier+=1
                else:
                    matchingSCs[0].adjusted_length = min( matchingSCs[0].adjusted_length,
                                                          newSC.adjusted_length ) # keep the minimum perimeter length...
                    orderOfSCsByValue[v].append(-matchingSCs[0].identifier) # Negative identifier means the subcountour is backwards for this cell!
    scList.sort(cmpGen(lambda x: x.values)) # was just cVLS.sort()... this works, I hope?
    IDs = [sc.identifier for sc in scList]
    for sc in scList:      # scrub the id's, probably not necessary... 
        sc.identifier=None
    
    # Reindex after sorting...
    for v in allValsByFrame[frame]:       #cmp(x,0) -> sign(x) in python
        orderOfSCsByValue[v] = [ cmp(ID,0)*IDs.index(abs(ID)) for ID in orderOfSCsByValue[v] ]
    
    return CellNetwork( subContours=scList , orderOfSCsByValue=orderOfSCsByValue , allValues=allValues )

def GetCellNetworksByFrame(watershed,allValsByFrame):
    '''Get a list of CellNetworks based on a watershed segmentation'''
    return [ GetCellNetwork(watershed[i],allValsByFrame[i])
            for i in range(len(watershed)) ]

def GetXYListAndPolyListFromCellNetworkList(cellNetworkList,closeLoops=True):
    return [ cn.GetXYListAndPolyList(closeLoops=closeLoops) for cn in cellNetworkList ]


def GetCellNetworkListWithLimitedPointsBetweenNodes(cellNetworkList,splitLength=1,fixedNumInteriorPoints=None,interpolate=True):
    '''Based on matching subcontours by value pair, this function defines a fixed number of interior points for each subcontour
       and then applies this "trimming" procedure equitably to each frame in the cellNetworkList (uses LimitPointsBetweenNodes)'''
    allValues = sorted(list(set( [ v for cn in cellNetworkList for v in cn.allValues ] )))
    allPairs = sorted(list(set( [ tuple(sc.values) for cn in cellNetworkList for sc in cn.subContours ] ))) # Value pairs...
    
    # Build the numInteriorPointsDict:
    if fixedNumInteriorPoints:
        numInteriorPointsDict = {p:fixedNumInteriorPoints for p in allPairs}
    else:
        # minLength is the number of points of the shortest subcountour between cells p[0] and p[1] from all frames
        minLength = { p : min( [ sc.numPoints
                                for cn in cellNetworkList
                                for sc in cn.subContours
                                if tuple(sc.points)==p ] )
                     for p in allPairs }
        numInteriorPointsDict = { p:(minLength[p]//splitLength) for p in allPairs }
    
    cellNetworkListNew = deepcopy(cellNetworkList) # otherwise, we'd also change the input argument scListByFrame in the outside world!
    for cn in cellNetworkListNew:
        cn.LimitPointsBetweenNodes(numInteriorPointsDict,interpolate=interpolate)
    
    return cellNetworkListNew

def GetXYListAndPolyListWithLimitedPointsBetweenNodes(cellNetworkList,splitLength=1,fixedNumInteriorPoints=None,interpolate=True):
    '''Get an list of points and a connection network from a cVLS and orderOfSCsByValueByFrame; limit points between triple junctions
       (Applies GetCellNetworkListWithLimitedPointsBetweenNodes and then GetXYListAndPolyListFromCellNetworkList)'''
    return GetXYListAndPolyListFromCellNetworkList(
             GetCellNetworkListWithLimitedPointsBetweenNodes(cellNetworkList,splitLength,fixedNumInteriorPoints,interpolate) )

def GetMatchedCellNetworksCollapsing(cnA,cnB):
    '''Make 2 simplified cell networks, making sure that there is a 1-to-1 mapping between all subcontours
       This function removes values that are not common to both networks and collapses
       pairs of subcontours that do not match but are in between the same 4 cells'''
    
    if cnA==cnB: # if we got the same object for some reason, just return 2 shallow clones
        return cnA,cnB
    
    sharedVals = sorted(list(set(cnA.allValues+cnB.allValues)))
    valsNotInA = sorted(list(set(cnB.allValues).difference(cnA.allValues)))
    valsNotInB = sorted(list(set(cnA.allValues).difference(cnB.allValues)))
    
    # Delete any values that are not in both, replacing with background...
    cnA,cnB = deepcopy(cnA),deepcopy(cnB) # Make copies so we don't modify the originals
    cnA.RemoveValues(valsNotInB)
    cnB.RemoveValues(valsNotInA)
    
    matchedInB,removeFromA_a,removeFromB_a = cnA.FindMatchesAndRemovals(cnB)
    unMatchedInB = [ i for i in range(len(cnB)) if i not in matchedInB ] # This lets us skip the indexes that already matched
    
    _,removeFromB_b,removeFromA_b = cnB.FindMatchesAndRemovals(cnA,searchInds = unMatchedInB) # FLIP
    
    cnA.RemoveMultipleSubContours(removeFromA_a + removeFromA_b)
    cnB.RemoveMultipleSubContours(removeFromB_a + removeFromB_b)
    
    return cnA,cnB

def GetMatchedCellNetworksCollapsingWithLimitedPoints(cnA,cnB,splitLength=1,fixedNumInteriorPoints=None,interpolate=True):
    '''Make 2 simplified cell networks, making sure that there is a 1-to-1 mapping between all points; this function collapses
       pairs of subcontours that do not match but are in between the same 4 cells'''
    
    cnANew,cnBNew = GetMatchedSubContourListsCollapsing(cnA,cnB)
    cnALim,cnBLim = GetCellNetworkListWithLimitedPointsBetweenNodes([cnANew,cnBNew],splitLength,fixedNumInteriorPoints,interpolate)
    
    return cnALim,cnBLim

def SaveXYListAndPolyListToMMAFormat(xyList,polyList,filename,bumpIndsUp1=True,removeLastPoint=True):
    '''xyList: nested list of xy pairs for each time point.
       polyList: nested list of dictionaries for each time point where
                 each entry is like: {cellID: [ <indexes into xyList> ]}
       Exports a MMA compatible dataStructure also called "polyList" which looks like:
           {xyList,{{k,listOfIndicesTopointXYList}...}}
           where listOfIndicesTopointXYList is of course 1-indexed'''
    
    outputStr='polyList = {'
    for t,polyDict in enumerate(polyList):
        outputStr+='\n{\n'
        outputStr+=repr(xyList[t]).replace('[','{').replace(']','}').replace('(','{').replace(')','}')
        outputStr+=',\n{'
        for k in sorted(polyDict.keys()):
            inds = polyDict[k]
            if bumpIndsUp1:
                inds = [i+1 for i in inds]
            if removeLastPoint:
                inds=inds[:-1]
            outputStr+='{'+str(k)+','+repr(inds).replace('[','{').replace(']','}').replace('(','{').replace(')','}')+'}, '
        outputStr=outputStr[:-2]
        outputStr+='}\n},'
    outputStr=outputStr[:-1]+'\n}'
    open(filename,'w').write(outputStr)

def SaveCellNetworkListToMMAFormat(cellNetworkList,filename,bumpIndsUp1=True,removeLastPoint=True):
    xyAndPolyList = [ cn.GetXYListAndPolyList() for cn cellNetworkList ]
    xyList,polyList = zip(*xyAndPolyList)
    SaveXYListAndPolyListToMMAFormat(xyList,polyList,filename,bumpIndsUp1,removeLastPoint)

############################################################################################################################################################

####OLD -- a barely simpler way of doing this (aka, has no regard for ordering)
def GetSubContoursByFrame(watershed,allValsByFrame):
    '''Gets lists of lists of SubContour objects from a watershed array using ImageContour'''
    scListByFrame = []
    for frame in range(len(watershed)):
        scList = []
        for v in allValsByFrame[frame]:
            boundingRect=ImageContour.GetBoundingRect(watershed[frame],v)
            # No longer needed: #contour,turns,vals = ImageContour.GetContour(watershed[0],v,boundingRect=boundingRect,byNeighbor=True)
            perimeterVals,perimeterList,scPoints = ImageContour.GetPerimeterByNeighborVal(watershed[frame],v,boundingRect=boundingRect,getSubContours=True)
            scPointsAdj = [ (np.array(scp)+[boundingRect[0][0],boundingRect[1][0]]).tolist()
                           for scp in scPoints ] # Will need to - 0.5 to line up on an overlay
            if len(perimeterList)>0:
                scList += [ SubContour( points           = scPointsAdj[i],
                                        numPoints        = len(scPointsAdj[i]),
                                        adjusted_length  = perimeterList[i],
                                        values           = tuple(sorted([v,perimeterVals[i]])),
                                        startPointValues = GetValuesAroundSCPoint( watershed[frame], scPointsAdj[i][0] ),
                                        endPointValues   = GetValuesAroundSCPoint( watershed[frame], scPointsAdj[i][-1] ) )
                           for i in range(len(perimeterVals)) ]
        scList.sort(cmpGen(lambda x: x.values))
        for i in range(len(scList)-1,0,-1):
            # if 2 subcoutours are the same, keep only the one with the minimum length computation
            if scList[i-1].values==scList[i].values:
                scList[i-1].adjusted_length = min(scList[i-1].adjusted_length,scList[i].adjusted_length)
                del(scList[i])
        scListByFrame.append(scList)
    return scListByFrame

####OLD
def GetContourValuesLengthsAndSubContoursByFrame(watershed,allValsByFrame):
    '''Gets lists of lists of CVLS's from a watershed array -- use GetSubContoursByFrame instead!'''
    return [ [ sc.cVLS()
              for sc in scList ]
            for scList in GetSubContoursByFrame(watershed,allValsByFrame) ]

####OLD
def GetSubContoursAndOrdering(watershed2D,allValsByFrame=None):
    cellNetwork = GetCellNetwork(watershed2D,allValsByFrame)
    return cellNetwork.subContours,cellNetwork.orderOfSCsByValue

####OLD
def GetSubContoursAndOrderingByFrame(watershed,allValsByFrame):
    '''Identify every contour in the watershed segmentation for each frame as defined by each pair of cell IDs (vals).
       Also, return the order of contours that is needed to reconstruct the border surrounding each cell.'''
    cellNetworkList = GetCellNetworksByFrame(watershed,allValsByFrame)
    scListByFrame = [ cellNetworkList[i].subContours
                     for i in range(len(watershed)) ]
    orderOfSCsByValueByFrame = [ cellNetworkList[i].orderOfSubContoursDict
                                for i in range(len(watershed)) ]
    return scListByFrame,orderOfSCsByValueByFrame
    
####OLD
def GetContourValuesLengthsAndSubContoursAndOrderOfSubContoursByFrame(watershed,allValsByFrame):
    '''Identify every contour in the watershed segmentation for each frame as defined by each pair of cell IDs (vals).
       Also, return the order of contours that is needed to reconstruct the border surrounding each cell.
       DEPRECATED, use GetSubContoursAndOrderingByFrame instead!'''
    scListByFrame,orderOfSCsByValueByFrame = GetSubContoursAndOrderingByFrame(watershed,allValsByFrame)
    cVLSByFrame =  [ [ sc.cVLS()
                      for sc in scList ]
                    for scList in scListByFrame ]
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

#Old
def GetXYListAndPolyListFromCVLS(cVLS,allValsByFrame,orderOfSCsByValueByFrame):
    '''Turn a cVLS into a list of points (xyList) and a dictionary of index lists (into xyList) with cellID keys (polyList)
       polyList basically contains the information that reconstructs each individual contour from points
       Appends the first point to the end of each list to close the loop and makes plotting cleaner, but be cautious of this'''
    numFrames = len(cVLS)
    xyList = [ sorted(list(set(  tuple(pt) for c in cvlsByVal for pt in c[2]  ))) for cvlsByVal in cVLS ]
    polyList = []
    
    for t in range(numFrames):
        polyList.append({})
        for v in allValsByFrame[t]:
            subContours = [ ( cVLS[t][-index][2][::-1] if index<0 else cVLS[t][index][2] )
                           for index in orderOfSCsByValueByFrame[t][v] ] # reconstruct the sc's, flipping if index is negative
            polyList[-1][v] = [ xyList[t].index(totuple(pt)) for sc in subContours for pt in sc[:-1] ]
            polyList[-1][v] = polyList[-1][v]+[polyList[-1][v][0]] # Tack on the first point at the end to close the loop
                                                                   # VFMin doesn't like this format; make sure to remove this last point before saving to a file or passing to VFM...
            #polyList[-1][v] = removeDuplicates(polyList[-1][v])+[polyList[-1][v][0]] # Remove interior duplication...
    return xyList,polyList

####OLD
def LimitPointsBetweenSubContourNodes(scList,numInteriorPointsDict,interpolate=True):
    '''Operates IN-PLACE, so use cautiously...'''
    limIntPtsFunc = limitInteriorPointsInterpolating if interpolate else limitInteriorPoints
    for sc in scList:
        sc.points = limIntPtsFunc(sc.points,numInteriorPointsDict[tuple(sc.values)])

####OLD
def GetSubContoursByFrameWithLimitedPointsBetweenNodes(scListByFrame,allValsByFrame,splitLength=1,fixedNumInteriorPoints=None,interpolate=True):
    # IN PROGRESS...
    allValues = list(set(tuple(flatten(allValsByFrame))))
    allPairs = sorted(list(set([tuple(sc.values) for scList in scListByFrame for sc in scList]))) # Value pairs...
    
    # Build the numInteriorPointsDict:
    if fixedNumInteriorPoints:
        numInteriorPointsDict = {p:fixedNumInteriorPoints for p in allPairs}
    else:
        minLength = {}
        for p in allPairs:
            #minLength is the number of points of the shortest subcountour between cells p[0] and p[1] from all frames
            minLength[p] = min([ len(sc.adjusted_length) for scList in scListByFrame for sc in scList if tuple(sc.points)==p ])
            #                    length of subcontour
        numInteriorPointsDict = { p:(minLength[p]//splitLength) for p in allPairs }
    
    scListByFrameNew = deepcopy(scListByFrame) # otherwise, we'd also change the input argument scListByFrame in the outside world!
    for scList in scListByFrameNew:
        LimitPointsBetweenSubContourNodes(scList,numInteriorPointsDict,interpolate=interpolate)
    
    return scListByFrameNew

####OLD
def GetCVLSWithLimitedPointsBetweenNodes(cVLS,allValsByFrame,splitLength=1,fixedNumInteriorPoints=None,interpolate=True):
    '''Return a copy of the cVLS, but trim some points between triple junctions using one of several approaches.'''
    allValues = list(set(tuple(flatten(allValsByFrame))))
    allPairs = sorted(list(set([tuple(c[0]) for cVLSByFrame in cVLS for c in cVLSByFrame]))) # Value pairs...

    if fixedNumInteriorPoints:
        numInteriorPoints = {p:fixedNumInteriorPoints for p in allPairs}
    else:
        minLength = {}
        for p in allPairs:
            #minLength is the number of points of the shortest subcountour between cells p[0] and p[1] from all frames
            minLength[p] = min([ len(c[2]) for cVLSByFrame in cVLS for c in cVLSByFrame if tuple(c[0])==p ])
            #                    length of subcontour
        numInteriorPoints = {}
        for p in allPairs:
            numInteriorPoints[p] = minLength[p]//splitLength
    
    cVLS2 = deepcopy(cVLS) # otherwise, we'd also change the input argument cVLS in the outside world!
    limIntPtsFunc = limitInteriorPointsInterpolating if interpolate else limitInteriorPoints
    
    for cvlsByFrame in cVLS2:
        for c in cvlsByFrame:
            c[2] = limIntPtsFunc(c[2],numInteriorPoints[tuple(c[0])])
    
    return cVLS2

####OLD
def GetXYListAndPolyListWithLimitedPointsBetweenNodes_CVLS(cVLS,allValsByFrame,orderOfSCsByValueByFrame,splitLength=1,fixedNumInteriorPoints=None,interpolate=True):
    '''Get an list of points and a connection network from a cVLS and orderOfSCsByValueByFrame; limit points between triple junctions
       (Applies GetCVLSWithLimitedPointsBetweenNodes and then GetXYListAndPolyListFromCVLS)'''
    cVLS2 = GetCVLSWithLimitedPointsBetweenNodes(cVLS,allValsByFrame,splitLength,fixedNumInteriorPoints,interpolate=interpolate)
    return GetXYListAndPolyListFromCVLS(cVLS2,allValsByFrame,orderOfSCsByValueByFrame)

####OLD, never finished
def GetMatchedSubContourLists(scListRef,scList,allValsByFrame,orderOfSCsByValue,splitLength=1,fixedNumInteriorPoints=None):
    '''Make 2 simplified subcontour networks, one for scListRef and another for scList with a 1-to-1 mapping to the first'''
    ## NOT DONE! MERGE LATER??
    return simplifiedSCListRef,simplifiedSCList

####OLD
def RemoveSubContour(scList,index,useSimpleRemoval=True,leaveTinyFlipFlopContour=False):
    '''This removes a subcontour from sc network (scList) by one of 3 methods:
    
    useSimple=True:   Take any connecting contours and shift the connecting points all to the midpoint of the contour to be deleted
    useSimple=False:                    THESE TWO ARE NOT IMPLEMENTED
        leaveTinyFlipFlopContour=False: A lot like simple, except that the contour will contribute multiple points to the connecting contours
        leaveTinyFlipFlopContour=True:  A lot like above except that 2 3-junctions are created instead of a 4-junction;
                                        the two parallel contours are connected at the midpoint by a very tiny contour instead
    
    In any of the 3 cases, the last step is to delete the old sc
    
    THIS OPERATES ON DATA IN-PLACE, so be careful!'''
    
    scDel = scList[index]
    
    if useSimpleRemoval:
        # Find the center point of the problematic subcontour:
        #npD2 = scDel.numPoints//2
        scDelMidpoint = interpGen(scDel.points,scDel.numPoints*0.5)
                     #( scDel.points[npD2] if scDel.numPoints%2==1 else
                     #  shallowMul(shallowAdd(scDel.points[npD2-1],scDel.points[npD2]),0.5) )
        
        # Find all the subcountours that share triple junctions with the start and end points with:
        connectedSCsToStart = [ (i,sc) for i,sc in enumerate(scList)
                                       if sc!=scDel and (scDel.points[0] in (sc.points[0],sc.points[-1])) ]
        connectedSCsToEnd   = [ (i,sc) for i,sc in enumerate(scList)
                                       if sc!=scDel and (scDel.points[-1] in (sc.points[0],sc.points[-1])) ]
        #print len(connectedSCsToStart),len(connectedSCsToEnd)
        
        for scDelPtInd,connectedSCs in ((0,connectedSCsToStart),(-1,connectedSCsToEnd)):
            for i,s in connectedSCs:
                connPtInd = ( 0 if s.points[0]==scDel.points[scDelPtInd] else -1 ) # it has to be either the start or the end of s
                scList[i].points[connPtInd] = scDelMidpoint
    else:
        print 'NOT IMPLEMENTED'
        return
        if leaveTinyFlipFlopContour:
            pass
        else:
            pass
        ########################################################
        ## Can always try something like this if the simple skip-the-contour solution doesn't work...
        ## This is the more complex, but more flexible way to do this:
        #import ImageContour.ImageContour
        #reload(ImageContour.ImageContour)
        #ImageContour.AdjustPointsAwayFromLine = ImageContour.ImageContour.AdjustPointsAwayFromLine
        #
        #cL, cR = ImageContour.AdjustPointsAwayFromLine(np.array(scDel.points),0.2,pinch=True,usePlot=True)
        #print ind,scDelPtInd.points
        #print cL
        #print cR
        #scTmp=SWHelpers.SubContour(points=cL)
        #scTmp=SWHelpers.SubContour(points=cR)
        #del scTmp
        ########################################################
        
    del(scList[index])

####OLD
def FindMatchesAndRemovals(scListA,scListB,searchInds=None):
    '''Check all the subcontours in A to see if they have a match in B and return them
       Also check for flipped matches; in the event of a flip, both sc's are flagged for removal'''
    scsMatchedInB = []
    removeFromA = []
    removeFromB = []
    
    if searchInds==None:
        searchInds = range(len(scListA))
    
    for ind in searchInds:
        sca = scListA[ind]
        # Look for each sc from A in B:
        matchTuples = [(i,scb) for i,scb in enumerate(scListB) if sca.values==scb.values]
        matchInds,matches = zip(*matchTuples)
        
        if len(matches)==0:
            print 'sc in A but not in B:', sca.values, matches
            sp,ep = sca.startPointValues,sca.endPointValues
            # get the values connected to the subcontour only at the corners:
            opposingValues = tuple(sorted(list(set(sp).union(ep).difference(sca.values))))
            # get the sc's from B that have these values as their main values (aka, sca switched to this/these)
            matchOppTuples = [(i,scb) for i,scb in enumerate(scListB) if scb.values==opposingValues]
            matchOppInds,matchesOpp = zip(*matchOppTuples)
            
            if matchesOpp==[]:
                print 'Not Recoverable!'
            else:
                print 'Recoverable: sc in A at index',ind,sca.values,'matches to sc in B at index',matchOppInds[0],matchesOpp[0].values
                print scListA[ind]
                removeFromA.append(ind)            # actually DO the removals later so we don't muck up the indexing!
                removeFromB.append(matchOppInds[0])
                scsMatchedInB.append(matchOppInds[0])
                
                if len(matchesOpp)>1:
                    print 'More than 1 match!'
            
        elif len(matches)>1:
            print "sc in A matches multiple sc's in B:",sca.values,matches
        else:
            scsMatchedInB.append(matchInds[0])
            
            # Also check all the start and end points:
            # This stuff isn't really being used right now...
            sp1, sp2 = sca.startPointValues, matches[0].startPointValues
            if sp1!=sp2:
                print "start points don't match",sp1,sp2
            ep1, ep2 = sca.endPointValues, matches[0].endPointValues
            if ep1!=ep2:
                print "end points don't match",ep1,ep2
    
    return scsMatchedInB,removeFromA,removeFromB

####OLD
def GetMatchedSubContourListsCollapsing(scListA,scListB):
    '''Make 2 simplified subcontour networks, making sure that there is a 1-to-1 mapping between all points; this function collapses
       pairs of subcontours that do not match but are in between the same 4 cells '''
    
    if scListA==scListB: # if we got the same object for some reason, just return 2 shallow clones
        return scListA,scListB
    
    scsMatchedInB = [] # This keeps us from looping over both lists

    scsMatchedInB,removeFromA,removeFromB = FindMatchesAndRemovals(scListA,scListB)
    unMatchedInB = [i for i in range(len(scListB)) if i not in scsMatchedInB] # This lets us skip the indexes that already matched
    
    _,removeFromB_2,removeFromA_2 = FindMatchesAndRemovals(scListB,scListA,searchInds = unMatchedInB) # FLIP
    
    removeFromA = sorted(list(set(removeFromA + removeFromA_2)))
    removeFromB = sorted(list(set(removeFromB + removeFromB_2)))
    
    scListANew = deepcopy(scListA)
    scListBNew = deepcopy(scListB)
    
    for i in removeFromA[::-1]:
        RemoveSubContour(scListA,i)
    for i in removeFromB[::-1]:
        RemoveSubContour(scListB,i)
    
    return scListANew,scListBNew

####OLD
def GetMatchedSubContourListsCollapsingWithLimitedPoints(scListA,scListB,allValsA,allValsB,orderOfSCsByValue,splitLength=1,fixedNumInteriorPoints=None,interpolate=True):
    sharedVals = sorted(list(set(allValsA+allValsB)))
    valsNotInA = sorted(list(set(allValsB).difference(allValsA)))
    valsNotInB = sorted(list(set(allValsA).difference(allValsB)))
    
    # Remove cells not shared by both by deleting all the sc's with those cells' id's
    scListANew = [sc for sc in scListA if len(set(sc.values).intersection(valsNotInB))==0]
    scListBNew = [sc for sc in scListB if len(set(sc.values).intersection(valsNotInA))==0]
    
    scListANew,scListBNew = GetMatchedSubContourListsCollapsing(scListA,scListB)
    scListALim,scListBLim = GetSubContoursByFrameWithLimitedPointsBetweenNodes([scListANew,scListBNew],[sharedVals,sharedVals],splitLength,fixedNumInteriorPoints,interpolate)
    
    return scListALim,scListBLim

####OLD
def MakeMatchedCVLSFrames(scListByFrame,allValsByFrame,orderOfSCsByValueByFrame,frames=[0,1],splitLength=1,fixedNumInteriorPoints=None):
    # This is still a work in progress!
    if len(frames)!=2:
        raise ValueError('Frames needs to have 2 elements!')
    cVLSPrev,cVLSNext = [ cVLS[f] for f in frames ]
    valsPrev,valsNext = [ allValsByFrame[f] for f in frames ]
    
    # Dump cells not shared between both frames:
    cVLSPrev = [c for c in cVLSPrev if len(set(c[0]).intersection(valsNext))==2] # keep a cvls if the values are both found in the other frame
    cVLSNext = [c for c in cVLSNext if len(set(c[0]).intersection(valsPrev))==2] # ""
    
    pairsPrev = [tuple(c[0]) for c in cVLSPrev]
    pairsNext = [tuple(c[0]) for c in cVLSNext]
    
    if len(set(pairsPrev))!=len(pairsPrev):
        print 'Frame',frames[0],'has more than one subcontour between the same cells!'
        #return
    if len(set(pairsNext))!=len(pairsNext):
        print 'Frame',frames[1],'has more than one subcontour between the same cells!'
        #return
    
    differentSCs = set(pairsPrev).difference(pairsNext)
    if len(differentSCs)>0:
        print "Different subcontours! The following sc's are not in both frames",tuple(frames),":"
        print differentSCs
        #return
    
    if orderOfSCsByValueByFrame[frames[0]]!=orderOfSCsByValueByFrame[frames[1]]:
        print "Subcontours are in a different order between frames",tuple(frames),"!"
        #return
    
    
    cVLS = [cVLSPrev,cVLSNext]
    return GetCVLSWithLimitedPointsBetweenNodes(cVLS,allValsByFrame,splitLength,fixedNumInteriorPoints,interpolate=True) # always interpolate

####OLD
def MakeMatchedXYListPolyListFrames(cVLS,allValsByFrame,orderOfSCsByValueByFrame,frames=[0,1],splitLength=1,fixedNumInteriorPoints=None):
    cVLS2 = MakeMatchedCVLSFrames(cVLS,allValsByFrame,orderOfSCsByValueByFrame,frames,fixedNumInteriorPoints)
    return GetXYListAndPolyListFromCVLS(cVLS2,allValsByFrame,orderOfSCsByValueByFrame)

def ContourPlotFromImage(im,neighborPairs,colors=['b','g','r','c','m','y','k']):
    '''Plot an array as a grayscale image (im)
       and then plot the sub contours from an array (im) based on a set of pixel diffs
       Needs a precomputed set of neighbor pairs, but works WITHOUT ever using ImageContour
       Very useful for plotting specific contours an inspecing them (adjust neighborPairs)'''
    from ValueReceived import imshow_vr # external
    
    if len(colors)<len(neighborPairs): # Make sure there are enough colors!
        lenC = len(colors)
        for i in range(lenC,len(neighborPairs)):
            colors.append( colors[i%lenC] )
    
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
    '''Plot a cVLS'''
    for cvls in cVLSByFrame[frame]:
        cvls=np.array(cvls[2])
        _=plt.plot( cvls[:,0], cvls[:,1] )

####OLD
# I honestly don't even get what this is trying to do...
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

####OLD
# If you need to use the cVLS's anyway, why compute them twice??
# BTW, this works fine for static but DOES NOT give good connection from frame to frame
#   1: No connection from point to point through time
#   2: In-between points can come and go willy-nilly!!
def MakePolygonNetworkFromWaterSeg(waterArr,minSplit=30,allValsByFrame=None,cVLS=None):
    '''Go all the way from raw 2d+t array to an pointList,indexes format'''
    if allValsByFrame==None:
        allValsByFrame = [ list(set(w.flat)) for w in waterArr ]
    if cVLS==None:
        cVLS = GetContourValuesLengthsAndSubContoursByFrame(waterArr,allValsByFrame)
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
