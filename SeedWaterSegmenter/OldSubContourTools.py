#!/usr/bin/env python
'''Just the functions that are being kept for legacy purposes only...
   A number of the functions also don't work in the state they are in.'''

from copy import deepcopy

import numpy as np

import ImageContour


#from list_utils:
from np_utils import flatten,totuple,interpGen
#from np_utils:
from np_utils import limitInteriorPoints,limitInteriorPointsInterpolating


# a (barely) simpler way of doing this (aka, has no regard for ordering)
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
        scList.sort( key = lambda x: x.values )
        for i in range(len(scList)-1,0,-1):
            # if 2 subcoutours are the same, keep only the one with the minimum length computation
            if scList[i-1].values==scList[i].values:
                scList[i-1].adjusted_length = min(scList[i-1].adjusted_length,scList[i].adjusted_length)
                del(scList[i])
        scListByFrame.append(scList)
    return scListByFrame

def GetContourValuesLengthsAndSubContoursByFrame(watershed,allValsByFrame):
    '''Gets lists of lists of CVLS's from a watershed array -- use GetSubContoursByFrame instead!'''
    return [ [ sc.cVLS()
              for sc in scList ]
            for scList in GetSubContoursByFrame(watershed,allValsByFrame) ]

def GetSubContoursAndOrdering(watershed2D,allValsByFrame=None):
    cellNetwork = GetCellNetwork(watershed2D,allValsByFrame)
    return cellNetwork.subContours,cellNetwork.orderOfSCsByValue

def GetSubContoursAndOrderingByFrame(watershed,allValsByFrame):
    '''Identify every contour in the watershed segmentation for each frame as defined by each pair of cell IDs (vals).
       Also, return the order of contours that is needed to reconstruct the border surrounding each cell.'''
    cellNetworkList = GetCellNetworksByFrame(watershed,allValsByFrame)
    scListByFrame = [ cellNetworkList[i].subContours
                     for i in range(len(watershed)) ]
    orderOfSCsByValueByFrame = [ cellNetworkList[i].orderOfSubContoursDict
                                for i in range(len(watershed)) ]
    return scListByFrame,orderOfSCsByValueByFrame
    
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

def LimitPointsBetweenSubContourNodes(scList,numInteriorPointsDict,interpolate=True):
    '''Operates IN-PLACE, so use cautiously...'''
    limIntPtsFunc = limitInteriorPointsInterpolating if interpolate else limitInteriorPoints
    for sc in scList:
        sc.points = limIntPtsFunc(sc.points,numInteriorPointsDict[tuple(sc.values)])

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

def GetXYListAndPolyListWithLimitedPointsBetweenNodes_CVLS(cVLS,allValsByFrame,orderOfSCsByValueByFrame,splitLength=1,fixedNumInteriorPoints=None,interpolate=True):
    '''Get an list of points and a connection network from a cVLS and orderOfSCsByValueByFrame; limit points between triple junctions
       (Applies GetCVLSWithLimitedPointsBetweenNodes and then GetXYListAndPolyListFromCVLS)'''
    cVLS2 = GetCVLSWithLimitedPointsBetweenNodes(cVLS,allValsByFrame,splitLength,fixedNumInteriorPoints,interpolate=interpolate)
    return GetXYListAndPolyListFromCVLS(cVLS2,allValsByFrame,orderOfSCsByValueByFrame)

# never finished
def GetMatchedSubContourLists(scListRef,scList,allValsByFrame,orderOfSCsByValue,splitLength=1,fixedNumInteriorPoints=None):
    '''Make 2 simplified subcontour networks, one for scListRef and another for scList with a 1-to-1 mapping to the first'''
    ## NOT DONE! MERGE LATER??
    return simplifiedSCListRef,simplifiedSCList

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

def MakeMatchedXYListPolyListFrames(cVLS,allValsByFrame,orderOfSCsByValueByFrame,frames=[0,1],splitLength=1,fixedNumInteriorPoints=None):
    cVLS2 = MakeMatchedCVLSFrames(cVLS,allValsByFrame,orderOfSCsByValueByFrame,frames,fixedNumInteriorPoints)
    return GetXYListAndPolyListFromCVLS(cVLS2,allValsByFrame,orderOfSCsByValueByFrame)

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
