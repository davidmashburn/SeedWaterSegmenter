
def GetSubContoursByFrame(watershed,allValsByFrame):
    '''Gets lists of lists of SubContour objects from a watershed array'''
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

def GetContourValuesLengthsAndSubContoursByFrame(watershed,allValsByFrame):
    '''Gets lists of lists of CVLS's from a watershed array -- use GetSubContoursByFrame instead!'''
    return [ [ sc.cVLS()
              for sc in scList ]
            for scList in GetSubContoursByFrame(watershed,allValsByFrame) ]

def GetSubContoursAndOrderingByFrame(watershed,allValsByFrame):
    '''Identify every contour in the watershed segmentation for each frame as defined by each pair of cell IDs (vals).
       Also, return the order of contours that is needed to reconstruct the border surrounding each cell.'''
    scListByFrame = []
    orderOfSCsByValueByFrame = [] # List of dicts
    for frame in range(len(watershed)):
        identifier=1 # unique id for each subContour -- have to start with 1 b/c we use negatives and 0 can't be negative
        scList = []
        orderOfSCsByValue = {} # For each cellID, an ordered list of indexes to the scList that reconstruct the full contour
                               # (negative values mean the sc needs to be flipped)
        for v in allValsByFrame[frame]:
            boundingRect=ImageContour.GetBoundingRect(watershed[frame],v)
            # No longer needed: #contour,turns,vals = ImageContour.GetContour(watershed[0],v,boundingRect=boundingRect,byNeighbor=True)
            perimeterVals,perimeterList,scPoints = ImageContour.GetPerimeterByNeighborVal(watershed[frame],v,boundingRect=boundingRect,getSubContours=True)
            numSCs=len(perimeterVals)
            scPointsAdj = [ (np.array(scp)+[boundingRect[0][0],boundingRect[1][0]]).tolist()
                           for scp in scPoints ] # Will need to - 0.5 to line up on an overlay
            if len(perimeterList)>0:
                orderOfSCsByValue[v] = []
                for i in range(numSCs):
                    newSC = SubContour( points           = scPointsAdj[i],
                                        numPoints        = len(scPointsAdj[i]),
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
        scListByFrame.append(scList)
        orderOfSCsByValueByFrame.append(orderOfSCsByValue)
    return scListByFrame,orderOfSCsByValueByFrame

###OLD
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

def GetCVLSWithLimitedPointsBetweenNodes(cVLS,allValsByFrame,splitLength=1,fixedNumInteriorPoints=None,interpolate=True):
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

def GetXYListAndPolyListWithLimitedPointsBetweenNodes(cVLS,allValsByFrame,orderOfSCsByValueByFrame,splitLength=1,fixedNumInteriorPoints=None,interpolate=True):
    cVLS2 = GetCVLSWithLimitedPointsBetweenNodes(cVLS,allValsByFrame,splitLength,fixedNumInteriorPoints,interpolate=interpolate)
    return GetXYListAndPolyListFromCVLS(cVLS2,allValsByFrame,orderOfSCsByValueByFrame)

def MakeMatchedCVLSFrames(cVLS,allValsByFrame,orderOfSCsByValueByFrame,frames=[0,1],splitLength=1,fixedNumInteriorPoints=None):
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

