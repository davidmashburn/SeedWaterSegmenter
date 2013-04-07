#!/usr/bin/env python
''''A useful collection of helper functions for SWS'''
import os
from copy import deepcopy
import cPickle

import numpy as np
import matplotlib.pyplot as plt

from cmpGen import cmpGen
import ImageContour


#from list_utils:
from np_utils import totuple,interpGen,flatten,ziptranspose
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
       State: This class does not manage its own state past construction
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
    def cVLS(self):
        '''For legacy purposes, returns a list'''
        return [self.values,self.adjusted_length,self.points]
    
    def cornerConnections(self):
        '''Returns connection values at the endpoints that are not one of the 2 main values'''
        return tuple(sorted( set(self.startPointValues+self.endPointValues).difference(self.values) ))
    
    def plot(self,*args,**kwds):
        x = [ p[0]-0.5 for p in self.points ]
        y = [ p[1]-0.5 for p in self.points ]
        return plt.plot( y,x, *args, **kwds )
    def plotT(self,*args,**kwds):
        x = [ p[0]-0.5 for p in self.points ]
        y = [ p[1]-0.5 for p in self.points ]
        return plt.plot( x,y, *args, **kwds )

class QuadPoint:
    values = None
    point = None
    def __init__(self,values,point):
        if values.__class__!=tuple:
            raise TypeError('values must be a tuple')
        if not len(values)==4:
            raise TypeError('values must be a 4-tuple')
        self.values=values
        self.point=point

class DegenerateNetworkException(Exception):
    pass

class CellNetwork:
    '''Holds the critical information for a single frame in order to reconstruct any subcontour of full contour
       State: This object manages it's state and the state of all it's members like SubContour
              External functions, however should NOT mutate its state if it is passed as an argument;
              use a deepcopy-modify-return solution instead'''
    subContours = [] # list of SubContours
    quadPoints = [] # list of QuadPoints
    contourOrderingByValue = {} # information about reconstructing full contours; these should each be a tuple like:
                                # (<index into subContours>,<boolean that determines if this is forward (True) or backwards (False)>)
                                # Don't use this directly, use the GetContourPoints method instead
                                # no more negative values... (used to mean reverse the contour when inserting)
    allValues = []
    def __init__(self,**kwds):
        for k in kwds:
            self.__dict__[k] = kwds[k]
        if 'quadPoints' not in kwds.keys():
            if 'subContours' in kwds.keys():
                self.UpdateQuadPoints()
    
    def UpdateQuadPoints(self):
        '''Update quadPoints from subContours
           State: Changes state of quadPoints variable, but only to be more consistent'''
        quadPoints = sorted(set( [ (sc.startPointValues,tuple(sc.points[0])) for sc in self.subContours
                                                                             if len(sc.startPointValues)==4 ] +
                                 [ (sc.endPointValues,tuple(sc.points[-1])) for sc in self.subContours
                                                                            if len(sc.endPointValues)==4 ] ))
        self.quadPoints = [ QuadPoint(vals,pt) for vals,pt in quadPoints ]
    def GetContourPoints(self,v,closeLoop=True):
        '''Get Contour points around a value v
           State: Access only'''
        def reverseIfFalse(l,direction):
            return ( l if direction else l[::-1] )
        scPointsList = [ reverseIfFalse( self.subContours[index].points, direction ) # reconstruct sc's, flipping if direction is False
                        for index,direction in self.contourOrderingByValue[v] ] # for each index & direction tuple
        contourPoints = [ totuple(pt) for scp in scPointsList for pt in scp[:-1] ] # start point is not end point; assumed to be cyclical...
        if closeLoop:
            contourPoints.append(contourPoints[0]) # Tack on the first point back on at the end to close the loop
        return contourPoints
    
    def GetCvlsListAndOrdering(self):
        '''For legacy purposes, returns a list
           State: Access only'''
        cVLS_List = [ sc.cVLS() for sc in self.subContours ]
        return cVLS_List,self.contourOrderingByValue
    
    def GetCellularVoronoiDiagram(self):
        '''Creates a VFMin frame input structure (CellularVoronoiDiagram (cvd)), of the form:
           [
              [(x1,y1),(x2,y2)...],                      # xyList
              [
                [v,[<list of indexes into xyList...>]],  # polygon
                [v,[<list of indexes into xyList...>]],
                ...
              ]
           ]
           This is really, really similar to GetXYListAndPolyList except that here polyList is a tuple and not a dict...
           State: Access only'''
        xyList = self.GetAllPoints()
        indexMap = { tuple(pt):i for i,pt in enumerate(xyList) } # This should massively speed this up over calling xyList.index over and over
        polyList = []
        for v in self.allValues:
            contourList = self.GetContourPoints(v,closeLoop=False)
            #indexList = [ xyList.index(totuple(pt)) + 1 for pt in contourList ] # THIS IS 1-INDEXED! # this is really slow
            indexList = [ indexMap[totuple(pt)] + 1 for pt in contourList ] # THIS IS 1-INDEXED!
            polyList.append( [ v,indexList ] )
        return [xyList,polyList]
        
    def UpdateAllValues(self):
        '''Go through the values in all the subContours and collect a list of all of them
           State: Changes state of allValues, but only to to be more consistent'''
        self.allValues = sorted(set( [ v for sc in self.subContours for v in sc.values if v!=1 ] ))
    
    def GetAllPoints(self):
        '''Get a sorted set of all points in the subContours
           State: Access only'''
        return sorted(set( [ tuple(pt) for sc in self.subContours for pt in sc.points ] ))
    
    def GetXYListAndPolyList(self,closeLoops=True):
        '''Get a list of points (xyList) and a dictionary of index lists (into xyList) with cellID keys (polyList)
           polyList contains the information that reconstructs each individual contour from points' indices
               (much like contourOrderingByValue does using scs' indices)
           'closeLoops' determines if the first point is also appended to the end of each list to close the loop
                        and makes plotting cleaner, but be cautious of this
           State: Access only'''
        xyList = self.GetAllPoints()
        polyList = {}
        
        for v in self.allValues:
            contourPoints = self.GetContourPoints(v,closeLoop=False)
            polyList[v] = [ xyList.index(totuple(pt)) for pt in contourPoints ] # skip each endpoint
            if closeLoops:
                polyList[v] = polyList[v]+[polyList[v][0]] # Tack on the first point back on at the end to close the loop
                                                           # VFMin doesn't like this format; make sure to remove this last
                                                           # point before saving to a file or passing to VFM...
            #polyList[-1][v] = removeDuplicates(polyList[-1][v])+[polyList[-1][v][0]] # Remove interior duplication... bad idea!
        
        return xyList,polyList
    
    def CheckForDegenerateContours(self):
        '''Get a list off subcontours (indexes) that belong to degenerate contours (either points or lines)'''
        problemSCs = set()
        problemVals = set()
        for v,scInds in self.contourOrderingByValue.iteritems():
            if len(scInds)<3:
                scInds = [i[0] for i in scInds] # before, this also contained subcontour direction information
                if sum([ self.subContours[i].numPoints-1 for i in scInds ]) < 3:
                    # Aka, there are not enough points to make a triangle...
                    problemVals.add(v)
                    problemSCs.update(scInds)
        problemValuePairs = [self.subContours[i].values for i in problemSCs]
        return sorted(problemVals),sorted(problemValuePairs)
    
    def LimitPointsBetweenNodes(self,numInteriorPointsDict,interpolate=True,checkForDegeneracy=True):
        '''Operates IN-PLACE, so use cautiously...
           State: Changes subContours and numPoints; leaves quadPoints in an improper state
                  (allValues and contourOrderingByValue are unaffected) '''
        limIntPtsFunc = limitInteriorPointsInterpolating if interpolate else limitInteriorPoints
        for sc in self.subContours:
            sc.points = limIntPtsFunc(sc.points,numInteriorPointsDict[tuple(sc.values)])
            sc.numPoints = len(sc.points)
        
        if checkForDegeneracy:
            problemVals,problemValuePairs = self.CheckForDegenerateContours()
            if problemVals!=[]:
                print 'Values with degenerate contours:',problemVals
                raise DegenerateNetworkException('Degenerate contours (zero area) found between these cellIDs!'+repr(problemValuePairs))
    
    def CleanUpEmptySubContours(self):
        '''If we deleted a bunch of contours, this reindexes everything.
           State: Changes the variables subContours,numPoints,contourOrderingByValue
                  (allValues is not changed, but maybe it should be)
                  This should take most dirty states and clean them up'''
        # First things first, make a mapping from old indices to new:
        scIndexMap = {}
        count = 0
        for i in range(len(self.subContours)):
            if self.subContours[i]!=None:
                scIndexMap[i] = count
                count+=1
                
        # Now go through and delete all the dead sc's
        self.subContours = [ sc for sc in self.subContours if sc!=None ]
        
        # Now, go in and reindex contourOrderingByValue
        for v in self.contourOrderingByValue.keys():
            self.contourOrderingByValue[v] = [ (scIndexMap[i],d) for i,d in self.contourOrderingByValue[v] ]
        
        # And, last-but not least, update the quad points
        self.UpdateQuadPoints()
    
    def RemoveValues(self,valuesToRemove):
        '''Remove all the values from all relevant attributes -- until further notice, don't use this;
           just go back to the raw array and remake a network with fewer values...
           State: (could leave object in an inconsistent state)'''
        ### THIS NEEDS SERIOUS RETHINKING! When you leave segments behind, two or more may become part of a single loop:
        # The steps to remedy this:
        # ID ordered groups of sc's with like value pairs (after cell removal)
        # ID all the places that the sc's are used in contourOrderingByValue
        # join these all together in the correct order:
        #      add the points to end of the first segment, set the other to None
        # remove links to the dead sc and ensure that the remaining one is in the places it should be
        # do the standard cleanup; remove the None's from subContours and reindex conourOrderingByValue
        
        # So, the real key is ID'ing segs that are linked:
        #  * must have same value pair
        #  * must have an end-to-end linkage (end or start) = (end or start) = ...
        #  * must be following one another in the same (or reverse) order in ALL contourOrderings
        #  * connecting points have to be singly-connected; no internal triple junctions
        
        if valuesToRemove==[]:
            return # Do nothing...
        
        if (1 in valuesToRemove):
            raise ValueError("You can't eliminate the background (value=1)!")
        
        # Collect a list of all SC's to be outright removed:
        scsToRemoveInternal = [ (i,sc) for i,sc in enumerate(self.subContours)
                               if len(set(sc.values).intersection(valuesToRemove))==2 ] # aka, this sc is is between 2 values we're removing
        scsToRemoveByBackground = [ (i,sc) for i,sc in enumerate(self.subContours)
                                   if sc.values in [(1,v) for v in valuesToRemove] ] # aka, this sc is between the background and a value to be removed
        # and remove them...
        for i,sc in (scsToRemoveInternal + scsToRemoveByBackground):
            self.subContours[i]=None
        
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
                
                startNew = set(sc.startPointValues).difference(valuesToRemove)
                endNew = set(sc.endPointValues).difference(valuesToRemove)
                
                if len(startNew) < len(sc.startPointValues): startNew.add(1)
                if len(endNew) < len(sc.endPointValues):     endNew.add(1)
                
                sc.startPointValues = tuple(sorted(startNew))
                sc.endPointValues = tuple(sorted(endNew))
        
        '''
        # Now, go through and make a list of the junctions and the number of connections for each:
        junctions = [] # list of tuples; format: ( pt , [(scInd,True if startPt or False if endPt), ...] )
        for scInd,sc in enumerate(self.subContours):
            start,end = sc.points[0],sc.points[-1]
            jpts = [ pt for pt,connSCList in junctions ]
            tup = (scInd,True)
            if start in jpts:
                junctions[jpts.index(start)].append(tup)
            else:
                junctions.append( (start,[tup]) )
            jpts = [ pt for pt,connSCList in junctions ]
            tup = (scInd,False)
            if end in jpts:
                junctions[jpts.index(end)].append(tup)
            else:
                junctions.append( (end,[tup]) )
        singleJunctions = [ j for j in junctions if len(j[1])==1] # Shouldn't occur
        selfJunctions = [ j for j in junctions if len(j[1])==2 and j[1][0][0] == j[1][1][0]] # same scInd
        doubleJunctions = [ j for j in junctions if len(j[1])==2 and j[1][0][0] != j[1][1][0]] # diff scInd
        tripleJunctions = [ j for j in junctions if len(j[1])==3 ]
        quadJunctions = [ j for j in junctions if len(j[1])==4 ]
        
        doubleJunctionsByValuePair = {}
        for pt,scIndTupleList in doubleJunctions:
            scInds = [ scInd for scInd,startOrEnd in scIndTupleList ]
            
            if self.subContours[scInds[0]].values == self.subContours[scInds[1]].values:
                # this means the 2 sc's are between the same 2 values; should always be the case
                doubleJunctionsByValuePair[self.subContours[scInds[0]].values].append( scIndTupleList )
            else:
                print "Something went wrong; double junction between sc's with different values!!"
        
        for vals in doubleJunctionsByValuePair:
            djs_v = doubleJunctionsByValues[vals] # All the double junctions between pairs of values
                                                  # Aka, what should be internal points in a single subContour
            for scIndTuple0,scIndTuple1 in djs:
                pass # now somehow order the djs based on start and end points with the same scInd
        
        # update quad points from the quadJunctions instead??
        
        '''
        
        # Remove the values from contourOrderingByValue and allValues
        for v in valuesToRemove:
            del(self.contourOrderingByValue[v])
        
        self.allValues = [ v for v in self.allValues if v not in [1]+valuesToRemove ]
        
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
        
        THIS OPERATES ON DATA IN-PLACE, so be careful!
        State: Leaves everything in a temporary dirty state
               The clean up phase (CleanUpEmptySubContours) can then happen efficiently all at once (after multiple sc removals)'''
        
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
                    connPtInd = ( 0 if s.points[0]==scDel.points[scDelPtInd] else -1 ) # it has to be either the start or the end of a sc
                    self.subContours[i].points[connPtInd] = scDelMidpoint
        else:
            print 'NOT IMPLEMENTED'
            if leaveTinyFlipFlopContour:
                pass
            else:
                pass
            
            return
            
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
            
        for v in sorted(set(scDel.startPointValues + scDel.endPointValues)): # Luckily, we only have to check values that were touching the deleted sc
            if v!=1: # as usual, skip the background...
                contourIndices = [ i for i,d in self.contourOrderingByValue[v] ]
                if index in contourIndices:
                    self.contourOrderingByValue[v] = [ (i,d) for i,d in self.contourOrderingByValue[v] if i!=index ]
        
        self.subContours[index] = None # This saves us from having to reindex contourOrderingByValue until later...
                                       # use CleanUpEmptySubContours to clean up
    
    def RemoveMultipleSubContours(self,indexList,useSimpleRemoval=True,leaveTinyFlipFlopContour=False):
        '''Remove a bunch of subcontours and then clean up after ourselves
           State: Changes the state of all variables (except possibly allValues)
                  Should leave everything in a clean state'''
        for i in sorted(set(indexList))[::-1]:
            self.RemoveSubContour(i,useSimpleRemoval,leaveTinyFlipFlopContour)
        self.CleanUpEmptySubContours()
    
    def FindMatchesAndRemovals(self,other,searchInds=None): # "other" is a different CellNetwork
        '''Check all the subContours to see if they have a match in 'other' and return them
           Also check for flipped matches; in the event of a flip, both sc's are flagged for
           removal and these lists are also returned
           State: Access only'''
        if other.__class__!=self.__class__:
            raise TypeError('other must be a CellNetwork!')
        
        if other==self:
            raise ValueError('other must be a different object!')
        
        matchedInOther = []
        removeFromSelf = []
        removeFromOther = []
        notRecoverable = []
        
        if searchInds==None:
            searchInds = range(len(self.subContours))
        
        for ind in searchInds:
            sc = self.subContours[ind]
            # get the values connected to the subcontour only at the corners
            opposingValues = sc.cornerConnections()
            
            # Look for matches between this sc in A and a sc (or possibly a 4-junction) in B with various properties...
            keepSearching = True
            collapse = False
            
            # 1: look for a subContour that has the same values and the same start/end point values
            if keepSearching:
                matchTuples = [ (i,scOther)
                               for i,scOther in enumerate(other.subContours)
                               if ( sc.values==scOther.values and
                                    ( [sc.startPointValues,sc.endPointValues]==[scOther.startPointValues,scOther.endPointValues] or
                                      [sc.startPointValues,sc.endPointValues]==[scOther.endPointValues,scOther.startPointValues] ) ) ]
                matchInds,matches = ([],[]) if matchTuples==[] else zip(*matchTuples)
                if len(matches)==1:
                    # silent match
                    matchedInOther.append(matchInds[0])
                    keepSearching = False
                elif len(matches)>1:
                    print 'Not recoverable'
                    print "sc in A matches multiple sc's in B with the same start/end point values:",sc.values,matches
                    keepSearching = False
            # 2: look for a subContour that has the same values and either one or the other same start/end point values
            if keepSearching:
                matchTuples = [ (i,scOther)
                               for i,scOther in enumerate(other.subContours)
                               if ( sc.values==scOther.values and
                                    ( (sc.startPointValues in [scOther.startPointValues,scOther.endPointValues]) or
                                      (sc.endPointValues in [scOther.endPointValues,scOther.startPointValues]) ) ) ]
                matchInds,matches = ([],[]) if matchTuples==[] else zip(*matchTuples)
                if len(matches)==1:
                    # silent match
                    matchedInOther.append(matchInds[0])
                    keepSearching = False
                elif len(matches)>1:
                    print 'Not recoverable'
                    print "sc in A matches multiple sc's in B with one common endpoint:",sc.values,matches
                    keepSearching = False
            # 3: look for a subContour that has the same values without the same start/endpoints
            if keepSearching:
                matchTuples = [ (i,scOther)
                               for i,scOther in enumerate(other.subContours)
                               if sc.values==scOther.values ]
                matchInds,matches = ([],[]) if matchTuples==[] else zip(*matchTuples)
                if len(matches)==1:
                    # silent match
                    matchedInOther.append(matchInds[0])
                    keepSearching = False
                elif len(matches)>1:
                    print 'Not recoverable'
                    print "sc in A matches multiple sc's in B without common start/endpoints:",sc.values,matches
                    keepSearching = False
            
            # 4: look for a subContour that has opposing values (values = other's cornerValues and vice versa)
            if keepSearching:
                matchTuples = [ (i,scOther)
                               for i,scOther in enumerate(other.subContours)
                               if sc.values==scOther.cornerConnections() and opposingValues==scOther.values ]
                matchInds,matches = ([],[]) if matchTuples==[] else zip(*matchTuples)
                if len(matches)==1:
                    print 'Recoverable: sc in A at index',ind,sc.values,'matches to sc in B at index',matchInds[0],matches[0].values
                    matchedInOther.append(matchInds[0])
                    removeFromSelf.append(ind)        # actually DO the removals later so we don't muck up the indexing!
                    removeFromOther.append(matchInds[0])
                    keepSearching = False
                    collapse=True
                elif len(matches)>1:
                    print 'Not recoverable'
                    print "sc in A matches multiple sc's in B with flip-flopped values:",sc.values,matches
                    keepSearching = False
            
            # 5: look for a subContour that has either of the opposing values (values = other cornerValues OR vice versa)
            if keepSearching:
                matchTuples = [ (i,scOther)
                               for i,scOther in enumerate(other.subContours)
                               if sc.values==scOther.cornerConnections() or opposingValues==scOther.values ]
                matchInds,matches = ([],[]) if matchTuples==[] else zip(*matchTuples)
                if len(matches)==1:
                    print 'Recoverable: sc in A at index',ind,sc.values,'matches to sc in B at index',matchInds[0],matches[0].values
                    matchedInOther.append(matchInds[0])
                    removeFromSelf.append(ind)        # actually DO the removals later so we don't muck up the indexing!
                    removeFromOther.append(matchInds[0])
                    keepSearching = False
                    collapse=True
                elif len(matches)>1:
                    print 'Not recoverable'
                    print "sc in A matches multiple sc's in B with flip-flopped values (OR):",sc.values,matches
                    keepSearching = False
            
            # 6: look for a subContour that collapsed to a 4-junction
            if keepSearching:
                allValuesAroundSC = tuple(sorted(set(sc.startPointValues+sc.endPointValues)))
                matches = [ q for q in other.quadPoints
                              if q.values==allValuesAroundSC ] # Check to see if this sc collapsed into a 4-junction
                if len(matches)==1:
                    #print 'Recoverable: sc in A at index',ind,sc.values,'collapsed into a 4-junction with values',quadPointMatches[0].values,'at',quadPointMatches[0].point
                    removeFromSelf.append(ind)        # actually DO the removals later so we don't muck up the indexing!
                    keepSearching = False
                    collapse=True
                elif len(matches)>1:
                    print 'Not recoverable'
                    print "sc in A matches multiple 4-junctions in B:",sc.values,matches
                    keepSearching = False
            
            # We failed to find the sc in B
            if keepSearching:
                print 'Not Recoverable! No matches found'
                notRecoverable.append(sc.values)
            
            # Also check all the start and end points:
            # This stuff isn't really being used right now...
            #sp1, sp2 = sc.startPointValues, matches[0].startPointValues
            #if sp1!=sp2:
            #    print "start points don't match",sp1,sp2
            #ep1, ep2 = sc.endPointValues, matches[0].endPointValues
            #if ep1!=ep2:
            #    print "end points don't match",ep1,ep2
        
        return matchedInOther,removeFromSelf,removeFromOther,notRecoverable
    
    def scPlot(self,*args,**kwds):
        '''Plot the subContours in a way that can be overlaid on an imshow
           State: Access only'''
        for sc in self.subContours:
            _=sc.plot(*args,**kwds)
    def cellPlot(self,*args,**kwds):
        '''Plot the full contours in a way that can be overlaid on an imshow
           State: Access only'''
        contourPoints = { v:self.GetContourPoints(v) for v in self.contourOrderingByValue.keys() }
        for v in contourPoints.keys():
            x = [ p[0]-0.5 for p in contourPoints[v] ]
            y = [ p[1]-0.5 for p in contourPoints[v] ]
            _=plt.plot( y,x, *args, **kwds )
    def scPlotT(self,*args,**kwds):
        '''Plot the subContours in a way that can be overlaid on a transposed imshow
           (matches closely with diagramPlot in VFMLite)
           State: Access only'''
        for sc in self.subContours:
            _=sc.plotT(*args,**kwds)
    def cellPlotT(self,*args,**kwds):
        '''Plot the full contours in a way that can be overlaid on a transposed imshow
           (matches closely with diagramPlot in VFMLite)
           State: Access only'''
        contourPoints = { v:self.GetContourPoints(v) for v in self.contourOrderingByValue.keys() }
        for v in contourPoints.keys():
            x = [ p[0]-0.5 for p in contourPoints[v] ]
            y = [ p[1]-0.5 for p in contourPoints[v] ]
            _=plt.plot( x,y, *args, **kwds )
        

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

def GetCellNetwork(watershed2d,allValues=None):
    '''Basically a constructor for CellNetwork based on a watershed array'''
    if allValues==None:
        allValues = np.unique(watershed2d)
    allValues = np.array(allValues).tolist() # force numpy arrays to lists
    allValues = [ v for v in allValues if v!=1 ] # skip the background
    
    identifier=0 # unique id for each subContour
    scList = []
    contourOrderingByValue = {} # For each cellID, an ordered list of index to the scList/direction pairs that reconstruct the full contour
    for v in allValues:
        boundingRect=ImageContour.GetBoundingRect(watershed2d,v)
        # No longer needed: #contour,turns,vals = ImageContour.GetContour(watershed[0],v,boundingRect=boundingRect,byNeighbor=True)
        perimeterVals,perimeterList,scPointsList = ImageContour.GetPerimeterByNeighborVal(watershed2d,v,boundingRect=boundingRect,getSubContours=True)
        numSCs=len(perimeterVals)
        scPointsListAdj = [ (np.array(scp)+[boundingRect[0][0],boundingRect[1][0]]).tolist()
                       for scp in scPointsList ] # Will need to - 0.5 to line up on an overlay
        if len(perimeterList)>0:
            contourOrderingByValue[v] = []
            for i in range(numSCs):
                newSC = SubContour( points           = scPointsListAdj[i],
                                  # numPoints        = len(scPointsAdj[i]), # happens automatically
                                    adjusted_length  = perimeterList[i],
                                    values           = tuple(sorted([v,perimeterVals[i]])),
                                    startPointValues = GetValuesAroundSCPoint( watershed2d, scPointsListAdj[i][0] ),
                                    endPointValues   = GetValuesAroundSCPoint( watershed2d, scPointsListAdj[i][-1] ),
                                    identifier=identifier )
                matchingSCs = [ sc for sc in scList if sc.values==newSC.values ] # match any subcoutours in cVLS so far that are for the same pair of cells
                matchingSCs = [ sc for sc in matchingSCs if totuple(sc.points[::-1])==totuple(newSC.points) ] # Only keep subcoutours where the points match the reverse of the points in newSC
                                #sorted([newSC.points[0],newSC.points[-1]]) == sorted([sc.points[0],sc.points[-1]]) ] # Should only possibly find 1 match...
                if matchingSCs==[]: # This is a new subContour, not a duplicate!
                    scList.append(newSC)
                    contourOrderingByValue[v].append( (identifier,True) )
                    identifier+=1
                else:
                    matchingSCs[0].adjusted_length = min( matchingSCs[0].adjusted_length,
                                                          newSC.adjusted_length ) # keep the minimum perimeter length...
                    contourOrderingByValue[v].append( (matchingSCs[0].identifier,False) ) # False means the subcountour is backwards for this cell!
    scList.sort(cmpGen(lambda x: x.values)) # was just cVLS.sort()... this works, I hope?
    IDs = [sc.identifier for sc in scList]
    for sc in scList:      # scrub the id's, probably not necessary... 
        sc.identifier=None
    
    # Reindex after sorting...
    for v in allValues:
        contourOrderingByValue[v] = [ (IDs.index(i),d) for i,d in contourOrderingByValue[v] ]
    
    return CellNetwork( subContours=scList , contourOrderingByValue=contourOrderingByValue , allValues=allValues )

def GetCellNetworksByFrame(watershed,allValsByFrame):
    '''Get a list of CellNetworks based on a watershed segmentation'''
    return [ GetCellNetwork(watershed[i],allValsByFrame[i])
            for i in range(len(watershed)) ]

def GetXYListAndPolyListFromCellNetworkList(cellNetworkList,closeLoops=True):
    '''Get a multi-frame xyList and polyList'''
    ret = [ cn.GetXYListAndPolyList(closeLoops=closeLoops) for cn in cellNetworkList ]
    xyList,polyList = zip(*ret) # effectively like transpose...
    return xyList,polyList

def GetCVDListFromCellNetworkList(cellNetworkList):
    return [ cn.GetCellularVoronoiDiagram() for cn in cellNetworkList ]

def GetCellNetworkListWithLimitedPointsBetweenNodes(cellNetworkList,splitLength=1,fixedNumInteriorPoints=None,interpolate=True,checkForDegeneracy=True):
    '''Based on matching subcontours by value pair, this function defines a fixed number of interior points for each subcontour
       and then applies this "trimming" procedure equitably to each frame in cellNetworkList (uses LimitPointsBetweenNodes)'''
    #allValues = sorted(set( [ v for cn in cellNetworkList for v in cn.allValues ] )) # not used...
    allPairs = sorted(set( [ tuple(sc.values) for cn in cellNetworkList for sc in cn.subContours ] )) # Value pairs...
    
    # Build the numInteriorPointsDict:
    if fixedNumInteriorPoints:
        numInteriorPointsDict = { p:fixedNumInteriorPoints for p in allPairs }
    else:
        # minLength is the number of points of the shortest subcountour between cells p[0] and p[1] from all frames
        minLength = { p : min( [ sc.numPoints
                                for cn in cellNetworkList
                                for sc in cn.subContours
                                if tuple(sc.values)==p ] )
                     for p in allPairs }
        numInteriorPointsDict = { p:(minLength[p]//splitLength) for p in allPairs }
    
    cellNetworkListNew = deepcopy(cellNetworkList) # otherwise, we'd also change the input argument in the outside world!
    for cn in cellNetworkListNew:
        cn.LimitPointsBetweenNodes(numInteriorPointsDict,interpolate=interpolate,checkForDegeneracy=False)
    
    if checkForDegeneracy:
        problemValsList,problemValuePairsList = ziptranspose([ cn.CheckForDegenerateContours() for cn in cellNetworkListNew ])
        if flatten(problemValsList)!=[]:
            print 'All degenerate values:',sorted(set(flatten(problemValsList)))
            print 'Degenerate values by frame:'
            for i,problemVals in enumerate(problemValsList):
                if problemVals!=[]:
                    print ' ',i,':',problemVals
            print 'Degenerate contours (zero area) found between these cellIDs on these frames!'
            for i,problemValuePairs in enumerate(problemValuePairsList):
                if problemValuePairs!=[]:
                    print ' ',i,':',problemValuePairs
            raise DegenerateNetworkException('Degeneracy check failed!')
    
    return cellNetworkListNew

def collectPreviousNextResults(prv,nxt):
    '''A generic function that joins elements from prv and nxt like so: [prv[0],...,prv[i]+nxt[i-1],...,nxt[-1]]
       prv and nxt are lists of lists and were mapped (with some function)
       from elements of an original list (length n) with this relation:
           prv[0 -> n-1]  maps to  original[0 -> n-1]
           nxt[0 -> n-1]  maps to  original[1 -> n]
       again, prv and nxt must be lists of lists'''
    return [prv[0]] + [prv[i]+nxt[i-1] for i in range(1,len(prv))] + [nxt[-1]]

def GetMatchedCellNetworkListsWithLimitedPointsBetweenNodes(cellNetworkListPrev,cellNetworkListNext,splitLength=1,fixedNumInteriorPoints=None,interpolate=True,checkForDegeneracy=True):
    '''Uses GetCellNetworkListWithLimitedPointsBetweenNodes to match each pair between cellNetworkListPrev and cellNetworkListNext'''
    cnListPrevNew = []
    cnListNextNew = []
    for i in range(len(cellNetworkListPrev)):
        cnl = GetCellNetworkListWithLimitedPointsBetweenNodes([cellNetworkListPrev[i],cellNetworkListNext[i]],splitLength,fixedNumInteriorPoints,interpolate,checkForDegeneracy=False)
        cnListPrevNew.append(cnl[0])
        cnListNextNew.append(cnl[1])
    
    if checkForDegeneracy:
        problemValsListPrev,problemValuePairsListPrev = ziptranspose([ cn.CheckForDegenerateContours() for cn in cnListPrevNew ])
        problemValsListNext,problemValuePairsListNext = ziptranspose([ cn.CheckForDegenerateContours() for cn in cnListNextNew ])
        problemValsList = [ sorted(set(i))
                           for i in collectPreviousNextResults(problemValsListPrev,problemValsListNext) ]
        problemValuePairsList = [ sorted(set(i))
                                 for i in collectPreviousNextResults(problemValuePairsListPrev,problemValuePairsListNext) ]
        
        if flatten(problemValsList)!=[]:
            print 'All degenerate values:',sorted(set(flatten(problemValsList)))
            print 'Degenerate values by frame:'
            for i,problemVals in enumerate(problemValsList):
                if problemVals!=[]:
                    print ' ',i,':',problemVals
            print 'Degenerate contours (zero area) found between these cellIDs on these frames!'
            for i,problemValuePairs in enumerate(problemValuePairsList):
                if problemValuePairs!=[]:
                    print ' ',i,':',problemValuePairs
            raise DegenerateNetworkException('Degeneracy check failed!')
    
    return cnListPrevNew,cnListNextNew

def GetXYListAndPolyListWithLimitedPointsBetweenNodes(cellNetworkList,splitLength=1,fixedNumInteriorPoints=None,interpolate=True):
    '''Get a list of points and a set of polygons network from a cellNetwork limit points between triple junctions
       (Applies GetCellNetworkListWithLimitedPointsBetweenNodes and then GetXYListAndPolyListFromCellNetworkList)'''
    return GetXYListAndPolyListFromCellNetworkList(
             GetCellNetworkListWithLimitedPointsBetweenNodes(cellNetworkList,splitLength,fixedNumInteriorPoints,interpolate) )

def GetMatchedCellNetworksCollapsing(cnA,cnB):
    '''Make 2 simplified cell networks, making sure that there is a 1-to-1 mapping between all subcontours
       This function removes values that are not common to both networks and collapses
       pairs of subcontours that do not match but are in between the same 4 cells'''
    
    if cnA==cnB: # if we got the same object for some reason, just return 2 shallow clones
        return cnA,cnB
    
    # sharedVals = sorted(set(cnA.allValues+cnB.allValues)) # not used...
    valsNotInA = sorted(set(cnB.allValues).difference(cnA.allValues))
    valsNotInB = sorted(set(cnA.allValues).difference(cnB.allValues))
    
    # Delete any values that are not in both, replacing with background...
    cnA,cnB = deepcopy(cnA),deepcopy(cnB) # Make copies so we don't modify the originals
    cnA.RemoveValues(valsNotInB)
    cnB.RemoveValues(valsNotInA)
    
    matchedInB,removeFromA_a,removeFromB_a,notRecoverableInA = cnA.FindMatchesAndRemovals(cnB)
    unMatchedInB = [ i for i in range(len(cnB.subContours)) if i not in matchedInB ] # This lets us skip the indexes that already matched
    
    _,removeFromB_b,removeFromA_b,notRecoverableInB = cnB.FindMatchesAndRemovals(cnA,searchInds = unMatchedInB) # FLIP
    
    print 'Summary of pairs that are not recoverable from A:',notRecoverableInA
    print 'Summary of pairs that are not recoverable from B:',notRecoverableInB
    
    cnA.RemoveMultipleSubContours(removeFromA_a + removeFromA_b)
    cnB.RemoveMultipleSubContours(removeFromB_a + removeFromB_b)
    
    return cnA,cnB

def GetMatchedCellNetworksCollapsingWithLimitedPoints(cnA,cnB,splitLength=1,fixedNumInteriorPoints=None,interpolate=True):
    '''Make 2 simplified cell networks, making sure that there is a 1-to-1 mapping between all points; this function collapses
       pairs of subcontours that do not match but are in between the same 4 cells'''
    
    cnANew,cnBNew = GetMatchedCellNetworksCollapsing(cnA,cnB)
    cnALim,cnBLim = GetCellNetworkListWithLimitedPointsBetweenNodes( [cnANew,cnBNew],splitLength,
                                                                     fixedNumInteriorPoints,interpolate)
    
    return cnALim,cnBLim

def GetCVDListStatic( waterArr,d,vfmParameters,
                      extraRemoveValsByFrame=[],splitLength=20, fixedNumInteriorPoints=None,
                      usePlot=False,forceRemake=False ):
    '''Get cvdLists from a waterArr, ignoring any differences between frames. This function:
        * removes values in extraRemoveValsByFrame,
        * limits points between nodes
        * compresses cellID's to begin at 1 (excludes background)
        * returns a cvdList
       
       Throws error if static flag is not set; in this case,
        use GetMatchedCVDListPrevNext instead'''
    
    allValsByFrame = [ np.unique(i)[1:] for i in waterArr ] # Skip background
    
    # Ensure we're doing a static analysis
    assert vfmParameters.useStaticAnalysis, 'useStaticAnalysis is not set! Did you mean to use the function GetMatchedCVDListPrevNext?'
    # Ensure we got a list-of-lists for extraRemovalsByFrame
    assert all([ hasattr(vals,'__iter__') for vals in extraRemoveValsByFrame ])
    # Ensure that this has enough elements, if not, add more empty lists
    extraRemoveValsByFrame += [[] for i in range(len(waterArr)-len(extraRemoveValsByFrame))]
    
    cnListStaticFile = os.path.join(d,'cellNetworksListStatic.pickle') # Saved cellNetworks file
    
    loadCnListFromFile = os.path.exists(cnListStaticFile) and not any(extraRemoveValsByFrame) and not forceRemake
    
    if loadCnListFromFile:
        print 'Reloading cnLists from file:',cnListStaticFile
        cnList = cPickle.load(open(cnListStaticFile,'r'))
        
        if len(cnList)!=len(waterArr):
            print 'cnList is the wrong length in the file:',cnListStaticFile
            print 'Will remake them'
            loadCnListFromFile=False
    if not loadCnListFromFile:
        cnList = []
        for i in range(len(waterArr)):
            water = np.array(waterArr[i])
            valsToKeep = sorted( set(allValsByFrame[i]).difference(extraRemoveValsByFrame[i]) )
            for v in sorted(extraRemoveValsByFrame[i]):
                water[np.where(water==v)] = 1
            
            # next, create the CellNetwork and append it to the list:
            cnList.append( GetCellNetwork(water,valsToKeep) )
        
        # Only save this if we're using all the values; otherwise it gets confusing!
        if not any(extraRemoveValsByFrame):
            print 'Saving cnList to file:',cnListStaticFile
            cPickle.dump(cnList,open(cnListStaticFile,'w'))
    
    cnListLim = GetCellNetworkListWithLimitedPointsBetweenNodes(cnList,splitLength,fixedNumInteriorPoints,interpolate=True)
    cvdList = GetCVDListFromCellNetworkList(cnListLim)
    
    if usePlot:
        print 'Plotting the subContours:'
        cnList[0].scPlotT('r-')
        cnListLim[0].scPlotT('ro-')
    
    return cvdList

def GetMatchedCVDListPrevNext( waterArr,d,vfmParameters,
                               extraRemoveValsByFrame=[],splitLength=20, fixedNumInteriorPoints=None,
                               usePlot=False,forceRemake=False ):
    '''Get matched before and after cvdLists from a waterArr. This function:
        * removes values in extraRemoveValsByFrame,
        * matches subContours between CellNetworks
            - since this is the slowest step, these 2 matched networks are saved to a file in d and reloaded if this is run again
              (with the caveat that save/load will not occur if extraRemoveValsByFrame are specified)
              load can be overridden with forceRemake
        * limits points between nodes
        * compresses cellID's to begin at 1 (excludes background)
        * returns 2 cvdLists
       
       Throws error if static flag is set; in this case,
        use GetCVDListStatic instead'''
    
    allValsByFrame = [ np.unique(i)[1:] for i in waterArr ] # Skip background
    
    # Ensure we're doing a dynamic analysis (if you just want to get rid of viscous effects, set viscosityTimeStepRatio to 0)
    assert not vfmParameters.useStaticAnalysis, 'useStaticAnalysis is set! Did you mean to use the function GetCVDListStatic?'
    # Ensure we got a list-of-lists for extraRemovalsByFrame
    assert all([ hasattr(vals,'__iter__') for vals in extraRemoveValsByFrame ])
    # Ensure that this has enough elements, if not, add more empty lists
    extraRemoveValsByFrame += [[] for i in range(len(waterArr)-1-len(extraRemoveValsByFrame))]
    
    cnListPrevAndNextFile = os.path.join(d,'cellNetworksListPrevAndNext.pickle') # Saved cellNetworks file
    
    loadCnListFromFile = os.path.exists(cnListPrevAndNextFile) and not any(extraRemoveValsByFrame) and not forceRemake
    
    if loadCnListFromFile:
        print 'Reloading cnLists from file:',cnListPrevAndNextFile
        cnListPrev,cnListNext = cPickle.load(open(cnListPrevAndNextFile,'r'))
        
        if len(cnListPrev)!=len(waterArr)-1:
            print 'cnLists are the wrong length in the file:',cnListPrevAndNextFile
            print 'Will remake them'
            loadCnListFromFile=False
    if not loadCnListFromFile:
        cnListPrev = []
        cnListNext = []
        for i in range(len(waterArr)-1):
            # create matched arrays first:
            waterA,waterB = np.array(waterArr[i]) , np.array(waterArr[i+1])
            extraRemovalsAB = set(extraRemoveValsByFrame[i]).union(extraRemoveValsByFrame[i+1])
            commonVals = sorted( set(allValsByFrame[i]).intersection(allValsByFrame[i+1]).difference(extraRemovalsAB) )
            valsRemoveA = set(allValsByFrame[i]).difference(commonVals)
            valsRemoveB = set(allValsByFrame[i+1]).difference(commonVals)
            for v in sorted(valsRemoveA):
                waterA[np.where(waterA==v)] = 1
            for v in sorted(valsRemoveB):
                waterB[np.where(waterB==v)] = 1
            
            # next, create matched CellNetworks:
            cnA = GetCellNetwork(waterA,commonVals)
            cnB = GetCellNetwork(waterB,commonVals)
            
            cnA,cnB = GetMatchedCellNetworksCollapsing(cnA,cnB)
            cnListPrev.append(cnA)
            cnListNext.append(cnB)
        # Only save this if we're using all the values; otherwise it gets confusing!
        if not any(extraRemoveValsByFrame):
            print 'Saving cnLists to file:',cnListPrevAndNextFile
            cPickle.dump([cnListPrev,cnListNext],open(cnListPrevAndNextFile,'w'))
    
    cnListPrevLim,cnListNextLim = GetMatchedCellNetworkListsWithLimitedPointsBetweenNodes(cnListPrev,cnListNext,splitLength,fixedNumInteriorPoints,interpolate=True)
    
    cvdListPrev = GetCVDListFromCellNetworkList(cnListPrevLim)
    cvdListNext = GetCVDListFromCellNetworkList(cnListNextLim)
    
    if usePlot:
        print 'Plotting the subContours:'
        cnListPrev[0].scPlotT('r-')
        cnListNext[0].scPlotT('b-')
        cnListPrevLim[0].scPlotT('ro-')
        cnListNextLim[0].scPlotT('bo-')
    
    return cvdListPrev,cvdListNext

def SaveXYListAndPolyListToMMAFormat(xyList,polyList,filename,bumpIndsUp1=True,removeLastPoint=True):
    '''xyList: nested list of xy pairs for each time point.
       polyList: nested list of dictionaries for each time point where
                 each entry is like: {cellID: [ <indexes into xyList> ]}
       Exports a MMA compatible dataStructure also called "polyList" which looks like lists of:
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
    '''Save a cellNetwork to the MMA format
       (basically just GetXYListAndPolyListFromCellNetworkList followed by SaveXYListAndPolyListToMMAFormat)'''
    xyList,polyList = GetXYListAndPolyListFromCellNetworkList(cellNetworkList)
    SaveXYListAndPolyListToMMAFormat(xyList,polyList,filename,bumpIndsUp1,removeLastPoint)


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
