from __future__ import division
import numpy as np
cimport numpy as np
cimport cython
from string import upper
from heapq import heapify, heappush, heappop

def cwloopC32(np.ndarray[np.uint16_t, ndim=1] f not None,
            int i,
            np.ndarray[np.float64_t, ndim=1] costM not None,
            np.ndarray[np.int32_t, ndim=1] yflat not None,
            list hqueue,
            np.ndarray[np.uint8_t, ndim=1] status not None,
            np.ndarray[np.int32_t, ndim=1] Bi not None,
            int is_gvoronoi,
            int return_lines,
            y1flat=None):
    
    cdef float ncost
    cdef float cost
    cdef int _
    cdef int pi
    cdef int qi
    
    while hqueue:
        cost,_,pi = heappop(hqueue)
        status[pi]=1                                            # make it a permanent label
        for qi in Bi+pi:                                        # for each neighbor of pi
            if (status[qi] != 3):                               # not image border
                if (status[qi] != 1):                           # if not permanent
                    if is_gvoronoi:
                        ncost=costM[pi]
                    else:
                        ncost=f[qi]
                    if ncost < costM[qi]:
                        costM[qi]=ncost
                        yflat[qi] = yflat[pi]                   # propagate the label
                        heappush(hqueue,(ncost,i,qi))
                        i += 1
                elif (return_lines  and
                     (yflat[qi] != yflat[pi]) and
                     (y1flat[qi] == 0)     ):
                    y1flat[pi] = 1
    # Modifies  hqueue  and arrays:  yflat, y1flat, status, costM

def cwloopC64(np.ndarray[np.uint16_t, ndim=1] f not None,
            int i,
            np.ndarray[np.float64_t, ndim=1] costM not None,
            np.ndarray[np.int64_t, ndim=1] yflat not None,
            list hqueue,
            np.ndarray[np.uint8_t, ndim=1] status not None,
            np.ndarray[np.int64_t, ndim=1] Bi not None,
            int is_gvoronoi,
            int return_lines,
            y1flat=None):
    
    cdef float ncost
    cdef float cost
    cdef int _
    cdef int pi
    cdef int qi
    
    while hqueue:
        cost,_,pi = heappop(hqueue)
        status[pi]=1                                            # make it a permanent label
        for qi in Bi+pi:                                        # for each neighbor of pi
            if (status[qi] != 3):                               # not image border
                if (status[qi] != 1):                           # if not permanent
                    if is_gvoronoi:
                        ncost=costM[pi]
                    else:
                        ncost=f[qi]
                    if ncost < costM[qi]:
                        costM[qi]=ncost
                        yflat[qi] = yflat[pi]                   # propagate the label
                        heappush(hqueue,(ncost,i,qi))
                        i += 1
                elif (return_lines  and
                     (yflat[qi] != yflat[pi]) and
                     (y1flat[qi] == 0)     ):
                    y1flat[pi] = 1
    # Modifies  hqueue  and arrays:  yflat, y1flat, status, costM

def PointsToArray32( list seedList,
                     list seedVals,
                     np.ndarray[np.int32_t, ndim=2] seedArray not None):
    '''Update seedArray from seedList'''
    cdef int i
    cdef list s
    seedArray[:,:] = 0
    for i,s in enumerate(seedList):
        seedArray[s[0],s[1]]=seedVals[i]

def PointsToArray64( list seedList,
                     list seedVals,
                     np.ndarray[np.int64_t, ndim=2] seedArray not None):
    '''Update seedArray from seedList'''
    cdef int i
    cdef list s
    seedArray[:,:] = 0
    for i,s in enumerate(seedList):
        seedArray[s[0],s[1]]=seedVals[i]

def convToRandColors32(np.ndarray[np.uint8_t, ndim=2] mapPlotRandomArray not None,
                     np.ndarray[np.uint8_t, ndim=3] rgbM not None,
                     np.ndarray[np.int32_t, ndim=2] water not None,
                     int x, int y):
    cdef int i,j,k
    for i in range(x):
        for j in range(y):
            for k in range(3):
                rgbM[i,j,k]=mapPlotRandomArray[k,water[i,j]]

def convToRandColors64(np.ndarray[np.uint8_t, ndim=2] mapPlotRandomArray not None,
                     np.ndarray[np.uint8_t, ndim=3] rgbM not None,
                     np.ndarray[np.int64_t, ndim=2] water not None,
                     int x, int y):
    cdef int i,j,k
    for i in range(x):
        for j in range(y):
            for k in range(3):
                rgbM[i,j,k]=mapPlotRandomArray[k,water[i,j]]
