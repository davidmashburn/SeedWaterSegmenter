#from numpy import ones, zeros, nonzero, array, put, take, argmin, \
#                  transpose, compress, concatenate, where, uint8, \
#                  asarray, minimum, maximum, float64, bool, uint8, \
#                  uint16, int32, ravel, nonzero, newaxis
import numpy as np
from string import upper
from heapq import heapify, heappush, heappop
import Cpyx
from Pyx.SeedWaterImports import *

# Functions I don't know:
# secross,limits,binary,intersec,pad4n,label,isbinary,se2flatidx,to_int32,sedilate,sesum,mat2set,set2mat

def add4dilate(f, c):
    if not c:
        return f
    y = np.asarray(f,np.float64) + c
    k1,k2 = limits(f)
    y = ((f==k1) * k1) + ((f!=k1) * y)
    y = np.maximum(np.minimum(y,k2),k1)
    return y.astype(f.dtype)

def limits(f):
    code = f.dtype
    if code == np.bool: return np.array([0,1])
    if code == np.uint8: return np.array([0,255])
    if code == np.uint16: return np.array([0,65535])
    if code == np.int32: return np.array([-2147483647,2147483647])

    raise ValueError('pymorph.limits: does not accept this typecode: %s' % code)

def setrans(Bi, t):
    x,v=mat2set(Bi)
    Bo = set2mat((x+t,v))
    return Bo.astype(Bi.dtype)

def mat2set(A):
    if len(A.shape) == 1: A = A[np.newaxis,:]
    offsets = np.nonzero(np.ravel(A) - limits(A)[0])[0]
    if len(offsets) == 0: return ([],[])
    h,w = A.shape
    x = [0,1]
    x[0] = offsets//w - (h-1)//2
    x[1] = offsets%w - (w-1)//2
    x = np.transpose(x)
    return x,np.take(np.ravel(A),offsets)


def set2mat(A):
    if len(A) == 2:
        x, v = A
        v = np.asarray(v)
    elif len(A) == 1:
        x = A[0]
        v = np.ones((len(x),), np.uint8)
    else:
        raise TypeError, 'pymorph.set2mat: argument must be a tuple of length 1 or 2'
    if len(x) == 0:  return np.array([0]).astype(v.dtype)
    if len(x.shape) == 1: x = x[np.newaxis,:]
    dh,dw = abs(x).max(0)
    h,w = (2*dh)+1, (2*dw)+1
    M=np.ones((h,w),np.int32) * limits(v)[0]
    offset = x[:,0] * w + x[:,1] + (dh*w + dw)
    np.put(M,offset,v)
    return M.astype(v.dtype)


def isbinary(f):
    return f.dtype == bool

def binary(f, k=1):
    f = np.asanyarray(f)
    return (f >= k)

def se2flatidx(f,Bc):
    h,w=Bc.shape
    h //= 2
    w //= 2
    Bi=[]
    for i,j in zip(*np.where(Bc)):
        Bi.append( (j-w)+(i-h)*f.shape[1] )
    return np.array(Bi)

def seunion(B1, B2):
    assert B1.dtype == B2.dtype, 'Cannot have different datatypes:'
    type1 = B1.dtype
    if len(B1) == 0: return B2
    if len(B1.shape) == 1: B1 = B1[np.newaxis,:]
    if len(B2.shape) == 1: B2 = B2[np.newaxis,:]
    if B1.shape != B2.shape:
        inf = limits(B1)[0]
        h1,w1 = B1.shape
        h2,w2 = B2.shape
        H,W = max(h1,h2),max(w1,w2)
        Hc,Wc = (H-1)//2,(W-1)//2    # center
        BB1,BB2 = np.asarray(B1),np.asarray(B2)
        B1, B2  = inf * np.ones((H,W),np.int32), inf *np.ones((H,W),np.int32)
        dh1s , dh1e = (h1-1)//2 , (h1-1)//2 + (h1+1)%2 # deal with even and odd dimensions
        dw1s , dw1e = (w1-1)//2 , (w1-1)//2 + (w1+1)%2
        dh2s , dh2e = (h2-1)//2 , (h2-1)//2 + (h2+1)%2
        dw2s , dw2e = (w2-1)//2 , (w2-1)//2 + (w2+1)%2
        B1[ Hc-dh1s : Hc+dh1e+1  ,  Wc-dw1s : Wc+dw1e+1 ] = BB1
        B2[ Hc-dh2s : Hc+dh2e+1  ,  Wc-dw2s : Wc+dw2e+1 ] = BB2
    B = np.maximum(B1,B2).astype(type1)
    return B

def sedilate(B1, B2):
    assert (isbinary(B1) or (B1.dtype == np.int32)) and (isbinary(B2) or B2.dtype == np.int32), 'pymorph.sedilate: s.e. must be binary or int32'
    if len(B1.shape) == 1: B1 = B1[np.newaxis,:]
    if len(B2.shape) == 1: B2 = B2[np.newaxis,:]
    if B1.dtype==np.int32 or B2.dtype == np.int32:
       Bo = to_int32([limits(to_int32([0]))[0]])
       if isbinary(B1):
          B1 = gray(B1,'int32',0)
       if isbinary(B2):
          B2 = gray(B2,'int32',0)
    else:
       Bo = binary([0])
    x,v = mat2set(B2)
    if len(x):
        for i in xrange(x.shape[0]):
            s = add4dilate(B1,v[i])
            st= setrans(s,x[i])
            Bo = seunion(Bo,st)
    return Bo

def sesum(B=None, N=1):
    if B is None: B = secross()
    if N==0:
        if isbinary(B): return binary([[1]])
        else:           return to_int32([[0]]) # identity
    NB = B
    for i in xrange(N-1):
        NB = sedilate(NB,B)
    return NB


def secross(r=1):
    return sesum(binary([[0,1,0],
                         [1,1,1],
                         [0,1,0]]),
                 r)

def seshow(b, option="normal"):
    option = upper(option)
    if option=='NON-FLAT':
        y = to_int32([0])
        if isbinary(b):
            b = intersec(gray(b,'int32'),0)
    elif option=='NORMAL':
        if isbinary(b):
            y = binary([1])
        else:
           y = to_int32([0])
    elif option=='EXPAND':
        assert isbinary(b), 'pymorph.seshow: \'expand\' option is only available with flat SE'
        y = sedilate(binary([1]),b)
        b1 = (y>=0)
        b0 = erode(y,b)
        return bshow(b1,y,b0)
    else:
        assert False, 'pymorph.seshow: not a valid flag: NORMAL, EXPAND or NON-FLAT'
    return sedilate(y,b)

def pad4n(f, Bc, value, scale=1):
    if type(Bc) is not np.array:
      Bc = seshow(Bc)
    Bh, Bw = Bc.shape
    assert Bh%2 and Bw%2, 'structuring element must be odd sized'
    ch, cw = scale * Bh/2, scale * Bw/2
    g = value * np.ones( f.shape + scale * (np.array(Bc.shape) - 1))
    g[ ch: -ch, cw: -cw] = f
    return g.astype(f.dtype)

def cwloop(f,i,costM,yflat,hqueue,status,Bi,is_gvoronoi,return_lines,y1flat=None):
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

#exec(Cpyx.CythonInline('''
#from __future__ import division
#import numpy as np
#cimport numpy as np
#cimport cython
#from string import upper
#from heapq import heapify, heappush, heappop

#def cwloopC32(np.ndarray[np.uint16_t, ndim=1] f not None,
            #int i,
            #np.ndarray[np.float64_t, ndim=1] costM not None,
            #np.ndarray[np.int32_t, ndim=1] yflat not None,
            #list hqueue,
            #np.ndarray[np.uint8_t, ndim=1] status not None,
            #np.ndarray[np.int32_t, ndim=1] Bi not None,
            #bool is_gvoronoi,
            #bool return_lines,
            #y1flat=None):
    
    #cdef float ncost
    #cdef float cost
    #cdef int _
    #cdef int pi
    #cdef int qi
    
    #while hqueue:
        #cost,_,pi = heappop(hqueue)
        #status[pi]=1                                            # make it a permanent label
        #for qi in Bi+pi:                                        # for each neighbor of pi
            #if (status[qi] != 3):                               # not image border
                #if (status[qi] != 1):                           # if not permanent
                    #if is_gvoronoi:
                        #ncost=costM[pi]
                    #else:
                        #ncost=f[qi]
                    #if ncost < costM[qi]:
                        #costM[qi]=ncost
                        #yflat[qi] = yflat[pi]                   # propagate the label
                        #heappush(hqueue,(ncost,i,qi))
                        #i += 1
                #elif (return_lines  and
                     #(yflat[qi] != yflat[pi]) and
                     #(y1flat[qi] == 0)     ):
                    #y1flat[pi] = 1
    ## Modifies  hqueue  and arrays:  yflat, y1flat, status, costM

#def cwloopC64(np.ndarray[np.uint16_t, ndim=1] f not None,
            #int i,
            #np.ndarray[np.float64_t, ndim=1] costM not None,
            #np.ndarray[np.int64_t, ndim=1] yflat not None,
            #list hqueue,
            #np.ndarray[np.uint8_t, ndim=1] status not None,
            #np.ndarray[np.int64_t, ndim=1] Bi not None,
            #bool is_gvoronoi,
            #bool return_lines,
            #y1flat=None):
    
    #cdef float ncost
    #cdef float cost
    #cdef int _
    #cdef int pi
    #cdef int qi
    
    #while hqueue:
        #cost,_,pi = heappop(hqueue)
        #status[pi]=1                                            # make it a permanent label
        #for qi in Bi+pi:                                        # for each neighbor of pi
            #if (status[qi] != 3):                               # not image border
                #if (status[qi] != 1):                           # if not permanent
                    #if is_gvoronoi:
                        #ncost=costM[pi]
                    #else:
                        #ncost=f[qi]
                    #if ncost < costM[qi]:
                        #costM[qi]=ncost
                        #yflat[qi] = yflat[pi]                   # propagate the label
                        #heappush(hqueue,(ncost,i,qi))
                        #i += 1
                #elif (return_lines  and
                     #(yflat[qi] != yflat[pi]) and
                     #(y1flat[qi] == 0)     ):
                    #y1flat[pi] = 1
    ## Modifies  hqueue  and arrays:  yflat, y1flat, status, costM
#'''))

def cwatershed(f, markers, Bc=None, return_lines=False,is_gvoronoi=False):
    if Bc is None: Bc = secross()
    if isbinary(markers):
        markers = label(markers,Bc)
    status = pad4n(np.zeros(f.shape,np.uint8),Bc, 3)
    f = pad4n(f,Bc,0)                       # pad input image with 0
    y = pad4n(markers,Bc,0)                 # pad marker image with 0
    if return_lines:
        y1 = intersec(binary(y), 0)
        y1flat=y1.ravel()
    else:
        y1flat=None
    costM = limits(f)[1] * np.ones(f.shape)  # cummulative cost function image
    # get 1D displacement neighborhood pixels
    Bi=se2flatidx(f,Bc)
    status=status.ravel()
    f=f.ravel()
    costM=costM.ravel()
    yflat=y.ravel()

    initial, = np.where(yflat > 0)
    costM[initial]=f[initial]
    # I sort in order of insertion for the simple reason that that's what
    # the original code intended to do (although that code was not functional)
    hqueue=[(costM[idx],i,idx) for i,idx in enumerate(initial)]
    heapify(hqueue)
    
    if np.array(1).dtype==np.int32:
        cwloopC32(f,i,costM,yflat,hqueue,status,Bi,is_gvoronoi,return_lines,y1flat=y1flat)
    elif np.array(1).dtype==np.int64:
        cwloopC64(f,i,costM,yflat,hqueue,status,Bi,is_gvoronoi,return_lines,y1flat=y1flat)
    
    #while hqueue:
    #    cost,_,pi = heappop(hqueue)
    #    status[pi]=1                                            # make it a permanent label
    #    for qi in Bi+pi:                                        # for each neighbor of pi
    #        if (status[qi] != 3):                               # not image border
    #            if (status[qi] != 1):                           # if not permanent
    #                if is_gvoronoi:
    #                    ncost=costM[pi]
    #                else:
    #                    ncost=f[qi]
    #                if ncost < costM[qi]:
    #                    costM[qi]=ncost
    #                    yflat[qi] = yflat[pi]                   # propagate the label
    #                    heappush(hqueue,(ncost,i,qi))
    #                    i += 1
    #            elif (return_lines  and
    #                 (yflat[qi] != yflat[pi]) and
    #                 (y1flat[qi] == 0)     ):
    #                y1flat[pi] = 1
    y=y[1:-1,1:-1]
    if return_lines:
        y1=y1[1:-1,1:-1]
        return y,y1
    return y
