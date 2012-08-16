import numpy as np
import Cpyx

CREATE_CYTHON_VERSION = False

# Original Documentation:
'''/*
Best-fitting ellipse routines by:

  Bob Rodieck
  Department of Ophthalmology, RJ-10
  University of Washington, 
  Seattle, WA, 98195

Notes on best-fitting ellipse:

  Consider some arbitrarily shaped closed profile, which we wish to
  characterize in a quantitative manner by a series of terms, each 
  term providing a better approximation to the shape of the profile.  
  Assume also that we wish to include the orientation of the profile 
  (i.e. which way is up) in our characterization. 

  One approach is to view the profile as formed by a series harmonic 
  components, much in the same way that one can decompose a waveform
  over a fixed interval into a series of Fourier harmonics over that 
  interval. From this perspective the first term is the mean radius,
  or some related value (i.e. the area).  The second term is the 
  magnitude and phase of the first harmonic, which is equivalent to the
  best-fitting ellipse.  

  What constitutes the best-fitting ellipse?  First, it should have the
  same area.  In statistics, the measure that attempts to characterize some
  two-dimensional distribution of data points is the 'ellipse of 
  concentration' (see Cramer, Mathematical Methods of Statistics, 
  Princeton Univ. Press, 945, page 283).  This measure equates the second
  order central moments of the ellipse to those of the distribution, 
  and thereby effectively defines both the shape and size of the ellipse. 

  This technique can be applied to a profile by assuming that it constitutes
  a uniform distribution of points bounded by the perimeter of the profile.
  For most 'blob-like' shapes the area of the ellipse is close to that
  of the profile, differing by no more than about 4%. We can then make
  a small adjustment to the size of the ellipse, so as to give it the 
  same area as that of the profile.  This is what is done here, and 
  therefore this is what we mean by 'best-fitting'. 

  For a real pathologic case, consider a dumbell shape formed by two small
  circles separated by a thin line. Changing the distance between the
  circles alters the second order moments, and thus the size of the ellipse 
  of concentration, without altering the area of the profile. 

public class Ellipse_Fitter implements PlugInFilter {
	public int setup(String arg, ImagePlus imp) {
		return DOES_ALL;
	}
	public void run(ImageProcessor ip) {
		EllipseFitter ef = new EllipseFitter();
		ef.fit(ip);
		IJ.write(IJ.d2s(ef.major)+" "+IJ.d2s(ef.minor)+" "+IJ.d2s(ef.angle)+" "+IJ.d2s(ef.xCenter)+" "+IJ.d2s(ef.yCenter));
		ef.drawEllipse(ip);
	}
}
*/'''

# I decided just to gank the algorithm straight from ImageJ...
# I changed it to assume that you know the bounding rectangle,
# and have sliced the image array accordingly (this is arr)
def EllipseFitter(arr,usePrint=False):
    left=0 # holdover from the ImageJ version
    top=0 # holdover from the ImageJ version
    width=arr.shape[0]
    height=arr.shape[1]
    
    HALFPI = np.pi/2

    #double xCenter # X centroid
    #double yCenter # Y centroid
    #double major # Length of major axis
    #double minor # Length of minor axis
    #double angle # Angle in degrees
    #double theta # Angle in radians
    #int[] xCoordinates # Initialized by makeRoi()
    #int[] yCoordinates # Initialized by makeRoi()
    nCoordinates = 0 #int # Initialized by makeRoi()
    bitCount = 0 #int
    #double  xsum, ysum, x2sum, y2sum, xysum
    #byte[] mask
    #int left, top, width, height
    #double   n
    #double   xm, ym   #mean values
    #double   u20, u02, u11  #central moments
    #ImageProcessor ip
    ##private double pw, ph
    #boolean record

    sqrtPi = np.sqrt(np.pi)
    #double    a11, a12, a22, m4, z, scale, tmp, xoffset, yoffset
    #double    RealAngle

    #if mask==None: # mask should never be None
    #    major = (width*2) / sqrtPi
    #    minor = (height*2) / sqrtPi # * Info->PixelAspectRatio
    #    angle = 0.0
    #    theta = 0.0
    #    if major < minor:
    #        tmp = major
    #        major = minor
    #        minor = tmp
    #        angle = 90.0
    #        theta = np.pi/2.0
    #    xCenter = left + width / 2.0
    #    yCenter = top + height / 2.0
    #    return
    
    # computeSums
    xsum = 0.0
    ysum = 0.0
    x2sum = 0.0
    y2sum = 0.0
    xysum = 0.0
    #int bitcountOfLine
    #double   xe, ye
    #int xSumOfLine
    for y in range(height):
        bitcountOfLine = 0
        xSumOfLine = 0
        #offset = y*width # int
        for x in range(width):
            if arr[x,y] != 0:
                bitcountOfLine+=1
                xSumOfLine += x
                x2sum += x * x
        
        xsum += xSumOfLine
        ysum += bitcountOfLine * y
        ye = y
        xe = xSumOfLine
        xysum += xe*ye
        y2sum += ye*ye*bitcountOfLine
        bitCount += bitcountOfLine

    # getMoments
    #double   x1, y1, x2, y2, xy
    if bitCount != 0:
        x2sum += 0.08333333 * bitCount
        y2sum += 0.08333333 * bitCount
        n = bitCount
        x1 = xsum/n
        y1 = ysum / n
        x2 = x2sum / n
        y2 = y2sum / n
        xy = xysum / n
        xm = x1
        ym = y1
        u20 = x2 - (x1 * x1)
        u02 = y2 - (y1 * y1)
        u11 = xy - x1 * y1

    # rest...
    m4 = 4.0 * np.abs(u02 * u20 - u11 * u11)
    if m4 < 0.000001:
        m4 = 0.000001
    a11 = u02 / m4
    a12 = u11 / m4
    a22 = u20 / m4
    xoffset = xm
    yoffset = ym

    tmp = a11 - a22
    if tmp == 0.0:
        tmp = 0.000001
    theta = 0.5 * np.arctan(2.0 * a12 / tmp)
    if theta < 0.0:
        theta += HALFPI
    if a12 > 0.0:
        theta += HALFPI
    elif a12 == 0.0:
        if a22 > a11:
            theta = 0.0
            tmp = a22
            a22 = a11
            a11 = tmp
        elif a11 != a22:
            theta = HALFPI
    tmp = np.sin(theta)
    if tmp == 0.0:
        tmp = 0.000001
    z = a12 * np.cos(theta) / tmp
    major = np.sqrt (1.0 / np.abs(a22 + z))
    minor = np.sqrt (1.0 / np.abs(a11 - z))
    scale = np.sqrt (bitCount / (np.pi * major * minor)) #equalize areas
    major = major*scale*2.0
    minor = minor*scale*2.0
    angle = 180.0 * theta / np.pi
    if angle == 180.0:
        angle = 0.0
    if major < minor:
        tmp = major
        major = minor
        minor = tmp
    xCenter = left + xoffset + 0.5
    yCenter = top + yoffset + 0.5
    
    if usePrint:
        print angle
        print major,minor
        print xCenter,yCenter
    
    return angle,major,minor,xCenter,yCenter

def EllipseFitterDNM(arr,cm=None):
    pts = np.where(arr==1)
    numPts = len(pts[0])
    
    if cm==None:
        cm=map(np.mean,pts)
    
    #cm -- from a python / matplotlib perspective, the 0,0 pixel is centered on 0,0
    [xCenter,yCenter] = [cm[0]+0.5 , cm[1]+0.5] # same meaning as IJ algorithm
    # from an ImageJ perspective, the pixel has its corner on 0,0
    # so that means it is centered on 0.5,0.5
    
    I_sq = np.zeros([2,2])
    I_cm = np.zeros([2,2])
    for i,p in enumerate(arr):
        for j,q in enumerate(arr):
            if arr[i,j]==1:
                # add up Ixx, Iyy, and Ixy
                # for this, m=1 and l=1 (mass and length of pixel square)
                Dx = i-cm[0]
                Dy = j-cm[1]
                I_cm += [[Dy**2, Dx*Dy],
                         [Dx*Dy,Dx**2]]
    I_sq += 1./12*numPts * np.identity(2)
    I_tot = I_cm + I_sq
    
    # Get Eigenvalues
    # I for ellipse with major and minor axis a,b is:
    # ((1/4 m b^2,0        ),
    #  (        0,1/4 m a^2))
    #I_adj = np.sqrt(I_tot*4./numPts) # radii
    #I_adj*2 # double for diameters...
    
    
    evals,evecs = np.linalg.eig(I_tot)
    major,minor = 2*np.sqrt(evals*4./numPts)
    print 'major,minor'
    major,minor
    
    Isqrt = 2*np.sqrt(I_tot*4./numPts)
    evalsB,evecsB = np.linalg.eig(Isqrt)
    majB,minB = evalsB
    print 'MM2'
    majB,minB
    # principle axis 1 angle
    an1 = np.arctan2(evecs[0][0],evecs[0][1])*180/np.pi
    an1 = 180+an1 if an1<0 else an1
    # principle axis 2 angle
    an2 = np.arctan2(evecs[1][0],evecs[1][1])*180/np.pi
    an2 = 180+an2 if an2<0 else an2
    # plot the major axis angle in 0 to 180 format
    angle = (an1 if evals[0]>evals[1] else an2)
    # now angle is the same as the IJ algorithm
    return [angle,major,minor,cm[1],cm[0]]

if CREATE_CYTHON_VERSION:
    exec(Cpyx.CythonInline('''
from __future__ import division
import numpy as np
cimport numpy as np
cimport cython

def EllipseFitterCython(np.ndarray[np.int_t, ndim=2] arr not None,usePrint=False):
    cdef int left=0 # holdover from the ImageJ version
    cdef int top=0 # holdover from the ImageJ version
    cdef int width=arr.shape[0]
    cdef int height=arr.shape[1]
    
    cdef double HALFPI = np.pi/2

    cdef double xCenter # X centroid
    cdef double yCenter # Y centroid
    cdef double major # Length of major axis
    cdef double minor # Length of minor axis
    cdef double angle # Angle in degrees
    cdef double theta # Angle in radians
    # cdef int[] xCoordinates # Initialized by makeRoi() # Not used
    # cdef int[] yCoordinates # Initialized by makeRoi() # Not used
    cdef int nCoordinates = 0 #int # Initialized by makeRoi()
    cdef int bitCount = 0 #int
    cdef double  xsum, ysum, x2sum, y2sum, xysum
    # cdef byte[] mask # Not used
    #cdef int left, top, width, height # already defined!
    cdef double   n
    cdef double   xm, ym   #mean values
    cdef double   u20, u02, u11  #central moments
    #ImageProcessor ip # Not used
    cdef double pw, ph
    cdef int record # bool

    cdef double sqrtPi = np.sqrt(np.pi)
    cdef double    a11, a12, a22, m4, z, scale, tmp, xoffset, yoffset
    cdef double    RealAngle

    #if mask==None: # mask should never be None
    #    major = (width*2) / sqrtPi
    #    minor = (height*2) / sqrtPi # * Info->PixelAspectRatio
    #    angle = 0.0
    #    theta = 0.0
    #    if major < minor:
    #        tmp = major
    #        major = minor
    #        minor = tmp
    #        angle = 90.0
    #        theta = np.pi/2.0
    #    xCenter = left + width / 2.0
    #    yCenter = top + height / 2.0
    #    return
    
    # computeSums
    xsum = 0.0
    ysum = 0.0
    x2sum = 0.0
    y2sum = 0.0
    xysum = 0.0
    cdef int bitcountOfLine
    cdef double   xe, ye
    cdef int xSumOfLine
    cdef int x,y
    for y in range(height):
        bitcountOfLine = 0
        xSumOfLine = 0
        #offset = y*width # int
        for x in range(width):
            if arr[x,y] != 0:
                bitcountOfLine+=1
                xSumOfLine += x
                x2sum += x * x
        
        xsum += xSumOfLine
        ysum += bitcountOfLine * y
        ye = y
        xe = xSumOfLine
        xysum += xe*ye
        y2sum += ye*ye*bitcountOfLine
        bitCount += bitcountOfLine

    # getMoments
    cdef double   x1, y1, x2, y2, xy
    if bitCount != 0:
        x2sum += 0.08333333 * bitCount
        y2sum += 0.08333333 * bitCount
        n = bitCount
        x1 = xsum/n
        y1 = ysum / n
        x2 = x2sum / n
        y2 = y2sum / n
        xy = xysum / n
        xm = x1
        ym = y1
        u20 = x2 - (x1 * x1)
        u02 = y2 - (y1 * y1)
        u11 = xy - x1 * y1

    # rest...
    m4 = 4.0 * np.abs(u02 * u20 - u11 * u11)
    if m4 < 0.000001:
        m4 = 0.000001
    a11 = u02 / m4
    a12 = u11 / m4
    a22 = u20 / m4
    xoffset = xm
    yoffset = ym

    tmp = a11 - a22
    if tmp == 0.0:
        tmp = 0.000001
    theta = 0.5 * np.arctan(2.0 * a12 / tmp)
    if theta < 0.0:
        theta += HALFPI
    if a12 > 0.0:
        theta += HALFPI
    elif a12 == 0.0:
        if a22 > a11:
            theta = 0.0
            tmp = a22
            a22 = a11
            a11 = tmp
        elif a11 != a22:
            theta = HALFPI
    tmp = np.sin(theta)
    if tmp == 0.0:
        tmp = 0.000001
    z = a12 * np.cos(theta) / tmp
    major = np.sqrt (1.0 / np.abs(a22 + z))
    minor = np.sqrt (1.0 / np.abs(a11 - z))
    scale = np.sqrt (bitCount / (np.pi * major * minor)) #equalize areas
    major = major*scale*2.0
    minor = minor*scale*2.0
    angle = 180.0 * theta / np.pi
    if angle == 180.0:
        angle = 0.0
    if major < minor:
        tmp = major
        major = minor
        minor = tmp
    xCenter = left + xoffset + 0.5
    yCenter = top + yoffset + 0.5
    
    if usePrint:
        print angle
        print major,minor
        print xCenter,yCenter
    
    return angle,major,minor,xCenter,yCenter
'''))

if __name__=='__main__':
    arr=np.array([[0,0,0,1,1,1,0,0,0],
                  [0,0,1,1,1,1,1,0,0],
                  [0,1,1,1,1,1,1,1,0],
                  [1,1,1,1,1,1,1,1,1],
                  [1,1,1,1,1,1,1,1,1],
                  [1,1,1,1,1,1,1,1,1],
                  [0,1,1,1,1,1,1,1,0],
                  [0,0,1,1,1,1,1,0,0],
                  [0,0,0,1,1,1,0,0,0]])
    import time
    t=time.time()
    a=EllipseFitter(arr,usePrint=False)
    print 'Py Time',time.time()-t ; t=time.time()
    if CREATE_CYTHON_VERSION:
        b=EllipseFitterCython(arr,usePrint=False)
        print 'Cy Time',time.time()-t ; t=time.time()# 3x speed up... not bad, but not wonderful...
    c=EllipseFitterDNM(arr)
    print 'DNM Time',time.time()-t ; t=time.time()
    # Test for consistency...
    print a
    #print b
    print c
    arr

def RotationM(theta):
    return np.array([[np.cos(theta),-np.sin(theta)],
                     [np.sin(theta), np.cos(theta)]]) # Forwards Rotation
def GetRotatedMajorMinor(major,minor,angle,woundAngle):
    th=(woundAngle - angle)*np.pi/180
    R=RotationM(th)
    Rinv=RotationM(-th)
    I = np.sqrt(np.dot(Rinv,np.dot([[major**2,0],
                                    [0,minor**2]],R)))
    return I[0][0],I[1][1]

# Still broken...
def drawEllipse(angle,major,minor,xCenter,yCenter,maxY):
    if major==0 and minor==0:
        return
    xc = int(np.round(xCenter))
    yc = int(np.round(yCenter))
    #int maxY = ip.getHeight()
    #int xmin, xmax
    #double sint, cost, rmajor2, rminor2, g11, g12, g22, k1, k2, k3
    #int x, xsave, ymin, ymax
    txmin = np.zeros(maxY)
    txmax = np.zeros(maxY)
    #double j1, j2, yr

    sint = np.sin(theta)
    cost = np.cos(theta)
    rmajor2 = 1.0 / (major/2)**2
    rminor2 = 1.0 / (minor/2)**2
    g11 = rmajor2 * (cost)**2 + rminor2 * (sint)**2
    g12 = (rmajor2 - rminor2) * sint * cost
    g22 = rmajor2 * sint**2 + rminor2 * cost**2
    k1 = -g12 / g11
    k2 = (g12**2 - g11 * g22) / g11**2
    k3 = 1.0 / g11;
    ymax = int(np.floor(np.sqrt(np.abs(k3 / k2))))
    if ymax>maxY:
        ymax = maxY
    if ymax<1:
        ymax = 1
    ymin = -ymax
    # Precalculation and use of symmetry speed things up
    for y in range(ymax):
        #GetMinMax(y, aMinMax);
        j2 = np.sqrt(k2 * sqr(y) + k3)
        j1 = k1 * y
        txmin[y] = int(np.round(j1 - j2))
        txmax[y] = int(np.round(j1 + j2))
    
    if record:
        xCoordinates[nCoordinates] = xc + txmin[ymax - 1]
        yCoordinates[nCoordinates] = yc + ymin
        nCoordinates+=1
    else:
        #DRAW HERE!!
        ip.moveTo(xc + txmin[ymax - 1], yc + ymin);
    for y in range(ymin,ymax):
        x = txmax[-y] if y<0 else -txmin[y];
        if record:
            xCoordinates[nCoordinates] = xc + x
            yCoordinates[nCoordinates] = yc + y
            nCoordinates+=1
        else:
            #DRAW HERE!!!
            ip.lineTo(xc + x, yc + y);
    
    for y in range(ymax,ymin,-1):
        x = txmin[-y] if y<0 else -txmax[y]
        if record:
            xCoordinates[nCoordinates] = xc + x
            yCoordinates[nCoordinates] = yc + y
            nCoordinates+=1
        else:
            #DRAW HERE!!!
            ip.lineTo(xc + x, yc + y);
