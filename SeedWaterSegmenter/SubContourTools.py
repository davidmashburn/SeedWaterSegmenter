import glob
from copy import copy,deepcopy
import wx
from cmpGen import cmpGen
import ExcelHelper
from ValueReceived import imshow_vr
from np_utils import flatten,totuple,removeDuplicates,deletecases,floatIntStringOrNone
#from np_utils:
from np_utils import limitInteriorPoints,limitInteriorPointsInterpolating

def GetValuesAroundSCPoint(watershed2d,point):
    x,y = point
    if  0 < x <= watershed2d.shape[0] and  0 < y <= watershed2d.shape[1]:
        return tuple(np.unique( watershed2d[x-1:x+1,y-1:y+1] ).tolist())
    else:
        print "THIS POINT IS NOT INTERIOR TO THE ARRAY; THERE ARE NOT 4 POINTS AROUND IT!"
        return (None,None,None)


class SubContour:
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
    def cVLS(self): # for legacy purposes, returns a list
        return [self.values,self.adjusted_length,self.points]

def SubContourListfromCVLSList(cVLS_list,startPointValues_List=[],endPointValues_List=[]):
    if startPointValues_List==[]:
        startPointValues_List = [[None,None,None] for c in cVLS_List]
    if endPointValues_List==[]:
       endPointValues_List = [[None,None,None] for c in cVLS_List]
    return [ SubContour(points = cvls[2],
                        numPoints = len(cvls[2]),
                        adjusted_length = cvls[1],
                        values = tuple(cvls[1]),
                        startPointValues = tuple(startPointValues_List[i]),
                        endPointValues = tuple(endPointValues_List[i]))
            for i,cvls in enumerate(cVLS_List)]

def SavexyListAndPolyListToMMAFormat(xyList,polyList,filename,bumpIndsUp1=True,removeLastPoint=True):
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

