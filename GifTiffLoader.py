"""GifTiffLoader is a module to automatically load Tiff and Gif files as
numpy arrays using PIL.  Specifically it is designed to deal with 8-bit
and 16-bit images typically used as the standard formats for microscope
images.  It is also designed to deal with image stacks (Animated Gifs and
Multi-page Tiffs) and image sequences (both 3D and 4D sequences).

GifTiffLoader relies on two other modules, FilenameSort and cmpGen.

Sponsored by the NSF and HFSP through the Shane Hutson Laboratory, part of
Vanderbilt Institute of Integrative Biosystems Research (VIIBRE)."""

__author__ = "David N. Mashburn <david.n.mashburn@gmail.com>"

import os
import numpy as np
import Image
import wx

from FilenameSort import cmp_fnames,cmp_fnames_A,getSortedListOfFiles, \
                         getSortedListOfNumericalEquivalentFiles

# This module can load and save .tiff and .gif formats as numpy arrays...

_gif_names=['gif','GIF','Gif']
_tif_names=['tif','TIF','Tif','tiff','TIFF','Tiff']

def DivideConvertType(arr,bits=16,maxVal=None,zeroMode='clip',maxMode='clip'):
    #Get datatype from bits...
    # Note that this function is only useful for unsigned ints!
    type_dict={8:np.uint8,16:np.uint16,32:np.uint32,64:np.uint64}
    
    if bits not in type_dict.keys():
        print "'bits' must be one of 8, 16, 32, 64"
        return
    
    maxClip = 2**bits-1
    if maxVal==None:
        maxVal=arr.max()
    elif type(maxVal)==int:
        if maxMode=='clip':
            arr.clip(arr.min(),maxVal,out=arr) # in place
        elif maxMode=='stretch':
            # Now maxVal is used at the end instead...
            maxClip = maxVal
        else:
            print "'maxMode' must be either 'clip' or 'stretch'"
    else:
        print "'maxVal' must be given as an integer (or None)!"
        return
    
    minVal=0
    
    if zeroMode=='clip':
        arr.clip(0,arr.max(),out=arr) # in place
    elif zeroMode=='abs':
        np.absolute(arr,arr) # in place
    elif zeroMode=='stretch':
        minVal=arr.min()
    else:
        print "'zeroMode' must be one of: 'clip' , 'abs' , 'stretch'"
    
    clipRatio = 1.*maxClip/(maxVal-minVal)
    
    if zeroMode=='stretch' or clipRatio<1:
        arr = np.double(arr)
        arr -= minVal # separate step cuts down on memory usage...
        if clipRatio<1:
            arr *= clipRatio
    
    return np.array(arr,dtype=type_dict[bits])

def ConvertTo8Bit(arr):
    return DivideConvertType(arr,8)    
def ConvertTo16Bit(arr):
    return DivideConvertType(arr,16)
def ConvertTo32Bit(arr):
    return arr.asarray(np.float32) #DivideConvertType(t,32) # b/c we want 32 bit float, not uint

def GetShape(filename=None):
    if filename==None:
        filename=wx.FileSelector()
        print filename
    im=Image.open(filename)
    
    if im.format in _tif_names:  datatype=np.uint16
    elif im.format in _gif_names: datatype=np.uint8
    
    h,w=im.size
    
    numFrames=0
    while 1:
        numFrames+=1
        try:
            im.seek(numFrames)
        except EOFError:
            break # end of sequence
    
    return numFrames,w,h

def GetDirectoryShape(dirname=None,wildcard='*[!.txt]',dims=1):
    '''For data stored as image sequences, get the number of them'''
    if dirname==None:
        dirname=wx.DirSelector()
        print dirname
    
    files = getSortedListOfFiles(dirname,globArg=wildcard)
    
    if len(files)==0:
        print 'Empty Directory!'
        return 0,[]
    
    if dims==1: # flat sequence (time series)
        return len(files),files
    elif dims==2: # series of z-stacks
        l=[]
        
        for i in files:
            e=os.path.split(i)[1]
            if len(e.split('_'))>=3:
                #print e,':', e.split('_')[0],e.split('_')[1],e.split('_')[2][:3]
                #l.append([i,int(e.split('_')[1]),int(e.split('_')[2][:3])])
                l.append([i,int(e.split('_')[-2]),int(e.split('_')[-1][:3])])
        
        if len(l)==0:
            print 'Empty Directory!'
            return [0,0],[]
        
        tdim=max([i[1] for i in l])+1
        zdim=max([i[2] for i in l])+1
        return [tdim,zdim],l
    else:
        print 'No code for 5D stacks... yet...'
        return # Should throw downstream errors... just sayin...

def GetDatatype(im):
    """Automatically determine the datatype using the PIL -> numpy conversion."""
    
    # Old way:
    #if im.format in _tif_names:  datatype=np.uint16
    #elif im.format in _gif_names: datatype=np.uint8
    
    t=np.array(im)
    
    # Work around for windows:
    if t.dtype==np.object:
        if im.mode=='L': # 8-bit gray scale
            datatype = np.uint8
        elif im.mode=='I;16': # 16-bit gray scale
            datatype = np.uint16
        elif im.mode in ['I','F']: # 32-bit gray scale (int and float)
            print '32-bit images not supported at this time!'
        else:
            print 'Unknown image type!  Assume to 16-bit...'
            datatype = np.uint16
    elif t.dtype in [np.uint8,np.uint16,np.float32]:
        datatype = t.dtype
    elif t.dtype==np.int16:
        datatype = np.uint16 # Convert from int16 to uint16 (same as what ImageJ does [I think])
    else:
        print 'What kind of image format is this anyway???'
        print 'Should be a 16 or 32 bit tiff or 8 bit unsigned gif or tiff'
        return
    
    return datatype

def GetShapeMonolithicOrSequence(filename=None):
    if filename==None:
        filename=wx.FileSelector()
        print filename
    
    numFrames,w,h = GetShape(filename)
    
    isSequence = False
    
    if numFrames == 1:
        numFrames = len(getSortedListOfNumericalEquivalentFiles(
                                    filename,os.path.split(filename)[0]
                                                               ))
        isSequence = True
    return numFrames,w,h,isSequence

def LoadSingle(filename=None):
    if filename==None:
        filename=wx.FileSelector()
        print filename
    im=Image.open(filename)
    
    # TODO: Need to add more smarts to this based on the mode 'I;16' vs. RGB, etc...
    #'1','P','RGB','RGBA','L','F','CMYK','YCbCr','I',
    # Can I simplify this by simply fixing the PIL -> numpy conversion?
    
    #if im.format in _tif_names:  datatype=np.uint16
    #elif im.format in _gif_names: datatype=np.uint8
    
    # Try to automatically determine the datatype using the
    # PIL -> numpy conversion.
    t=np.array(im)
    
    # Work around for windows:
    if t.dtype==np.object:
        if im.mode=='L': # 8-bit gray scale
            t=np.array(im.getdata(),np.uint8)
        elif im.mode=='I;16': # 16-bit gray scale
            t=np.array(im.getdata(),np.uint16)
        elif im.mode in ['I','F']: # 32-bit gray scale (int and float)
            print '32-bit images not supported at this time!'
        else:
            print 'Unknown image type!  Assume to 16-bit...'
            t=np.array(im.getdata(),np.uint16)
    
    if t.dtype in [np.uint8,np.uint16,np.float32]:
        pass
    elif t.dtype==np.int16:
        t=np.uint16(t) # Convert from int16 to uint16 (same as what ImageJ does [I think])
    else:
        print 'What kind of image format is this anyway???'
        print 'Should be a 16 or 32 bit tiff or 8 bit unsigned gif or tiff'
        return
    
    t.resize(im.size[1],im.size[0])
    
    # Patch around 16-bit endianness bug for PIL before version 1.1.7...
    if im.format in _tif_names:
        fid=open(filename,'rb')
        endianness=fid.read(2)
        fid.close()
        if endianness=='MM':
            # I assume no one will use earlier than 1.1.4
            # 1.1.4 to 1.1.6 had a bug; fixed in 1.1.7 and later
            if Image.VERSION in ['1.1.4','1.1.5','1.1.6']:
                t.byteswap(True)
        elif endianness=='II':
            pass
        else:
            print "Just so you know... that ain't no tiff!"
    
    return t

# I also need to check for filenames with the same base in the selected folder
# And then give a "stupidity check" for overwriting
# Same for Sequence and monolithic...
# PIL saved images look stupid...  test manually again
# PIL is too dubm to figure out TIF instead of TIFF -- Fix this...
def SaveSingle(arr,filename=None,format='gif',tiffBits=16):
    if format in _gif_names:
        formatName='Gif'
        datatype=np.uint8
        bits=8
    elif format in _tif_names:
        formatName='Tiff'
        if tiffBits==8:
            datatype=np.uint8
            bits=8
        elif tiffBits==16:
            datatype=np.uint16
            bits=16
        elif tiffBits==32:
            datatype=np.float32
            bits=32
        else:
            print 'Unsupported tiff format!  Choose 8, 16, or 32 bit.'
            return
    else:
        print 'Unsupported Format!  Choose Gif or Tiff format'
        return
    
    if filename==None:
        filename=wx.SaveFileSelector('Array as '+formatName,'.'+format)
    
    # Make the directory if it does not exist
    if not os.path.exists(os.path.split(filename)[0]):
        os.mkdir(os.path.split(filename)[0])
    
    if datatype!=arr.dtype:
        if datatype in [np.uint8, np.uint16]:
            arr=DivideConvertType(arr,bits=bits)
        elif datatype==np.float32:
            arr = np.array(arr,np.float32)
    if datatype==np.uint8:
        im=Image.fromarray(arr)
    elif datatype==np.uint16:
        im=Image.fromarray(arr,'I;16') # Should I ever use I;16B ???
    elif datatype==np.float32:
        im=Image.fromarray(arr,'F')
    else:
        print 'Unsupported tiff format!  Choose 8, 16, or 32 bit.'
    im.save(filename,formatName)

def LoadFileSequence(dirname=None,wildcard='*[!.txt]'):
    if dirname==None:
        dirname=wx.DirSelector()
        print dirname
    
    numFiles,files = GetDirectoryShape(dirname,wildcard)
    
    if numFiles==0:
        return
    
    t0=LoadSingle(files[0])
    
    t=np.zeros([numFiles,t0.shape[0],t0.shape[1]],t0.dtype)
    
    for i in range(numFiles):
        t[i]=LoadSingle(files[i])
    
    return t

def SaveFileSequence(arr,basename=None,format='gif',tiffBits=16,sparseSave=None,functionToRunOnFrames=None):
    '''Save a sequence of files with the format _t_zzz where t is the
stack number and zzz is the frame number.
"tiffBits" should be 8 or 16.
"sparseSave" should either be None (the default) or a boolean list the same shape as the z-t shape of the array.
If specified, it determines which files should be saved and which should not.
It is basically a way to save only SOME of the files in the sequence.'''
    if basename==None:
        basename=wx.SaveFileSelector('Array as Sequence of Gifs -- numbers automatically added to','.'+format)
    
    if format in _gif_names:
        pass
    elif format in _tif_names:
        if tiffBits not in [8,16,32]:
            print 'Unsupported tiff format!  Choose 8 or 16 bit.'
            return
    else:
        print 'Unsupported Format!  Choose Gif or Tiff format'
        return
    
    if arr.dtype==np.uint8 and arr.max()>255:
        arr=np.array(arr*255./arr.max(),np.uint8)
    
    if functionToRunOnFrames==None:
        functionToRunOnFrames=( lambda x:x )
    
    if len(arr.shape)==2:
        SaveSingle(functionToRunOnFrames(arr),os.path.splitext(basename)[0]+'_0_000'+'.'+format,format=format,tiffBits=tiffBits)
    elif len(arr.shape)==3:
        if sparseSave==None:
            sparseSave=[True]*arr.shape[0]
        for i in range(arr.shape[0]):
            if sparseSave[i]:
                SaveSingle(functionToRunOnFrames(arr[i]),os.path.splitext(basename)[0]+'_0_'+str('%03d' % i)+'.'+format,format=format,tiffBits=tiffBits)
    elif len(arr.shape)==4:
        if sparseSave==None:
            sparseSave=[True]*arr.shape[0]
            for i in range(arr.shape[0]):
                sparseSave[i]=[True]*arr.shape[1]
        for i in range(arr.shape[0]):
            for j in range(arr.shape[1]):
                if sparseSave[i][j]:
                    SaveSingle(functionToRunOnFrames(arr[i,j]),os.path.splitext(basename)[0]+'_'+str(i)+'_'+str('%03d' % j)+'.'+format,format=format,tiffBits=tiffBits)
    else:
        print 'This function does not support saving arrays with more than 4 dimensions!'

SaveFileSequence4D=SaveFileSequence

def LoadMonolithic(filename=None):
    #f='C:/Documents and Settings/Owner/Desktop/Stack_Zproject_GBR_DC.gif'
    if filename==None:
        filename=wx.FileSelector()
        print filename
    im=Image.open(filename)
    
    datatype = GetDatatype(im)
    
    l=[]
    numFrames=0
    while 1:
        l.append(np.asarray(im.getdata(),dtype=datatype))
        numFrames+=1
        try:
            im.seek(im.tell()+1)
        except EOFError:
            break # end of sequence
    
    t=np.array(l,dtype=datatype)
    
    t.resize(numFrames,im.size[1],im.size[0])
    
    # Patch around 16-bit endianness bug for PIL before version 1.1.7...
    if im.format=='TIFF':
        fid=open(filename,'rb')
        endianness=fid.read(2)
        fid.close()
        if endianness=='MM':
            # I assume no one will use earlier than 1.1.4
            # 1.1.4 to 1.1.6 had a bug; fixed in 1.1.7 and later
            if Image.VERSION in ['1.1.4','1.1.5','1.1.6']:
                t.byteswap(True)
        elif endianness=='II':
            pass
        else:
            print "Just so you know... that ain't no tiff!"
    
    return t

#NOT TESTED YET!!!
def LoadFrameFromMonolithic(filename=None,frameNum=0):
    #f='C:/Documents and Settings/Owner/Desktop/Stack_Zproject_GBR_DC.gif'
    if filename==None:
        filename=wx.FileSelector()
        print filename
    im=Image.open(filename)
    
    datatype = GetDatatype(im)
    
    t=None
    i=0
    while 1:
        if i==frameNum:
            t = np.asarray(im.getdata(),dtype=datatype)
            t.resize(numFrames,im.size[1],im.size[0])
            break
        i+=1
        try:
            im.seek(im.tell()+1)
        except EOFError:
            break # end of sequence
    
    if t==None:
        print 'Frame not able to be loaded'
        return t
    else:
        # Patch around 16-bit endianness bug for PIL before version 1.1.7...
        if im.format=='TIFF':
            fid=open(filename,'rb')
            endianness=fid.read(2)
            fid.close()
            if endianness=='MM':
                # I assume no one will use earlier than 1.1.4
                # 1.1.4 to 1.1.6 had a bug; fixed in 1.1.7 and later
                if Image.VERSION in ['1.1.4','1.1.5','1.1.6']:
                    t.byteswap(True)
            elif endianness=='II':
                pass
            else:
                print "Just so you know... that ain't no tiff!"
        
        return t

def SaveMonolithic(arr,filename=None):
    pass #Coming soon! -- see Examples/PIL Examples/gifmaker
    # Can I do monolithic Tiff, too??

def LoadMonolithicOrSequenceSpecial(filename=None):
    if filename==None:
        filename=wx.FileSelector()
        print filename
    
    numFrames,w,h,isSequence = GetShapeMonolithicOrSequence(filename)
    
    if isSequence:
        files = getSortedListOfNumericalEquivalentFiles(
                                    filename,os.path.split(filename)[0])
        if len(files)==0:
            print 'Empty Directory!'
            return
        t0=LoadSingle(files[0])
        
        t=np.zeros([len(files),t0.shape[0],t0.shape[1]],t0.dtype)
        
        for i in range(len(files)):
            t[i]=LoadSingle(files[i])
        
        return t
    else:
        return LoadMonolithic(filename)

def LoadSequence4D(dirname=None,wildcard='*[!.txt]'):
    if dirname==None:
        dirname=wx.DirSelector()
        print dirname
    
    [tdim,zdim],files = GetDirectoryShape(dirname,wildcard,dims=2)
    
    if len(files)==0:
        print 'Empty Directory!'
        return
    
    test=LoadSingle(files[0][0])
    
    [xdim,ydim]=test.shape
    
    t=np.zeros([tdim,zdim,xdim,ydim],test.dtype)
    
    for i in files:
        t[i[1],i[2]]=LoadSingle(i[0]) #load the data for each gif directly into one big array, indexing based on numbering
    
    return t

def LoadMonolithicSequence4D(dirname=None,wildcard='*[!.txt]'):
    if dirname==None:
        dirname=wx.DirSelector()
        print dirname
    files = getSortedListOfFiles(dirname,globArg=wildcard)
    if len(files)==0:
        print 'Empty Directory!'
        return
    t0=LoadMonolithic(files[0])
    
    t=np.zeros([len(files)]+list(t0.shape),t0.dtype)
    
    for i in range(len(files)):
        t[i]=LoadMonolithic(files[i])
    
    return t

def LoadGroupedZCroppedByTxtInput(dirname,StartStopTxt,mergeOperation=None):
    """StartStopTxt gives the range to load for each stack... "0: 4-9\n1: 5-10\n2: 5-11"
The stack number before the colon is NOT ACTUALLY USED, stacks are just processed in order.
For mergeOperation, select 'None','mean','max','min',or 'sum'"""
    StartStop = [map(int,i.replace(' ','').split(':')[1].split('-')) for i in StartStopTxt.split('\n')]
    stacks = []
    for i,p in enumerate(StartStop):
        stacks.append([])
        files = getSortedListOfFiles(dirname,globArg="*_"+str(i)+"_*PMT1.TIF")
        files = files[p[0]:p[1]+1] # only pick out some files
        for f in files:
            stacks[-1].append(LoadSingle(f))
        stacks[-1]=np.array(stacks[-1])
    
    if mergeOperation in ['mean','max','min','sum']:
        return np.array([i.__getattribute__(mergeOperation)(axis=0) for i in stacks])
    elif mergeOperation!=None:
        print 'Invalid mergeOperation!  Returning the stacks list instead'
    
    return stacks # has to be a list not an array b/c elements might not be the same length
        

if __name__=='__main__':
    app=wx.App(0)
    t=LoadFileSequence()
    print 'Loaded an array of type',t.dtype,'with shape',t.shape
    
    import DelayApp
