# Author: David Mashburn
# Created July 2006
# Last Modified December 2011
# License: (Apache 2, BSD, ???) -- whatever is most compatible with Cython...

# This module is for the automatic compilation (and also inlining) of
# Pyrex / Cython code...
# It can use distutils or manual compilation with gcc (or another compiler)
# It can work with a single existing C source and automatically compile it as well
# It has been tested on Windows, Mac, and Ubuntu Linux

# That said, I make no guarantees that it will work as expected!
# Numpy support is automatically enabled

# Unless the printCmds option is set to False, the script will output every action taken
# and command run

# My main goal for this is to aid people in learning how to compile cython code
# on their system, and give them a starting point so they can tweak what they want...

# My other goal is to automate the Cython compile process so I can do everything in
# one step after getting it set up.  I really like the inline feature a lot for testing!
# And try it with PySlices, the latest incarnation of the wxPython shell, PyCrust! (Shameless plug ;) )

import os
import sys
import glob
import random
import numpy

# Making this work in Vista and Windows 7...
# Download the latest mingw and add "C:\MinGW\bin" to the PATH environment variable

# Making this work on Mac...
# Download Xcode from the apple developer site (create a login) and install it:
# http://connect.apple.com

# Sample output for Cpyx on Windows w/Python2.5:
# Pieces:
# gcc -c -IC:/Python25/include PyrexExample.c -o PyrexExample.o
# gcc -shared PyrexExample.o -LC:/Python25/libs -lpython25 -o PyrexExample.pyd
# All-in-one:
# gcc -shared PyrexExample.c -IC:/Python25/include -LC:/Python25/libs -lpython25 -o PyrexExample.pyd
# All-in-one with linking dll...
# gcc -shared numpyTest.c -IC:/Python25/include -LC:/Python25/libs -LC:/Users/mashbudn/Programming/Python/Pyx -lpython25 -lnumpyTestC -o numpyTest.pyd

globalUseCython = True
globalUseDistutils = False

pythonName=sys.executable
verStr='.'.join([str(i) for i in sys.version_info[:2]]) # recently, 2.6 or 2.7
verStr2=''.join([str(i) for i in sys.version_info[:2]]) # recently, 26 or 27

if sys.platform=='win32':
    pyrexcName='"' + os.path.join(sys.prefix,'Scripts','pyrexc.py') + '"' # Full path to the Pyrex compiler script
    # this changed from previous versions...
    #cythonName='"' + os.path.join(sys.prefix,'Scripts','cython-script.py') + '"' # Full path to the Cython compiler script
    cythonName='"' + os.path.join(sys.prefix,'Lib','site-packages','cython.py') + '"' # Full path to the Cython compiler script
    pythonName=os.path.join(sys.prefix,'python.exe') # Full path to python.exe
    pythonInclude=os.path.join(sys.prefix,'include')
    pythonLibs=os.path.join(sys.prefix,'libs')
    #sitePackages=os.path.join(sys.prefix,'Lib','site-packages')
elif sys.platform=='darwin':
    pyrexcName='"' + os.path.join(sys.prefix,'bin','pyrexc') + '"' # Full path to the Pyrex compiler script
    cythonName='"' + os.path.join(sys.prefix,'bin','cython') + '"' # Full path to the Cython compiler script
    pythonName=os.path.join(sys.prefix,'bin','python') # Full path to python.exe
    pythonInclude=os.path.join(sys.prefix,'include')
    pythonLibs=os.path.join(sys.prefix,'lib','python'+verStr,'config/') # contains libpython2.5.so
    #sitePackages=os.path.join(sys.prefix,'lib','python'+verStr,'site-packages')
elif sys.platform=='linux2':
    pyrexcName='"' + os.path.join(sys.prefix,'bin','pyrexc') + '"' # Full path to the Pyrex compiler script
    cythonName='"' + os.path.join(sys.prefix,'bin','cython') + '"' # Full path to the Cython compiler script
    pythonName=os.path.join(sys.prefix,'bin','python'+verStr) # Full path to python.exe
    pythonInclude=os.path.join(sys.prefix,'include','python'+verStr)
    pythonLibs=os.path.join(sys.prefix,'lib') # contains libpython2.5.so
    #Debian / Ubuntu now put most everything in dist-packages...
    #Luckily, these variables are not actually used in this script
    #sitePackages=os.path.join(sys.prefix,'lib','python'+verStr,'site-packages')
    #distPackages=os.path.join(sys.prefix,'lib','python'+verStr,'dist-packages')
else:
    print 'Platform "' + sys.platform + '" not supported yet'

# Find numpy's arrayobject.h to include
arrayobjecthPath = os.path.join(numpy.get_include(),'numpy','arrayobject.h')
arrayObjectDir = numpy.get_include()

def ResolvePath(pth):
    '''Return the path if it is not an empty
    Otherwise try to use os.environ['MYPYREX']
    And fall back on the current directory'''
    if pth!='':
        return pth
    elif 'MYPYREX' in os.environ:
        return os.environ['MYPYREX']
    else:
        return os.getcwd()

def Cdll(cNameIn='',printCmds=True,gccOptions=''):
    '''Use gcc to compile a shared library (.dll on Windows, .so on *nix)'''
    cwd=os.getcwd()
    
    (cPath,cName)=os.path.split(cNameIn) # input path and input file name
    dllPath=cPath=ResolvePath(cPath)
    
    stripName=(os.path.splitext(cName))[0] # input file name without extension
    
    if sys.platform=='win32':      dllName='"' + os.path.join(dllPath,stripName+'.dll') + '"'
    elif sys.platform=='darwin':   dllName='"' + os.path.join(dllPath,'lib'+stripName+'.so') + '"'
    elif sys.platform=='linux2':   dllName='"' + os.path.join(dllPath,'lib'+stripName+'.so') + '"'
    else:                          print 'Platform "' + sys.platform + '" not supported yet'
    
    cName='"' + os.path.join(cPath,stripName+'.c') + '"' # redefine cName
    hName='"' + os.path.join(cPath,stripName+'.h') + '"'
    oName='"' + os.path.join(dllPath,stripName+'.o') + '"'
    
    os.chdir(cPath)
    
    cmd=' '.join(['gcc',gccOptions,'-fPIC','-c',cName,'-o',stripName+'.o'])
    if printCmds:
        print '\n', cmd
    os.system(cmd)
    
    cmd=' '.join(['gcc','-shared','-o',dllName,oName])
    if printCmds:
        print '\n', cmd
    os.system(cmd)
    
    os.chdir(cwd)

def Cpyx(pyxNameIn='CythonExample.pyx',useDistutils=globalUseDistutils,useCython=globalUseCython,gccOptions='',printCmds=True):
    '''Run Cython (or Pyrex) and then GCC to generate a Python extension module'''
    cwd=os.getcwd()
    
    (pyxPath,pyxName)=os.path.split(pyxNameIn) # input path and input file name
    
    pydPath=mainDir=ResolvePath(pyxPath)
    pyxStrip=(os.path.splitext(pyxName))[0] # input file name without extension
    
    extName='"' + pyxStrip + '"'
    pyxName='"' + os.path.join(mainDir,pyxStrip+'.pyx') + '"' # Full path to the PYX file (must be in Python/Pyx folder) redefine pyxName
    pyx2cName='"' + os.path.join(pydPath,pyxStrip+'.c') + '"' # Full path to the C file to be created
    pydName='"' + os.path.join(pydPath,pyxStrip+'.pyd') + '"' # Full path to the PYD file to be created
    soName='"' + os.path.join(pydPath,pyxStrip+'.so') + '"' # Full path to the lib*.so file to be created
    setupName='"' + os.path.join(pydPath,'setup.py') + '"' # Full path of the Setup File to be created
    
    if useDistutils:
        #write setup.py which will make a PYD file that can be imported
        setupText="""### This file is setup.py ###
from distutils.core import setup 
from distutils.extension import Extension 
from """ + ('Cython' if useCython else 'Pyrex') + """.Distutils import build_ext 

setup( 
  name = 'Lock module', 
  ext_modules=[ 
    Extension(""" + extName + ', [' + pyxName.replace('\\','\\\\') + ']' + """),
  ], 
  cmdclass = {'build_ext': build_ext}
)"""
        
        if printCmds:
            print 'Write Stuff to ', setupName[1:-1]
        fid = open(setupName[1:-1],'w') # [1:-1] removes quotes
        fid.write(setupText)
        fid.close()
        
        # run setup.py
        
        os.chdir(mainDir)
        
        if sys.platform=='win32':        cmd=' '.join([pythonName,setupName,'build_ext','--compiler=mingw32','--inplace'])
        elif sys.platform=='darwin':     cmd=' '.join([pythonName,setupName,'build_ext','--inplace'])
        elif sys.platform=='linux2':     cmd=' '.join([pythonName,setupName,'build_ext','--inplace'])
        else:                            print 'Platform "' + sys.platform + '" not supported yet'
    
    else:
        # run the main pyrex command to make the C file
        pyxCompiler = cythonName if useCython else pyrexcName
        cmd=' '.join([pythonName,pyxCompiler,pyxName,'-o',pyx2cName])
        if printCmds:
            print '\n', cmd
        os.system(cmd)        
        
        if sys.platform=='win32':        cmd=' '.join(['gcc',gccOptions,'-fPIC','-shared',pyx2cName,'-I'+pythonInclude,'-I'+arrayObjectDir,'-L'+pythonLibs,'-lpython'+verStr2,'-o',pydName])
        elif sys.platform=='darwin':     cmd=' '.join(['gcc',gccOptions,'-fno-strict-aliasing','-Wno-long-double','-no-cpp-precomp','-mno-fused-madd','-fno-common',
                                                       '-dynamic','-DNDEBUG','-g','-O3','-bundle','-undefined dynamic_lookup','-I'+pythonInclude,
                                                       '-I'+pythonInclude+'/python'+verStr,'-I'+arrayObjectDir,'-L'+pythonLibs,'-L/usr/local/lib',pyx2cName,'-o',soName])
        elif sys.platform=='linux2':     cmd=' '.join(['gcc',gccOptions,'-fPIC','-shared',pyx2cName,'-I'+pythonInclude,'-L'+pythonLibs,'-lpython'+verStr,'-o',soName])
        else:                            print 'Platform "' + sys.platform + '" not supported yet'
    
    if printCmds:
        print '\n', cmd
    os.system(cmd)
    
    os.chdir(cwd)

def CpyxLib(pyxNameIn='CythonExample.pyx',cNameIn='CTestC.c',recompile=True,useDistutils=globalUseDistutils,useCython=globalUseCython,gccOptions='',printCmds=True):
    '''Compile a C source and then run Cython/GCC on a pyx file, linking the C source
    Alternatively, use Cython.Distutils to do the same thing'''
    cwd=os.getcwd()
    
    (pyxPath,pyxName)=os.path.split(pyxNameIn) # input path and input file name
    (cPath,cName)=os.path.split(cNameIn) # input path and input file name
    
    pydPath=mainDir=ResolvePath(pyxPath)
    dllPath=cPath=ResolvePath(cPath)
    
    pyxStrip=(os.path.splitext(pyxName))[0] # input file name without extension
    cStrip=(os.path.splitext(cName))[0] # input file name without extension
    
    extName='"' + pyxStrip + '"'
    pyxName='"' + os.path.join(mainDir,pyxStrip+'.pyx') + '"' # Full path to the PYX file (must be in Python\\Pyrex folder)
    
    if sys.platform=='win32':
        dllName='"' + os.path.join(dllPath,cStrip+'.dll') + '"'
        libName='"' + os.path.join(dllPath,cStrip) + '"'
        library_dirs_txt=''
    elif sys.platform=='darwin':
        dllName='"' + os.path.join(dllPath,'lib'+cStrip+'.so') + '"'
        libName='"' + cStrip + '"'
        library_dirs_txt="""
    library_dirs=[""" + '"' + pydPath.replace('\\','\\\\') + '"' + """],
    runtime_library_dirs=[""" + '"' + pydPath.replace('\\','\\\\') + '"' + """],"""
    elif sys.platform=='linux2':
        dllName='"' + os.path.join(dllPath,'lib'+cStrip+'.so') + '"'
        libName='"' + cStrip + '"'
        library_dirs_txt="""
    library_dirs=[""" + '"' + pydPath.replace('\\','\\\\') + '"' + """],
    runtime_library_dirs=[""" + '"' + pydPath.replace('\\','\\\\') + '"' + """],"""
    else:
        print 'Platform "' + sys.platform + '" not supported yet'
    
    cName='"' + os.path.join(cPath,cStrip+'.c') + '"' # redefine cName
    hName='"' + os.path.join(cPath,cStrip+'.h') + '"'
    oName='"' + os.path.join(dllPath,cStrip+'.o') + '"'
    
    os.chdir(cPath)
    
    pyx2cName='"' + os.path.join(pydPath,pyxStrip+'.c') + '"' # Full path to the C file to be created
    setupName='"' + os.path.join(pydPath,'setup.py') + '"' # Full path of the Setup File to be created
    pydName='"' + os.path.join(pydPath,pyxStrip+'.pyd') + '"' # Full path to the PYD file to be created
    soName='"' + os.path.join(pydPath,pyxStrip+'.so') + '"' # Full path to the lib*.so file to be created
    
    if useDistutils:
        #write setup.py which will make a PYD file that can be imported
        setupText="""### This file is setup.py ###
from distutils.core import setup 
from distutils.extension import Extension 
from """+('Cython' if useCython else 'Pyrex')+""".Distutils import build_ext 

setup( 
  name = 'Lock module', 
  ext_modules=[ 
    Extension(""" + extName + ', [' + pyxName.replace('\\','\\\\') + """],"""+library_dirs_txt+"""
    libraries=[""" + libName.replace('\\','\\\\') + ']' + """),
  ], 
  cmdclass = {'build_ext': build_ext} 
)"""
        if printCmds:
            print 'Write Stuff to ', setupName[1:-1]
        fid = open(setupName[1:-1],'w') # [1:-1] removes quotes
        fid.write(setupText)
        fid.close()
        
        # run setup.py
        os.chdir(mainDir)
        
        if sys.platform=='win32':        cmd=' '.join([pythonName,setupName,'build_ext','--compiler=mingw32','--inplace'])
        elif sys.platform=='darwin':     cmd=' '.join([pythonName,setupName,'build_ext','--inplace'])
        elif sys.platform=='linux2':     cmd=' '.join([pythonName,setupName,'build_ext','--inplace'])
        else:                            print 'Platform "' + sys.platform + '" not supported yet'
    else:
        # compile the DLL needed for the link to the C file
        if recompile:
            Cdll(cName[1:-1],printCmds=printCmds,gccOptions=gccOptions) # [1:-1] to remove the quotes
        
        # run the main pyrex command to make the C file
        pyxCompiler = cythonName if useCython else pyrexcName
        cmd=' '.join([pythonName,pyxCompiler,pyxName,'-o',pyx2cName])
        if printCmds:
            print '\n', cmd
        os.system(cmd)        
        
        if sys.platform=='win32':        cmd=' '.join(['gcc',gccOptions,'-fPIC','-shared',pyx2cName,'-I'+pythonInclude,'-L'+pythonLibs,'-L'+cPath,'-Wl,-R'+cPath,'-lpython'+verStr2,'-l'+cStrip,'-o',pydName])
        elif sys.platform=='darwin':     cmd=' '.join(['gcc',gccOptions,'-fno-strict-aliasing','-Wno-long-double','-no-cpp-precomp','-mno-fused-madd','-fno-common',
                                                       '-dynamic','-DNDEBUG','-g','-O3','-bundle','-undefined dynamic_lookup','-I'+pythonInclude,
                                                       '-I'+pythonInclude+'/python'+verStr,'-I'+arrayObjectDir,'-L'+pythonLibs,'-L/usr/local/lib','-L'+cPath,'-Wl,-R'+cPath,
                                                       '-l'+cStrip,pyx2cName,'-o',soName])
        elif sys.platform=='linux2':     cmd=' '.join(['gcc',gccOptions,'-fPIC','-shared',pyx2cName,'-I'+pythonInclude,'-L'+pythonLibs,'-L'+cPath,'-Wl,-R'+cPath,'-lpython'+verStr,'-l'+cStrip,'-o',soName])
        else:                            print 'Platform "' + sys.platform + '" not supported yet'
    
    if printCmds:
        print '\n', cmd
    os.system(cmd)
    
    os.chdir(cwd)

# Currently unused...
def InlineTableCreate():
    tmpDir = os.path.expanduser('~/.Cpyx_tmp')
    tabFile = os.path.join(tmpDir,'InlineTable.txt')
    if not os.path.exists(tmpDir):
        os.mkdir(tmpDir)
    if not os.path.exists(tabFile):
        fid=open(tabFile,'w')
        fid.write('# This is the Table of .pyx files created for Inline use)\n' + \
                  '# It has entries like: module_name (w/o .pyx),num_lines,num_characters,md5sum\n')
        fid.close()

# Currently unused...
def InlineTableFind(codeIn):
    import hashlib
    code = codeIn.replace(os.linesep,'\n')
    entry = [len(code),len(code.split('\n')),hashlib.md5(code).hexdigest()]
    tmpDir = os.path.expanduser('~/.Cpyx_tmp')
    tabFile = os.path.join(tmpDir,'InlineTable.txt')
    
    if os.path.exists(tabFile):
        out=-1
        fid=open(tabFile,'r')
        for i,s in enumerate(fid):
            if i>2:
                l=s.split(',')
                l=[l[0],l[1],l[2],l[3]]
                if l[1:]==entry:
                    file=os.path.join(tmpDir,l[0]+'.pyx')
                    if os.path.exists(file):
                        fid2=open(file,'r')
                        if code==fid2.read():
                            out = file
                        fid2.close()
                    break
        fid.close()
        return out
    else:
        return None

# Currently unused...
def InlineTableAddEntry(codeIn):
    val=InlineTableFind(codeIn)
    if val==None:
        InlineTableCreate()
    
    # Now Table should definitely exist...
    if val==-1:
        code = codeIn.replace(os.linesep,'\n')
        entry = [len(code),len(code.split('\n')),hashlib.md5(code).hexdigest()]
        tmpDir = os.path.expanduser('~/.Cpyx_tmp')
        tabFile = os.path.join(tmpDir,'InlineTable.txt')
        
        moduleName='Pyrex'+str(random.randint(0,1e18))
        
        fid = open(tabFile,'a')
        fid.write(moduleName+','+str(entry[0])+','+str(entry[1])+','+entry[2]+'\n')
        fid.close()
        
        return moduleName
    else:
        return None

# Shamelessly steal the idea used by scipy.weave.inline but for Pyrex/Cython instead...
# In order to be able to import *, you have to use exec in the calling module...
def CythonInline(code,cleanUp=False,useDistutils=globalUseDistutils,useCython=True,gccOptions='',printCmds=True):
    '''Use Cpyx to compile Cython code inline
    CythonInline returns a string that is an import statement to the temporary Cython module
    Use this like: exec(CythonInline(r"""<somecode>""",<options>))'''
    
    testCode=r"""
cdef extern from "stdio.h":
    ctypedef struct FILE

    FILE * stdout
    int printf(char *format,...)
    int fflush( FILE *stream )

def CythonPrint(mystring):
    printf(mystring)
    fflush(stdout)

CythonPrint('HelloWorld!')
"""
    tmpPath=os.path.expanduser('~/.Cpyx_tmp')
    if not os.path.isdir(tmpPath):
        os.mkdir(tmpPath)
    if tmpPath not in sys.path:
        sys.path.append(tmpPath)
    if cleanUp:
        CleanTmp()
    
    # Ensure you always get a new module!
    # This means there is no reason to "reload"
    # Also means memory gets majorly eaten up!
    # Can't have everything!
    moduleName='Cython'+str(random.randint(0,1e18))
    file=os.path.join(tmpPath,moduleName+'.pyx')
    
    fid=open(file,'w')
    fid.write(code)
    fid.close()
    
    Cpyx(file,useDistutils=useDistutils,useCython=useCython,gccOptions=gccOptions,printCmds=printCmds)
    
    #cmd="""import """+moduleName+""" as LoadPyrexInline"""
    cmd="""from """+moduleName+""" import *"""
    if printCmds:
        print cmd
    return cmd

# Create a dummy function that defaults to using Pyrex instead for clarity...
def PyrexInline(code,cleanUp=False,useDistutils=globalUseDistutils,useCython=False,gccOptions='',printCmds=True):
    return CythonInline(code,cleanUp=cleanUp,useDistutils=useDistutils,useCython=useCython,gccOptions=gccOptions,printCmds=printCmds)
PyrexInline.__doc__ = CythonInline.__doc__.replace('Cython','Pyrex')

def CleanTmp():
    '''Remove temp files and built modules from '~/.Cpyx_tmp' '''
    tmpPath=os.path.expanduser('~/.Cpyx_tmp')
    for i in glob.glob(os.path.join(tmpPath,'*')):
        os.remove(i)
