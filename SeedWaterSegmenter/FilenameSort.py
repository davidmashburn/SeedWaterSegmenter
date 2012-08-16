"""FilenameSort is a utility to aid in "human-like" sorting of file names.

Normally using sort, ["file_1_10a.png","file_1_1a.png","file_1_5a.png"] would sort as:
["file_1_10a.png","file_1_1a.png","file_1_5a.png"]

Using the function getSortableList instead results in:
["file_1_1a.png","file_1_5a.png","file_1_10a.png"]

Which is more like what one would expect.

FilenameSort uses cmpGen to aid sorting."""

__author__ = "David N. Mashburn <david.n.mashburn@gmail.com>"

import numpy as np
import glob,copy,os
from time import time
from cmpGen import cmpGen

# The two dinky little functions sort files according to human-based sorting...
# They avoid the usual sorting problem with non-0 buffered integers in names...:
# 0a,1a,10a,11a,12a,13a,14a,15a,16a,17a,18a,19a,2a,20a ...
# with this function sorts properly (like follows) instead:
# 0a,1a,2a,3a,4a,5a,6a,7a,8a,9a,10a,11a,12a,13a,14a,15a,16a,17a,18a,19a,20a ...
def getSortableList(s):
    '''Turns a string into a list where numbers and non-numbers are separated'''
    currentNumber=[]
    currentNonNumber=[]
    sortableList=[os.path.split(s)[0]]
    for i in os.path.splitext(os.path.split(s)[1])[0]:
        if i.isdigit():
            currentNumber.append(i)
            if currentNonNumber!=[]:
                sortableList.append(''.join(currentNonNumber))
                currentNonNumber=[]
        else:
            currentNonNumber.append(i)
            if currentNumber!=[]:
                sortableList.append(int(''.join(currentNumber)))
                currentNumber=[]
    if currentNumber!=[]:
        sortableList.append(int(''.join(currentNumber)))
    elif currentNonNumber!=[]:
        sortableList.append(''.join(currentNonNumber))
    sortableList.append(os.path.splitext(s)[1])
    return sortableList

def getSortedListOfFiles(d,globArg='*[!.txt]'):
    files = glob.glob(os.path.join(d,globArg))
    
    files.sort(cmpGen(getSortableList))
    
    return files

def AreFilenamesNumericalIncrements(f1,f2):
    l1=getSortableList(f1)
    l2=getSortableList(f2)
    
    if len(l1)!=len(l2):
        return False
    
    for i in range(len(l1)):
        if not l1[i]==l2[i] and (not l1[i].__class__==int or not l2[i].__class__==int):
            return False
    
    return True

def getSortedListOfNumericalEquivalentFiles(f,d):
    sortedList=getSortedListOfFiles(d,globArg='*')
    for i in range(len(sortedList)-1,-1,-1):
        if not AreFilenamesNumericalIncrements(f,sortedList[i]):
            del sortedList[i]
    
    return sortedList

def getSortedListOfFilesOld(d,globArg='*[!.txt]'):# old attmpt at this using re... way too complicated...
    import re
    files = glob.glob(os.path.join(d,globArg))
    l=copy.deepcopy(files)

    start,end=0,-1
    # Find the first character that does not match in all strings
    done=False
    for i in range(len(l[0])):
        for f in l:
            if f[i]!=l[0][i]:
                start = i
                done=True
                break
        if done:
            break
    
    # Do the same from the end...
    done=False
    for i in range(len(l[0])-1,-1,-1):
        for f in l:
            if f[i]!=l[0][i]:
                end = i-len(l[0])
                done=True
                break
        if done:
            break
    
    # Find any non-numerical parts of the filename
    m=re.findall('\\D*',l[0][start:end])
    if m!=[]:
        print m
        for i in m:
            l[0].find()
    
    for i in range(len(l[0])-1,-1,-1):
        pass
    if not l[0].isdigit():
        pass
    
    l.sort(cmpGen(lambda x: float(x[start:end])))
    
    return l

# These were stolen from VolumeGif...

# Comparison to avoid the problem of improper file name sorting.
# This was not sorting properly if files are like: GBR_0_0etc, GBR_1_0etc, ..., GBR_13_0etc, ...
# Because '_' has a lower priority than digits, so 13 comes in before 1!!!
# ...Could potentially also pad the z-values w/zeros
def cmp_fnames_A(f1,f2):
    if f1==f2:
        return 0
    s1=os.path.split(f1)[1].split('_')
    s2=os.path.split(f2)[1].split('_')
    maxDigits=max(len(s1[1]),len(s2[1]))
    s1a='0'*(maxDigits-len(s1[1]))+s1[1]+'_'+s1[2]
    s2a='0'*(maxDigits-len(s2[1]))+s2[1]+'_'+s2[2]
    return (s1a>s2a)*2-1

def cmp_fnames(f1,f2):
    if f1==f2:
        return 0
    return (os.path.getmtime(f1)>os.path.getmtime(f2))*2-1
