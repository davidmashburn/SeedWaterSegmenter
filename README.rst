SeedWater Segmenter
===================

Seedwater Segmenter (SWS) is a graphical Python program to interactively segment image stacks of cells in tissue with edge-labels (aka. white outlines). The interactions are entirely based on the editing of seeds, which in turn are expanded by a watershed algorithm. The major difference between SWS and other tools is that you can place more than one seed per cell which can help you adjust the boundaries of difficult cells.

SWS is built on top of wxPython, matplotlib, numpy, scipy, and mahotas. It uses Cython and Cpyx for time-critical sections of code (namely updating the seed array from the seed list).

At its core, it uses a lightning-fast watershed algorithm (thanks to the mahotas project) and allows real-time updates. It has a simple (if cluttered) UI and is fully interactive, even including 1-level undo.

There is an upcoming publication about SWS that gives all these details and more in Cytometry Part A.

This requires manual compilation of the C extension with Cpyx.py (a hack-ish Cython/gcc compiler).
I realize there are better ways to do this, but this one always works.

Source is hosted on GitHub: https://github.com/davidmashburn/SeedWaterSegmenter/

The Google Code page has binary releases (http://code.google.com/p/seedwater/) that should work on 32-bit Windows or 64-bit Linux.

For more details, see Python Install Instructions 20XX.txt
You may also want to read the manual: SeedWater Segmenter V x.x Manual.txt

----

Installing and Running
----------------------
Installing is FINALLY easy!

In fact, you can now "easwy_install SeedWaterSegmenter"

Quick instructions:
^^^^^^^^^^^^^^^^^^^
Windows/Mac:
    Install Enthought Python Distribution
    http://www.enthought.com/products/edudownload.php
Ubuntu/Debian:
    sudo apt-get install python-setuptools python-wxtools python-numpy python-scipy python-matplotlib python-imaging python-xlrd python-xlwt

In a cmd window (must run as Administrator under Windows Vista/7) or Terminal (Mac/Linux), run:
    easy_install mahotas cmpGen EllipseFitter FilenameSort GifTiffLoader ImageContour SeedWaterSegmenter

Detailed Instructions:
^^^^^^^^^^^^^^^^^^^^^^
The latest version is called:

"Python Install Instructions August 2012.txt"

https://github.com/davidmashburn/SeedWaterSegmenter/blob/master/Python%20Install%20Instructions%20August%202012.txt

----

Screenshots
-----------

.. image:: http://seedwater.googlecode.com/svn/SeedwaterScreenshot.png

(This summary is still a work in progress...)

...
===
