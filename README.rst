SeedWater Segmenter
===================

Seedwater Segmenter (SWS) is a graphical Python program to interactively segment
image stacks of cells in tissue with edge-labels (aka. white outlines).
The interactions are entirely based on the editing of seeds,
which in turn are expanded by a watershed algorithm.
The major difference between SWS and other tools is that you can place more than one seed per cell,
which can help you adjust the boundaries of difficult cells.

SWS is built on top of wxPython, matplotlib, numpy, scipy, PIL, and mahotas.

At its core, it uses a lightning-fast watershed algorithm (thanks to the mahotas project) and allows real-time updates.
It has a simple (if cluttered) UI and is fully interactive, even including 1-level undo.

The publication about SWS that gives all these details and more in Cytometry Part A:

http://onlinelibrary.wiley.com/doi/10.1002/cyto.a.22034/abstract

Source code is mirrored to four repositories and to PyPI:

- GitHub:      http://github.com/davidmashburn/SeedWaterSegmenter/

- Bitbucket:   http://bitbucket.org/davidmashburn/seedwatersegmenter

- Gitorious:   http://gitorious.org/seedwatersegmenter

- Google Code: http://code.google.com/p/seedwater/

- PyPI:        http://pypi.python.org/pypi/SeedWaterSegmenter


You may also want to read the manual, but be aware that it needs updating: "SeedWaterSegmenter V x.x Manual.txt"

----

Installing and Running
----------------------
SeedWater is now pure-python. In fact, with just python and pip installed, you could install it with just::
    
    pip SeedWaterSegmenter

The big caveat to this is that SeedWater also depends on a number of binary/compiled dependencies which must be installed separately, the trickiest usually being mahotas.
It also depends on the standard scientific python packages (like numpy, scipy, matplotlib, PIL, etc), so those must also be installed.
Pretty much any method of obtaining those should work.

In short, you will need to install all of these dependencies either with binary installers, package managers, or by compiling from source.

Many great Python Distributions exist for installing Python and most of these dependencies all at once.
My current favorite is Anaconda which is free and works on both Windows and Mac.

One other important point is that SeedWater does NOT run on Python 3.x at this time, because not all of these dependencies have been ported over.
For the near future, please use Python 2.7 only.

Windows:
^^^^^^^^
Download and install Anaconda Python Distribution:
    http://continuum.io/downloads
    
    Download either 32 or 64 bit (to match your version of Windows).
    
    Double-click on the downloaded file and install it.
    I *highly* recommend doing a single-user (non-system-wide) install to save headaches with permissions.

Install wxPython and SeedWater Segmenter:
    Now, from the Start Menu (aka, the Windows Button) select the Anaconda Command Prompt inside the Anaconda folder.
    Here, run these two commands (and you may need to disable your antivirus before proceeding)::
        
        conda install wxpython
        pip install SeedWaterSegmenter

Run SeedWater:
    Now, from this same Anaconda command prompt, you can run SeedWater with the following command::
        
        python -c "from SeedWaterSegmenter import start_sws; start_sws()"

Make a desktop launcher (with an icon):
    To make running SeedWater easier, you can install a desktop shortcut with the included script.
    Just run one of these two command from the Anaconda Command Prompt
    (depending on whether you did a single-user or system-wide install)::
        
        python C:\Users\<your username>\Anaconda\Scripts\create_sws_shortcut.py -install
        python C:\Anaconda\Scripts\create_sws_shortcut.py -install

    That's it!

Mac OS X:
^^^^^^^^^
Obtain a C compiler:
    Download XCode from the Mac App Store or from https://developer.apple.com/xcode/
    
    Install it and run it.
    
    To get gcc, you must install command line tools, a package for XCode.
    You can access this from: XCode menu > Preferences > Downloads.
    Check "command line tools" and install.
    
    Reboot your system to make sure everything is loaded.

Download and install Anaconda Python Distribution:    
    http://continuum.io/downloads
    
    Download either the bash installer or the GUI installer.
    To run the bash installer, download the file to the Downloads folder and run this command in the Terminal::
        
        cd ~/Downloads
        bash Anaconda*.sh

    If you used the "pkg" download, just double-click it to install.
    
    In either case, you can then choose either a single-user of system-wide install.

Install wxPython:
    Now, open a NEW Terminal window, and run this command to install wxPython::
        
        conda install wxpython

    This failed for me on Mac OS X Snow Leopard, so I had to download this file
    (http://repo.continuum.io/pkgs/free/osx-64/wxpython-3.0-py27_0.tar.bz2)
    to the Downloads folder, and then run these commands::
        
        cd ~/Downloads
        conda install wxpython-3.0-py27_0.tar.bz2

Install SeedWater Segmenter:
    In the Terminal, run the following::
        
        pip install SeedWaterSegmenter

Run SeedWater:
    Now you can run SeedWater with the following command, noting that you HAVE to use "pythonw" and not just "python"::
        
        pythonw -c "from SeedWaterSegmenter import start_sws; start_sws()"

Download the App:
    Now also, thanks to Sveinbjorn Thordarson's Platypus tool, a packaged app is available for download:
    https://github.com/davidmashburn/SeedWaterSegmenter/blob/master/MacOSX/SeedWaterSegmenterApp.zip
    
    Just extract the zip file and place the App on the Desktop or in the Applications folder.

    Be aware that this is only a link to the python scripts and will not work by itself without the above installation.

    (There is also a ".command" file that can serve the same purpose if the App does not work at
    https://github.com/davidmashburn/SeedWaterSegmenter/blob/master/MacOSX/SeedWaterSegementer.command )

    That's it!


Ubuntu/Debian:
^^^^^^^^^^^^^^
Install:
    Run these two commands in the terminal::
        
        sudo apt-get install python-setuptools python-wxtools python-numpy python-scipy python-matplotlib python-imaging python-xlrd python-xlwt
        sudo easy_install -U SeedWaterSegmenter
    
    Run SeedWater:
    In the terminal, run::
        
        python2.7 -c "from SeedWaterSegmenter import start_sws; start_sws()"
    
    (just "python" may also work, depending on your system)
    
    Make a desktop launcher:
    Look at this to get you started:
    
    https://github.com/davidmashburn/SeedWaterSegmenter/blob/master/desktop/SeedWaterSegmenter.desktop

    This is how I created the symlinks that make this work::
        
        sudo ln -s /usr/local/lib/python2.7/dist-packages/SeedWaterSegmenter*/seedwatersegmenter/scripts/start_seedwatersegmenter.py /usr/local/bin/seedwatersegmenter
        sudo chmod +x /usr/local/bin/seedwatersegmenter
        sudo ln -s /usr/local/lib/python2.7/dist-packages/SeedWaterSegmenter*/seedwatersegmenter/icons/SeedWaterSegmenter.svg /usr/local/share/pixmaps/SeedWaterSegmenter.svg
    
    That's it!

----

Screenshot
-----------

.. image:: http://seedwater.googlecode.com/svn/SeedwaterScreenshot.png

.. image:: https://github.com/davidmashburn/SeedWaterSegmenter/blob/master/doc/SWS_Screenshot.png

