# SeedWater Segmenter

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

Source code is mirrored to three repositories and to PyPI:

  * GitHub:      http://github.com/davidmashburn/SeedWaterSegmenter/
  * Bitbucket:   http://bitbucket.org/davidmashburn/seedwatersegmenter
  * Gitlab:      http://gitlab.com/davidmashburn/seedwatersegmenter
  * PyPI:        http://pypi.python.org/pypi/SeedWaterSegmenter


You may also want to read the manual, but be aware that it needs updating: "SeedWaterSegmenter V x.x Manual.txt"


## Installing and Running
----

The quickest way to install SWS is using Anaconda. Once you install Anaconda from here:

https://www.anaconda.com/products/individual

open a Terminal window ([see here for help finding your terminal](https://docs.anaconda.com/anaconda/user-guide/getting-started/#cli-hello)).

and run this command:

```
pip install SeedWaterSegmenter
```

Run SWS with this command on Windows or Linus:

```
python -c "from SeedWaterSegmenter import start_sws; start_sws()"
```

and this command on macOS:

```
pythonw -c "from SeedWaterSegmenter import start_sws; start_sws()"
```

SWS now supports Python 3. If you have issues, please make open issue on Github [here](https://github.com/davidmashburn/SeedWaterSegmenter/issues) . For Python 2 instructions on macOS, see [here](old-python2-install-instructions-mac.md).


## Making a desktop launcher with icon

While I recommend running SWS from the terminal to see the console logs, you can also create a desktop launcher to make it more convenient to launch.

### Windows:
    
Run one of these two command from the Anaconda Command Prompt (depending on whether you did a single-user or system-wide install):

```
python C:\Users\<your username>\Anaconda\Scripts\create_sws_shortcut.py -install
python C:\Anaconda\Scripts\create_sws_shortcut.py -install
```

### Mac OS X
    (CURRENTLY BROKEN)

    Now also, thanks to Sveinbjorn Thordarson's Platypus tool, a packaged app is available for download at:
    https://github.com/davidmashburn/SeedWaterSegmenter/blob/master/MacOSX/SeedWaterSegmenterApp.zip
    Just extract the zip file and place the App on the Desktop or in the Applications folder
    
    Be aware that this is only a link to the python scripts and will not work by itself without the above installation.
    
    There is also a ".command" file that can serve the same purpose if the App does not work:
    https://github.com/davidmashburn/SeedWaterSegmenter/blob/master/MacOSX/SeedWaterSegementer.command


That's it!

### Linux

Look at this to get you started:

https://github.com/davidmashburn/SeedWaterSegmenter/blob/master/desktop/SeedWaterSegmenter.desktop

This is how I created the symlinks that make this work:

```
ln -s /usr/local/lib/python2.7/dist-packages/SeedWaterSegmenter*/seedwatersegmenter/SeedWaterSegmenter.py /usr/local/bin/seedwatersegmenter
ln -s /usr/local/lib/python2.7/dist-packages/SeedWaterSegmenter*/seedwatersegmenter/icons/SeedWaterSegmenter.svg /usr/local/share/pixmaps/SeedWaterSegmenter.svg
```

That's it!

----
= Screenshot =

https://github.com/davidmashburn/SeedWaterSegmenter/blob/master/doc/SWS_Screenshot.png

