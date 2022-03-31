## Obtain a C compiler:
Download XCode from the Mac App Store or from https://developer.apple.com/xcode/

Command line tools are now installed automatically with XCode.

Reboot your system to make sure everything is loaded.

## Download and install Anaconda Python Distribution:
http://continuum.io/downloads

Download either the bash installer or the GUI installer. To run the bash installer, download the file to the Downloads folder and run this command in the Terminal:

```
cd ~/Downloads
bash Anaconda*.sh
```

If you used the "pkg" download, just double-click it to install.

In either case, you can then choose either a single-user of system-wide install.

Anaconda will be installed with a python 3.x version, but SeedWaterSegmenter only runs properly under python 2.7, so you will need to create and activate a python 2.7 environment.

After installing Anaconda, open a terminal window and create a python 2.7 environment (named py2) by typing in:

```
conda create --name py2 python=2.7
```

Then activate the py2 environment by typing:

```
conda activate py2
```

## Install wxPython and mahotas:
Now, in the same Terminal window, run these commands to install wxPython and mahotas:

```
conda install wxpython
conda install -c conda-forge mahotas
```

## Install SeedWater Segmenter:

In the Terminal, run the following:

```
pip install SeedWaterSegmenter
```

## Run SeedWater:
Now you can run SeedWater with the following command, noting that you HAVE to use `pythonw` and not just `python`:

```
pythonw -c "from SeedWaterSegmenter import start_sws; start_sws()"
```