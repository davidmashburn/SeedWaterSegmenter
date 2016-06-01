from distutils.core import setup

# Read the version number
with open("SeedWaterSegmenter/_version.py") as f:
    exec(f.read())

setup(
    name='SeedWaterSegmenter',
    version=__version__, # use the same version that's in _version.py
    author='David N. Mashburn',
    author_email='david.n.mashburn@gmail.com',
    packages=['SeedWaterSegmenter'],
    package_dir={'SeedWaterSegmenter': 'SeedWaterSegmenter'},
    package_data={'SeedWaterSegmenter': ['icons/*']},
    scripts = ['create_sws_shortcut.py'],
    url='http://pypi.python.org/pypi/SeedWaterSegmenter/',
    license='LICENSE.txt',
    description='graphical program to interactively segment image stacks of cells in tissue with edge-labels (aka. white outlines)',
    long_description=open('README.rst').read(),
    install_requires=[
# Adding install requirements may make this quirk out on a new system
# If it does, I'll have to patch it!
                      #'wxPython>=2.8', # wxPython isn't being found correctly by setuptools -- please install it manually!
                      'numpy>=1.0',
                      'scipy>=0.8',
                      'matplotlib>=1.0',
                      'pillowfight', # Depend on either PIL or pillow if available (PIL must be at least version 1.1.5)
                      'xlrd>=0.7',
                      'xlwt>=0.7',
                      'mahotas>=0.5',
#                      # Projects that used to be internal
                      'EllipseFitter>=0.1',
                      'FilenameSort>=0.1.3',
                      'GifTiffLoader>=0.2.3',
                      'ImageContour>=0.1',
                      'np_utils>=0.2.1',
                      ],
)
