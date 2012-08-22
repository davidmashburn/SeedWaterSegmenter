from distutils.core import setup

setup(
    name='SeedWaterSegmenter',
    version='0.5.3.5',
    author='David N. Mashburn',
    author_email='david.n.mashburn@gmail.com',
    packages=['SeedWaterSegmenter'],
    package_dir={'SeedWaterSegmenter': 'SeedWaterSegmenter'},
    package_data={'SeedWaterSegmenter': ['icons/*']},   
    scripts = ['postinstall.py'],
    url='http://pypi.python.org/pypi/SeedWaterSegmenter/',
    license='LICENSE.txt',
    description='',
    long_description=open('README.rst').read(),
    install_requires=[
# Adding install requirements may make this quirk out on a new system
# If it does, I'll have to patch it!
                      'wxPython>=2.8',
                      'numpy>=1.0',
                      'scipy>=0.8',
                      'matplotlib>=1.0',
                      'PIL>=1.1.5',
                      'xlrd>=0.7',
                      'xlwt>=0.7',
                      'mahotas>=0.5',
#                      # Projects that used to be internal
                      'cmpGen>=0.1',
                      'EllipseFitter>=0.1',
                      'FilenameSort>=0.1',
                      'GifTiffLoader>=0.1',
                      'ImageContour>=0.1',
                      ],
)
