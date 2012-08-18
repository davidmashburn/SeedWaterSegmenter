from distutils.core import setup

setup(
    name='SeedWaterSegmenter',
    version='0.5.0',
    author='David N. Mashburn',
    author_email='david.n.mashburn@gmail.com',
    packages=['SeedWaterSegmenter'],
    scripts=[],
    url='http://pypi.python.org/pypi/SeedWaterSegmenter/',
    license='LICENSE.txt',
    description='',
    long_description=open('README.rst').read(),
    install_requires=[
# Adding install requirements may make this quirk out on a new system
# If it does, I'll have to patch it!
                      'wx>=2.8',
                      'numpy>=1.0',
                      'scipy>=0.8',
                      'matplotlib>=1.0',
                      'Image>=1.1.5',
                      'mahotas>=0.5',
#                      # For now, these are still housed internally...
#                      'GifTiffLoader>=0.1',
#                      'cmpGen>=0.1',
#                      'ImageContour>=0.1',
#                      'EllipseFitter>=0.1',
                      ],
)
