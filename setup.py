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
                      'wx>=2.8',
                      'numpy>=0.9',
                      ],
)
