# This file does not get installed and is merely a convenience so that everything acts the same if the whole directory is placed anywhere in sys.path

# This makes it so that import <name>.<module> works
# This is unnecessary for single-file packages
import os
__path__[-1] = os.path.join(__path__[-1],__name__)
del os

# This loads all the modules that should be in the top level and should be identical to the contents of __init__.py in the __name__ directory below this one...

from SeedWaterSegmenter import *
import SWHelpers
